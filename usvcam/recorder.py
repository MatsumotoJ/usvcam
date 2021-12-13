import numpy as np
import math
import time
import cv2
import copy

import nidaqmx
import nidaqmx.stream_readers
import nidaqmx.constants

import pyrealsense2 as rs

import scipy.io
import scipy.signal as signal
from scipy.ndimage import median_filter

import h5py
import sys
import datetime
import os
import yaml

import usvcam.tool as tool

config_path = sys.prefix + '/etc/usvcam/config.yaml'

def rs_start(color_mode, laser_power):

    pipe = rs.pipeline()

    config = rs.config()
    config.enable_stream(rs.stream.depth, 640, 480, rs.format.z16, 30)      # depth
    if color_mode == 0:
        config.enable_stream(rs.stream.infrared, 640, 480, rs.format.y8, 30)      # IR
    elif color_mode == 1:
        config.enable_stream(rs.stream.color, 640, 480, rs.format.rgb8, 30)    # RGB
    profile = pipe.start(config)
    depth_sensor = profile.get_device().first_depth_sensor()
    depth_sensor.set_option(rs.option.visual_preset, 5) # short range
    if laser_power is not None:
        depth_sensor.set_option(rs.option.laser_power, laser_power)

    ### get camera parameters
    frames = pipe.wait_for_frames()    
    depth_frame = frames.get_depth_frame()
    if color_mode == 0:
        color_frame = frames.get_infrared_frame()
    elif color_mode == 1:
        color_frame = frames.get_color_frame()

    depth_intrin = depth_frame.profile.as_video_stream_profile().intrinsics
    color_intrin = color_frame.profile.as_video_stream_profile().intrinsics
    depth_to_color_extrin = depth_frame.profile.get_extrinsics_to(color_frame.profile)
    ######

    return pipe, depth_intrin, color_intrin, depth_to_color_extrin

def update_spec(x, disp_spec, disp_spec_cur, daq_fs):

    f,t,s = signal.spectrogram(x, daq_fs, nperseg=256)
    s = np.log10(s) + 15
    s = (15 - s) / 15
    s = np.flip(s, axis=0)

    #s = cv2.resize(s, (int(s.shape[1] /2000 * disp_spec.shape[1]), disp_spec.shape[0]))
    s = cv2.resize(s, (max(int(s.shape[1] /2000 * disp_spec.shape[1]),1), disp_spec.shape[0]))

    n = disp_spec.shape[1] - 1 - disp_spec_cur

    if s.shape[1] > n:
        disp_spec[:,disp_spec_cur:-1] = s[:,0:n]
        disp_spec_cur = 0
    else:
        disp_spec[:,disp_spec_cur:disp_spec_cur+s.shape[1]] = s 
        disp_spec_cur += s.shape[1]

    return disp_spec, disp_spec_cur

def start_recording(data_dir, color_mode, record_depth, depth_intrin, color_intrin, depth_to_color_extrin, daq_fs, daq_n_ch, camera_height):
    
    # open output snd file
    
    fp_daq = open(data_dir + '/snd.dat', 'wb')
    global daq_total_sample
    daq_total_sample = 0

    # open output video file
    vid_size = (color_intrin.width, color_intrin.height)
    fmt = cv2.VideoWriter_fourcc('m', 'p', '4', 'v') 
    if color_mode == 0:
        writer = cv2.VideoWriter(data_dir + '/vid.mp4', fmt, 30, vid_size, isColor=False) 
    elif color_mode == 1:
        writer = cv2.VideoWriter(data_dir + '/vid.mp4', fmt, 30, vid_size, isColor=True) 

    # sync between camera and DAQ
    fp_sync = open(data_dir + '/sync.csv', 'w')

    # save recording parameters
    with open(config_path, 'r') as f:
        usvcam_cfg = yaml.safe_load(f)
    speedOfSound = float(usvcam_cfg['speed_of_sound'])
    pressure_calib = np.array(usvcam_cfg['sound_pressure_calibration'])
    mic0pos = np.asarray(usvcam_cfg['mic0pos'])/1000.0

    with h5py.File(data_dir + '/param.h5', mode='w') as f:
        tool.save_intrinsics(f, '/camera_param/depth_intrin', depth_intrin)
        tool.save_intrinsics(f, '/camera_param/color_intrin', color_intrin)
        tool.save_extrinsics(f, '/camera_param/depth_to_color_extrin', depth_to_color_extrin)
        f.create_dataset('/camera_param/color_mode', data = color_mode)
        f.create_dataset('/camera_param/camera_height', data = camera_height)
        f.create_dataset('/daq_param/fs', data = daq_fs)
        f.create_dataset('/daq_param/n_ch', data = daq_n_ch)
        f.create_dataset('/misc/speedOfSound', data=speedOfSound)
        f.create_dataset('/misc/pressure_calib', data=pressure_calib)
        f.create_dataset('/misc/mic0pos', data=mic0pos)
    
    if record_depth:
        fp_depth = h5py.File(data_dir + '/depth.h5', mode='w')
    else:
        fp_depth = -1

    return fp_daq, writer, fp_sync, fp_depth

is_recording = False
ch_preview = 0
daq_total_sample = 0
daq_buf_spec_n = 0
daq_n_ch = 4
daq_vmax = 1.0
daq_buf_spec_max = 10000
daq_buf_spec = np.zeros(daq_buf_spec_max)

def main():

    global is_recording, ch_preview, daq_total_sample, daq_buf_spec_n
    global daq_n_ch, daq_vmax, daq_buf_spec_max, daq_buf_spec

    with open(config_path, 'r') as f:
        usvcam_cfg = yaml.safe_load(f)

    color_mode = usvcam_cfg['color_mode']
    ch_preview = usvcam_cfg['ch_monitor']
    laser_power = usvcam_cfg['laser_power']
    record_depth = usvcam_cfg['record_depth']
    devname = usvcam_cfg['daq_dev_name']

    is_recording = False

    if len(sys.argv) > 1:
        outdir_prefix = sys.argv[1]
    else:
        outdir_prefix = ''

    ##### (1) initialize analog data acquisition
    daq_fs = usvcam_cfg['daq_fs']              
    daq_n_ch = 4
    daq_vmax = 1.0
    daq_nsample = int(daq_fs/1000)     # 1 msec
    daq_dev_list = devname + '/ai0, ' + devname + '/ai1, ' + devname + '/ai2, ' + devname + '/ai3'

    daq_buf_spec_max = daq_fs
    daq_buf_spec = np.zeros(daq_buf_spec_max)
    daq_buf_spec_n = 0

    daq_task = nidaqmx.Task()
    daq_task.ai_channels.add_ai_voltage_chan(daq_dev_list, min_val=-daq_vmax, max_val=daq_vmax)
    daq_task.timing.cfg_samp_clk_timing(rate = daq_fs,
                                        sample_mode = nidaqmx.constants.AcquisitionType.CONTINUOUS,
                                        active_edge = nidaqmx.constants.Edge.RISING,
                                        samps_per_chan = daq_nsample)
    daq_reader = nidaqmx.stream_readers.AnalogMultiChannelReader(daq_task.in_stream)
    daq_total_sample = 0
    daq_running = False

    def daq_callback(task_handle, every_n_samples_event_type, number_of_samples, callback_data):
            
        if not daq_running: 
            return 0

        buf = np.zeros([4, number_of_samples])
        daq_reader.read_many_sample(data=buf, number_of_samples_per_channel = number_of_samples)

        global daq_total_sample, is_recording
        if is_recording:
            r = 30000 / daq_vmax
            a = (buf.T * r).astype(np.int16)
            fp_daq.write(a.tobytes())
            daq_total_sample += number_of_samples

        global daq_buf_spec_n
        if daq_buf_spec_n < daq_buf_spec_max:
            daq_buf_spec[daq_buf_spec_n:(daq_buf_spec_n+number_of_samples)] = buf[ch_preview,:]
            daq_buf_spec_n += number_of_samples

        return 0

    daq_task.register_every_n_samples_acquired_into_buffer_event(daq_nsample, daq_callback)
    ##### end of (1)

    # initialize camera
    pipe, depth_intrin, color_intrin, depth_to_color_extrin = rs_start(color_mode, laser_power)
    vid_size = (color_intrin.width, color_intrin.height)

    # initialize spectrogram display
    disp_spec = np.zeros([200, vid_size[0]])
    disp_spec_cur = 0

    # daq start
    daq_task.start()
    daq_running = True

    # open monitor windows
    cv2.namedWindow('Camera', cv2.WINDOW_AUTOSIZE)
    cv2.namedWindow('Microphone', cv2.WINDOW_AUTOSIZE)
    cv2.moveWindow('Camera', 80, 10)
    cv2.moveWindow('Microphone', 80, 550)

    # initialize fps count
    fps = np.zeros(30)
    frames = pipe.wait_for_frames()
    t_pre = frames.get_timestamp()

    while True:
        # get new frame
        frames = pipe.wait_for_frames()

        # calc & show spectrogram
        if daq_buf_spec_n > 0:
            x = daq_buf_spec[0:daq_buf_spec_n]
            daq_buf_spec_n = 0

            disp_spec, disp_spec_cur = update_spec(x, disp_spec, disp_spec_cur, daq_fs)

            disp_spec2 = copy.deepcopy(disp_spec)
            disp_spec2[:, disp_spec_cur] = 0

            cv2.imshow('Microphone', disp_spec2)

        # update fps
        t = frames.get_timestamp()
        fps = np.roll(fps, -1)
        fps[-1] = t-t_pre
        t_pre = t

        # get color frame
        if color_mode == 0:
            color_frame = frames.get_infrared_frame()
        elif color_mode == 1:
            color_frame = frames.get_color_frame()
        color_image = np.asanyarray(color_frame.get_data())

        if color_mode == 1:
            color_image = cv2.cvtColor(color_image, cv2.COLOR_BGR2RGB)

        # record video
        if is_recording:
            fp_sync.write('{:.2f}, {:.0f}\n'.format(frames.get_timestamp(), daq_total_sample))
            writer.write(color_image)
            if record_depth:
                depth_frame = frames.get_depth_frame()
                depth_image = np.asanyarray(depth_frame.get_data())
                fp_depth.create_dataset('/f-{:06}'.format(framecnt), data=depth_image, compression='gzip')
            framecnt += 1

        # show camera image
        if color_mode == 0:
            disp_camera = cv2.cvtColor(color_image,cv2.COLOR_GRAY2BGR)
        else:
            disp_camera = copy.deepcopy(color_image)

        cv2.putText(disp_camera, '{:.1f} fps'.format(1000/np.mean(fps)) ,(10,30), cv2.FONT_HERSHEY_SIMPLEX, 1, (255,255,255), 2, cv2.LINE_AA)
        if is_recording:
            cv2.putText(disp_camera, '[R] {:.1f} sec'.format(daq_total_sample/daq_fs) ,(10,70), cv2.FONT_HERSHEY_SIMPLEX, 1, (50,50,255), 2, cv2.LINE_AA)
        cv2.imshow('Camera', disp_camera)

        # key inputs
        k = cv2.waitKey(1)
        if k== ord('s') and not is_recording:
            if len(outdir_prefix) > 0:
                outdir = './data/' + outdir_prefix + '_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            else:
                outdir = './data/' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
                
            os.makedirs(outdir, exist_ok=True)
            fp_daq, writer, fp_sync, fp_depth = start_recording(outdir, color_mode, record_depth, depth_intrin, color_intrin, depth_to_color_extrin, daq_fs, daq_n_ch, usvcam_cfg['camera_height'])
            is_recording = True
            framecnt = 0
        if k== 27:
            break

    # stop recording
    if is_recording:
        is_recording = False
        fp_sync.close()
        fp_daq.close()
        writer.release()
        if record_depth:
            fp_depth.close()
        
    # stop devices
    daq_running = False
    daq_task.stop()

    pipe.stop()
    daq_task.close()

if __name__ == '__main__':
    main()
