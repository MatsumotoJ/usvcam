import numpy as np
import cv2
import copy
import threading
import scipy.signal as signal
from scipy.ndimage import median_filter

import h5py
import sys
import datetime
import os
import yaml

import subprocess

import pyrealsense2 as rs
from .fadc import *
from . import tool
script_dir = os.path.dirname(__file__)
config_path = script_dir + '/config_fadc.yaml' 

cam_delay = 0.1

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
    #depth_sensor.set_option(rs.option.visual_preset, 5) # short range
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
    s = cv2.resize(s, (max(int(s.shape[1] / (daq_fs/256) * disp_spec.shape[1]),1), disp_spec.shape[0]))

    n = disp_spec.shape[1] - 1 - disp_spec_cur

    if s.shape[1] > n:
        disp_spec[:,disp_spec_cur:-1] = s[:,0:n]
        disp_spec_cur = 0
    else:
        disp_spec[:,disp_spec_cur:disp_spec_cur+s.shape[1]] = s 
        disp_spec_cur += s.shape[1]

    return disp_spec, disp_spec_cur

def start_recording(data_dir, color_mode, record_depth, depth_intrin, color_intrin, depth_to_color_extrin, daq_fs, daq_n_ch, camera_height):
    
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
        f.create_dataset('/camera_param/cam_delay', data = cam_delay)
        f.create_dataset('/daq_param/fs', data = daq_fs)
        f.create_dataset('/daq_param/n_ch', data = daq_n_ch)
        f.create_dataset('/misc/speedOfSound', data=speedOfSound)
        f.create_dataset('/misc/pressure_calib', data=pressure_calib)
        f.create_dataset('/misc/mic0pos', data=mic0pos)
    
    if record_depth:
        fp_depth = h5py.File(data_dir + '/depth.h5', mode='w')
    else:
        fp_depth = -1

    # ffmpeg start for real-time encoding
    ffmpeg_cmd = [
        'ffmpeg',
        '-f', 's16le',           # Input format: 16-bit signed little-endian PCM
        '-ar', str(daq_fs),          # Sample rate
        '-ac', str(daq_n_ch),        # Audio channels
        '-i', 'pipe:0',          # Input comes from stdin (pipe)
        '-f', 'flac',            # Output format (flac)
        '-compression_level', str(5),
        '-y',
        data_dir + '/snd.flac'               # Output file
    ]
    ffmpeg_process = subprocess.Popen(ffmpeg_cmd, stdin=subprocess.PIPE)

    return ffmpeg_process, writer, fp_sync, fp_depth

is_recording = False
ch_preview = 0
daq_total_sample = 0
daq_buf_spec_n = 0
daq_n_ch = 8
daq_buf_spec_max = 10000
daq_buf_spec = np.zeros(daq_buf_spec_max)

daq_devid = -1
daq_chunksize = int(384 * 4)

def daq_proc(dev_handle, chunk_size, stop_event):
    
    ch_preview = 0

    while not stop_event.is_set():

        buf = bytearray(chunk_size * daq_n_ch * 2)
        bytes_returned = ctypes.c_uint32()
        ret = fadcRead(dev_handle, (ctypes.c_ubyte*len(buf)).from_buffer(buf), len(buf), ctypes.pointer(bytes_returned))

        if bytes_returned.value != len(buf):
            print('length error', bytes_returned.value, '/', len(buf))

        buf = np.frombuffer(buf, dtype=np.int16)
        buf = np.reshape(buf, [-1, daq_n_ch])
        number_of_samples = buf.shape[0]

        global daq_total_sample, is_recording, ffmpeg_process
        if is_recording:
            #fp_daq.write(buf.tobytes())
            ffmpeg_process.stdin.write(buf.tobytes())
            daq_total_sample += number_of_samples

        global daq_buf_spec_n, daq_buf_spec
        if daq_buf_spec_n+number_of_samples < daq_buf_spec_max:
            daq_buf_spec[daq_buf_spec_n:(daq_buf_spec_n+number_of_samples)] = buf[:,ch_preview].astype(float)/32000.0
            daq_buf_spec_n += number_of_samples

    return 0

def main():

    global is_recording, ch_preview, daq_total_sample, daq_buf_spec_n
    global daq_n_ch, daq_vmax, daq_buf_spec_max, daq_buf_spec, ffmpeg_process

    with open(config_path, 'r') as f:
        usvcam_cfg = yaml.safe_load(f)

    color_mode = usvcam_cfg['color_mode']
    ch_preview = usvcam_cfg['ch_monitor']
    laser_power = usvcam_cfg['laser_power']
    record_depth = usvcam_cfg['record_depth']
    #devname = usvcam_cfg['daq_dev_name']

    is_recording = False

    if len(sys.argv) > 1:
        outdir_prefix = sys.argv[1]
    else:
        outdir_prefix = ''

    ##### (1) initialize sound acquisition
    #daq_fs = int(usvcam_cfg['daq_fs'])
    daq_fs = int(384000)    # fixed

    # open device
    dev_handle = ctypes.c_void_p()
    ret = fadcOpenByIndex(ctypes.pointer(dev_handle),0)

    # start ADC
    ret = fadcStart(dev_handle,1,0) # 1 = 384kHz

    daq_total_sample = 0

    daq_stopevent = threading.Event()
    daq_thread = threading.Thread(target=daq_proc, args=(
                dev_handle, daq_chunksize, daq_stopevent))
    
    ##### end of (1)

    # initialize camera
    pipe, depth_intrin, color_intrin, depth_to_color_extrin = rs_start(color_mode, laser_power)
    vid_size = (color_intrin.width, color_intrin.height)

    # initialize spectrogram display
    disp_spec = np.zeros([200, vid_size[0]])
    disp_spec_cur = 0

    # daq start
    daq_thread.start()

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
            ffmpeg_process, writer, fp_sync, fp_depth = start_recording(outdir, color_mode, record_depth, depth_intrin, color_intrin, depth_to_color_extrin, daq_fs, daq_n_ch, usvcam_cfg['camera_height'])
            is_recording = True
            framecnt = 0
        if k== 27:
            break

    # stop recording
    if is_recording:
        is_recording = False
        fp_sync.close()
        ffmpeg_process.stdin.close()
        ffmpeg_process.wait()
        writer.release()
        if record_depth:
            fp_depth.close()
        
    daq_stopevent.set()
    daq_thread.join()

    ret = fadcStop(dev_handle)
    ret = fadcClose(dev_handle)

    pipe.stop()

if __name__ == '__main__':
    main()
