# legacy app: mic pos calibration with a 40kHz ultrasonic distance sensor

import numpy as np
import math
import time
import cv2
import copy

import nidaqmx
import nidaqmx.stream_readers
import nidaqmx.constants

import pyrealsense2 as rs

from cv2 import aruco

import h5py

import platform
if platform.system() == "Windows":
    import winsound

import scipy.signal
import scipy.optimize

import datetime
import yaml

import usvcam.tool as tool

import sys
config_path = sys.prefix + '/etc/usvcam/config.yaml'

def rs_start():

    pipe = rs.pipeline()

    config = rs.config()
    config.enable_stream(rs.stream.color, 640, 480, rs.format.rgb8, 30)    # color
    config.enable_stream(rs.stream.depth, 640, 480, rs.format.z16, 30)      # depth
    profile = pipe.start(config)
    depth_sensor = profile.get_device().first_depth_sensor()
    depth_sensor.set_option(rs.option.visual_preset, 5) # short range

    pc = rs.pointcloud()

    ### get camera parameters
    frames = pipe.wait_for_frames()    
    depth_frame = frames.get_depth_frame()
    color_frame = frames.get_color_frame()

    depth_intrin = depth_frame.profile.as_video_stream_profile().intrinsics
    color_intrin = color_frame.profile.as_video_stream_profile().intrinsics
    depth_to_color_extrin = depth_frame.profile.get_extrinsics_to(color_frame.profile)
    ######

    return pipe, pc, depth_intrin, color_intrin, depth_to_color_extrin

def update_markerpos(color_image, dict_aruco, markerpos):

    markerpos_pre = copy.deepcopy(markerpos)

    corners, ids, rejectedImgPoints = aruco.detectMarkers(color_image, dict_aruco) 

    markerpos = np.zeros([2, 2])
    markerpos[:,:] = np.nan

    if ids is not None:
        for i in range(len(ids)):
            ii = ids[i][0]
            if ii > 2:
                continue
            #cv2.fillPoly(color_image, corners[i].astype(np.int), clr[ii])
            x = np.mean(corners[i][0,:,:], axis=0)
            markerpos[ii,:] = x

    if np.sum(np.isnan(markerpos)) > 0:
        d_markerpos = float('inf')
    else:
        d_markerpos = np.max(np.linalg.norm(markerpos - markerpos_pre, axis=1))

    return markerpos, d_markerpos, corners, ids

def draw_labels(color_image, ids, corners, px_checked, i_mode):
    clr = [[[0, 0, 255], [0, 0, 200]], [[0, 255, 255], [0, 200, 200]], [[0, 255, 0], [0, 200, 0]]]
    if ids is not None:
        for i in range(len(ids)):
            ii = ids[i][0]
            if ii > 2:
                continue
            cv2.polylines(color_image, corners[i].astype(np.int), True, clr[i_mode][ii], thickness=5)

    for px in px_checked:
        cv2.circle(color_image, px, 8, [0, 255, 0], thickness=1)

def draw_snd(x, w, v_dispmax):
    img_snd = np.zeros([200, w, 3], np.uint8)
    clr = [[0, 0, 255], [255, 255, 0], [255, 0, 255], [0, 255, 255]]
    for i in range(daq_n_ch):
        xx = x[i, :]
        pts = np.vstack([np.arange(len(xx))/len(xx)*img_snd.shape[1], xx/v_dispmax*img_snd.shape[0]/2+img_snd.shape[0]/2])
        pts = pts.T
        cv2.polylines(img_snd, [pts.astype(np.int32)], False, clr[i])

    return img_snd

def get_speaker_pos(ids, corners, pc, color_image, color_frame, depth_frame):
    mask_image = np.zeros((color_image.shape[0],color_image.shape[1],1), np.uint8)
    if ids is not None:
        for i in range(len(ids)):
            c = int(ids[i][0])+1
            cv2.fillPoly(mask_image, corners[i].astype(np.int), c)
            
    pc.map_to(color_frame)
    points = pc.calculate(depth_frame)
    v, t = points.get_vertices(), points.get_texture_coordinates()
    verts = np.asanyarray(v).view(np.float32).reshape(-1, 3)  # xyz
    texcoords = np.asanyarray(t).view(np.float32).reshape(-1, 2)  # uv

    #convert unit of texcoord map (0 to 1) to pixels 
    cw, ch = color_image.shape[:2][::-1]
    v, u = (texcoords * (cw, ch) + 0.5).astype(np.uint32).T
    np.clip(u, 0, ch-1, out=u)
    np.clip(v, 0, cw-1, out=v)

    # get point cloud corresponding aruco marker
    I = mask_image[u, v, 0]
    v_marker1 = verts[I==1,:]
    v_marker2 = verts[I==2,:]
    p = (np.median(v_marker1, axis=0) + np.median(v_marker2, axis=0))/2
    p[2] -= 0.0095 # speaker = 9.5 mm height
    
    p = rs.rs2_transform_point_to_point(depth_to_color_extrin, p.tolist())
    px = rs.rs2_project_point_to_pixel(color_intrin, p)

    return p, px

dtime_str = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
out_filename = './data/'+ dtime_str + '.calib.h5'

with open(config_path, 'r') as f:
    usvcam_cfg = yaml.safe_load(f)

devname = usvcam_cfg['daq_dev_name']
mic0pos = np.asarray(usvcam_cfg['mic0pos'])/1000.0

##### (1) initialize analog data acquisition
daq_fs = 3500000                # 3.5M Hz
daq_nsample = int(daq_fs*0.001)     # 0.001 sec
daq_n_ch = 4
daq_dev_list = devname + '/ai0, ' + devname + '/ai1, ' + devname + '/ai2, ' + devname + '/ai3'
daq_vmax = 10.0

daq_buf_max = daq_fs
daq_buf = np.zeros([daq_n_ch, daq_buf_max])
daq_buf_n = 0

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

        global daq_total_sample
        
        daq_total_sample += number_of_samples

        buf = np.zeros([4, number_of_samples])
        daq_reader.read_many_sample(data=buf, number_of_samples_per_channel = number_of_samples)

        global daq_buf_n
        if daq_buf_n < daq_buf_max - number_of_samples:
            daq_buf[:,daq_buf_n:(daq_buf_n+number_of_samples)] = buf
            daq_buf_n += number_of_samples

        return 0

daq_task.register_every_n_samples_acquired_into_buffer_event(daq_nsample, daq_callback)
##### end of (1)

# initialize camera
pipe, pc, depth_intrin, color_intrin, depth_to_color_extrin = rs_start()

# save camera parameters
with h5py.File(out_filename, mode='w') as f:
    tool.save_intrinsics(f, '/camera_param/depth_intrin', depth_intrin)
    tool.save_intrinsics(f, '/camera_param/color_intrin', color_intrin)
    tool.save_extrinsics(f, '/camera_param/depth_to_color_extrin', depth_to_color_extrin)
    f.create_dataset('/daq_param/fs', data = daq_fs)
    f.create_dataset('/daq_param/n_ch', data = daq_n_ch)

# initialize speaker pos detection
d_markerpos_thr = 1.0
dict_aruco = aruco.getPredefinedDictionary(aruco.DICT_4X4_50)
markerpos = np.zeros([2, 2])
px_checked = list()
cnt_stay = 0

# start daq
daq_task.start()
daq_running = True
daq_t0 = time.time()

# open monitor windows
app_winname = 'Calibrator'
cv2.namedWindow(app_winname, cv2.WINDOW_AUTOSIZE)
cv2.moveWindow(app_winname, 80, 10)

wsize = int(2048) #int(daq_fs*0.0005)
v_thr = 0.5
v_dispmax = 3

img_snd = np.zeros([200, color_intrin.width, 3], np.uint8)

cnt = 0
i_mode = 0
while True:
    # get new frame
    frames = pipe.wait_for_frames()
    depth_frame = frames.get_depth_frame()
    color_frame = frames.get_color_frame()

    color_image = np.asanyarray(color_frame.get_data())
    color_image = cv2.cvtColor(color_image, cv2.COLOR_BGR2RGB)

    # check speaker position
    markerpos, d_markerpos, corners, ids = update_markerpos(color_image, dict_aruco, markerpos)
   
    if d_markerpos > d_markerpos_thr:
        cnt_stay = 0
    else: 
        cnt_stay += 1
    
    if cnt_stay < 15:
        i_mode = 0
    elif i_mode == 0:
        i_mode = 1
    
    if daq_buf_n > 0:

        X = copy.deepcopy(daq_buf[:, 0:daq_buf_n])
        daq_buf_n = 0

        if i_mode == 1:

            x = X[0,:]
            i_max = np.argmax(x)

            if x[i_max] > v_thr and i_max - wsize > 0 and i_max + wsize < len(x):
                
                cnt += 1

                x = X[:, i_max-wsize:i_max+wsize]
                img_snd = draw_snd(x, color_intrin.width, v_dispmax)

                p, px = get_speaker_pos(ids, corners, pc, color_image, color_frame, depth_frame)
                px_checked.append((int(px[0]), int(px[1])))

                with h5py.File(out_filename, mode='a') as f:
                    group_name = '/sig_{:05}'.format(cnt)
                    f.create_dataset(group_name + '/color_image', data = color_image)
                    f.create_dataset(group_name + '/snd', data = x)
                    f.create_dataset(group_name + '/position', data = p)

                if platform.system() == "Windows":
                    winsound.Beep(1000,100)
                i_mode = 2
                

    
    draw_labels(color_image, ids, corners, px_checked, i_mode)
    img_disp = np.concatenate((color_image, img_snd))
    cv2.imshow(app_winname, img_disp)

    # exit with ESC key
    k = cv2.waitKey(1)
    if k== 27:
        break

daq_running = False
daq_task.stop()
daq_task.close()
pipe.stop()

f_target = 40000
tool.calc_micpos(out_filename, f_target, mic0pos, False, False)
