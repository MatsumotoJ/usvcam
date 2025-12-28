from intervaltree.intervaltree import IntervalTree
import numpy as np
import numexpr as ne
import scipy.optimize
import h5py
import itertools
import math
import glob
import copy
import scipy.interpolate
import scipy.signal
import scipy.io
import scipy.io.wavfile
import sklearn.cluster
import wave
import cv2
import os
import shutil
import joblib
import yaml
import intervaltree
import sys
from tqdm import tqdm
from multiprocessing import Pool, freeze_support, RLock
import matplotlib.pyplot as plt
import soundfile
import io
import usvseg
from pathlib import Path
from datetime import datetime

from collections import deque

script_dir = os.path.dirname(__file__)
config_path = script_dir + '/config_fadc.yaml' 

cam_delay_default = 0.1
v_fps_default = 30
only_max_default = True

usegpu = False
cupy = None

def set_default_analysis_params_for_device(dev_name):
    global only_max_default
    if dev_name == 'legacy4ch':  # for device used in Matsumoto et al. 2022, iScience, 
        only_max_default = False
    elif dev_name == 'katou8ch': # for BSA-USVCAM0808 by Katou Acoustics Consultant Office
        only_max_default = True
    else:
        raise ValueError(f'Unknown device name: {dev_name}')

def enable_gpu():
    print('gpu calculation enabled')
    global usegpu, cupy
    usegpu = True
    import cupy as cp
    cupy = cp

def disable_gpu():
    print('gpu calculation disabled')
    global usegpu
    usegpu = False

def change_param_h5(data_dir, key, data):
    """
    change param h5 with backup
    """

    filename = data_dir + '/param.h5'

    src = Path(filename)
    bak_dir = src.parent / "bak"
    bak_dir.mkdir(exist_ok=True)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    # Example: data.h5 -> bak/data.20251227_153012.h5
    dst = bak_dir / f"{src.stem}.{ts}{src.suffix}"

    # Safety check in case multiple backups are created within the same second
    i = 1
    while dst.exists():
        dst = bak_dir / f"{src.stem}.{ts}_{i:02d}{src.suffix}"
        i += 1

    # Back up the original file (raises if the source file does not exist)
    shutil.copy2(src, dst)

    # Replace the dataset
    with h5py.File(src, "a") as f:
        if key in f:
            del f[key]
        f.create_dataset(key, data=data)

def open_dat(data_dir):

    if os.path.isfile(data_dir + '/snd.dat'):
        return open(data_dir + '/snd.dat', 'rb')
    elif os.path.isfile(data_dir + '/snd.flac'):
        return soundfile.SoundFile(data_dir + '/snd.flac', 'r')
    else:
        raise ValueError("Unsupported file format. Only '.dat' and '.flac' files are supported.")

def read_dat(fp_dat, t_intv, fs, n_ch):
    simg = None
    if isinstance(fp_dat, io.BufferedReader):      # dat file
        spos = (int(t_intv[0]*fs)*2*n_ch).astype(np.int64)
        fp_dat.seek(spos, 0)
        simg = np.fromfile(fp_dat, np.int16, int(n_ch*int((t_intv[1]-t_intv[0])*fs)) )
        simg = simg.reshape([-1, n_ch])

    elif isinstance(fp_dat, soundfile.SoundFile):   # flac file
        start_sample = int(t_intv[0] * fs)
        end_sample = int(t_intv[1] * fs)

        fp_dat.seek(start_sample)
        simg = fp_dat.read(end_sample - start_sample, dtype='int16')

        # Ensure the correct shape
        if simg.shape[1] != n_ch:
            simg = simg.reshape([-1, n_ch])

    return simg

def angle2point(data_dir, calibfile, azimuth, elevation):

    with h5py.File(calibfile, mode='r') as f:
        micpos = f['/result/micpos'][()]
    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        r = f['/camera_param/depth_to_color_extrin/rotation'][()]
        t = f['/camera_param/depth_to_color_extrin/translation'][()]
        ppx = f['/camera_param/color_intrin/ppx'][()]
        ppy = f['/camera_param/color_intrin/ppy'][()]
        fx = f['/camera_param/color_intrin/fx'][()]
        fy = f['/camera_param/color_intrin/fy'][()]
        coeff = f['/camera_param/color_intrin/coeffs'][()]
        im_w = f['/camera_param/color_intrin/width'][()]
        im_h = f['/camera_param/color_intrin/height'][()]
        fs = f['/daq_param/fs'][()]
        camera_height = f['/camera_param/camera_height'][()]
        #if '/camera_param/camera_height' in f.keys():
        #    camera_height = f['/camera_param/camera_height'][()]
        #else:
        #    with open(config_path, 'r') as f_cfg:
        #        usvcam_cfg = yaml.safe_load(f_cfg)
        #        camera_height = usvcam_cfg['camera_height']

    r = np.reshape(r, [3, 3]).T
    cm = np.zeros([4,4])
    cm[:3, :3] = r
    cm[:3, 3] = t
    cm[3,3] = 1

    mtx = [[fx, 0, ppx],
           [0, fy, ppy],
           [0, 0, 1]]
    mtx = np.array(mtx)
    dist = coeff[np.newaxis, :]
    tvec = cm[0:3,3]
    rvec = cv2.Rodrigues(cm[0:3,0:3])[0]

    x = np.array([camera_height]*azimuth.shape[0])
    y = np.tan(azimuth)*x
    r = np.sqrt(x**2 + y**2)
    z = np.tan(elevation)*r

    P = np.array([y,z,x]).T
    mp0 = np.reshape(micpos[0,:], [1, 3])
    P = P+mp0

    Px = np.zeros([P.shape[0], 2])
    for i in range(P.shape[0]):
        p = P[i,:]

        imgp, jac = cv2.projectPoints(p, rvec, tvec, mtx, dist)
        imgp = np.squeeze(imgp)

        u = imgp[0]
        v = imgp[1]

        Px[i,0] = u
        Px[i,1] = v

    return Px

def assign_all_segs_simple(data_dir, calibfile, assignfile, conf_thr=0.99, d_max=10):
    """
    This algorithm is considered beta.
    It is grounded in existing literature, but the present method and its
    implementation have not yet been subject to peer review.
    """

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        speedOfSound = f['/misc/speedOfSound'][()]
        pressure_calib = f['/misc/pressure_calib'][()]
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    loc = np.genfromtxt(data_dir + '/loc.csv', delimiter=',')
    if loc.ndim == 1:
        loc = np.reshape(loc, [1, -1])

    snout_pos = np.genfromtxt(data_dir + '/snout.csv', delimiter=',')
    n_mice = int(snout_pos.shape[1]/2)
    snout_pos = snout_pos[:,:n_mice*2]

    with h5py.File(assignfile, mode='r') as f:   
        a = f['/a'][()]
        b = f['/b'][()]
        fitted_model = ErrDensityModel(a, b)

    seg = load_usvsegdata_ss(data_dir)
    _, I_ss = np.unique(seg[:,3:5], axis=0, return_inverse=True)

    mask = None
    if os.path.exists(data_dir + '/mask.csv'):
        mask = np.loadtxt(data_dir + '/mask.csv', delimiter=',')
        print('mask data is found', flush=True)
        tree_mask1 = intervaltree.IntervalTree.from_tuples(mask[mask[:,2]==1,:2])
        tree_mask2 = intervaltree.IntervalTree.from_tuples(mask[mask[:,2]==2,:2])

    outdata = []

    for i_ss in tqdm(range(np.max(I_ss)+1)):

        seg2 = seg[I_ss==i_ss,:]
        i_frame = time2vidframe(data_dir, (np.min(seg2[:,0])+np.max(seg2[:,0]))/2, T, fs, cam_delay)

        od = np.zeros([seg2.shape[0], seg2.shape[1]+6])
        od[:,:] = np.nan
        od[:,6:] = seg2
        od[:,0] = i_ss

        peaks = np.reshape(loc[i_ss,3:], [-1, 2])   

        if mask is not None and len(tree_mask2.overlap(np.min(seg2[:,0]), np.max(seg2[:,0]))):
            od[:,1] = -5 # masked - not being in the time of interest

        elif mask is not None and len(tree_mask1.overlap(np.min(seg2[:,0]), np.max(seg2[:,0]))):
            od[:,1] = -4 # masked - audible call existing
        
        elif np.isnan(peaks[0,0]):
            od[:,1] = -2 # not localized 

        elif i_frame >= snout_pos.shape[0]:
            od[:,1] = -3 # no tracking data 

        elif np.sum(np.isnan(snout_pos[i_frame, :])) == 0:  # if snout pos of all rats are detected

            snouts = np.reshape(snout_pos[i_frame, :], [-1, 2])
            az, el = point2angle(data_dir, calibfile, snouts)
            snouts_a = rad2deg(np.array([az, el]).T)
            
            az, el = point2angle(data_dir, calibfile, peaks)
            peaks_a = rad2deg(np.array([az, el]).T)
            
            D = np.zeros([snouts.shape[0]], dtype=float)
            D[:] = np.inf
            for i_mouse, s_a in enumerate(snouts_a):

                d = np.linalg.norm(peaks_a - s_a, axis=1)
                d = np.nanmin(d)
                D[i_mouse] = d
            
            i_min = np.argmin(D)
            d_min = D[i_min]

            if d_min > d_max:
                od[:,1] = -6 # too far from snout
            else:
                p = fitted_model(D)
                mpi = p[i_min] / np.sum(p)
                if mpi < conf_thr:
                    od[:,1] = -1 # not assigned (n.s.)
                else:
                    od[:,1] = i_min

                od[:,3] = mpi
                od[:,4] = np.nan
                od[:,5] = np.nan

        else:
            od[:,1] = -3 # no tracking data 

        outdata.append(od)

    np.savetxt(data_dir + '/assign.csv', np.concatenate(outdata, axis=0), delimiter=',',
                header=' segment_id,1st_place_mouse,2nd_place_mouse,confidence,snout_distance,p_value,time,freq,ampitude,ID_B,ID_A')

def assign_all_segs(data_dir, calibfile, assignfile, n_mice, conf_thr=0.99):

    #with open(config_path, 'r') as f:
    #    usvcam_cfg = yaml.safe_load(f)
    #speedOfSound = float(usvcam_cfg['speed_of_sound'])
    #pressure_calib = np.array(usvcam_cfg['sound_pressure_calibration'])

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        speedOfSound = f['/misc/speedOfSound'][()]
        pressure_calib = f['/misc/pressure_calib'][()]
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    loc = np.genfromtxt(data_dir + '/loc.csv', delimiter=',')
    if loc.ndim == 1:
        loc = np.reshape(loc, [1, -1])

    snout_pos = np.genfromtxt(data_dir + '/snout.csv', delimiter=',')
    snout_pos = snout_pos[:,:n_mice*2]

    with h5py.File(assignfile, mode='r') as f:   
        clut = []
        for i in range(1,n_mice): 
            clut.append(f['/conf/N_other_{:d}/lut'.format(i)][()])
        D_clut = f['/conf/N_other_1/d'][()]
        P_clut = f['/conf/N_other_1/p'][()]
        de = f['/dist_err_dist_in_rad'][()]
        d_thr = np.percentile(de, conf_thr*100)

    seg = load_usvsegdata_ss(data_dir)
    _, I_ss = np.unique(seg[:,3:5], axis=0, return_inverse=True)

    mask = None
    if os.path.exists(data_dir + '/mask.csv'):
        mask = np.loadtxt(data_dir + '/mask.csv', delimiter=',')
        print('mask data is found', flush=True)
        tree_mask1 = intervaltree.IntervalTree.from_tuples(mask[mask[:,2]==1,:2])
        tree_mask2 = intervaltree.IntervalTree.from_tuples(mask[mask[:,2]==2,:2])

    outdata = []
    with open_dat(data_dir) as fp_dat:
        for i_ss in tqdm(range(np.max(I_ss)+1)):

            seg2 = seg[I_ss==i_ss,:]

            i_frame = time2vidframe(data_dir, (np.min(seg2[:,0])+np.max(seg2[:,0]))/2, T, fs, cam_delay)

            od = np.zeros([seg2.shape[0], seg2.shape[1]+6])
            od[:,:] = np.nan
            od[:,6:] = seg2
            od[:,0] = i_ss

            peaks = np.reshape(loc[i_ss,3:], [-1, 2])

            if mask is not None and len(tree_mask2.overlap(np.min(seg2[:,0]), np.max(seg2[:,0]))):
                od[:,1] = -5 # masked - not being in the time of interest

            elif mask is not None and len(tree_mask1.overlap(np.min(seg2[:,0]), np.max(seg2[:,0]))):
                od[:,1] = -4 # masked - audible call existing
            
            elif np.isnan(peaks[0,0]):
                od[:,1] = -2 # not localized 

            elif np.sum(np.isnan(snout_pos[i_frame, :])) == 0:

                snouts = np.reshape(snout_pos[i_frame, :], [-1, 2])
                tau = get_tau(data_dir, calibfile, speedOfSound, points=snouts)
                s, _ = calc_seg_power(fp_dat, seg2, tau, fs, n_ch, pressure_calib, return_average=False)

                az, el = point2angle(data_dir, calibfile, snouts)
                snouts_a = np.array([az, el]).T
                
                az, el = point2angle(data_dir, calibfile, peaks)
                peaks_a = np.array([az, el]).T

                snouts_a_aligned = np.zeros(snouts_a.shape)
                for i in range(n_mice):
                    d = snouts_a[i,:] - peaks_a
                    d = np.linalg.norm(d, axis=1)
                    I = np.nanargmin(d)
                    snouts_a_aligned[i,:] = snouts_a[i,:] - peaks_a[I,:]

                d = np.linalg.norm(snouts_a_aligned, axis=1)
                n_in_area = np.sum(d<d_thr)

                if n_in_area == 0:
                    od[:,1] = -1 # not assigned (n.s.)
                elif n_in_area == 1:
                    od[:,1] = np.argmin(d)
                else:
                    s_mean = np.median(s, axis=0)
                    I_b2 = np.argsort(s_mean)
                    I_b2 = I_b2[-2:]   # best 2 id

                    d = np.linalg.norm(snouts_a_aligned[I_b2[0],:] - snouts_a_aligned[I_b2[1],:])
                    _, p = scipy.stats.wilcoxon(s[:,I_b2[0]], s[:,I_b2[1]])

                    f = scipy.interpolate.interp2d(D_clut, P_clut, clut[n_in_area-2], kind='linear')
                    assign_confidence = f(d,p)

                    od[:,1] = I_b2[1]
                    od[:,2] = I_b2[0]
                    od[:,3] = assign_confidence
                    od[:,4] = d
                    od[:,5] = p

            else:
                od[:,1] = -3 # no tracking data 

            outdata.append(od)

    np.savetxt(data_dir + '/assign.csv', np.concatenate(outdata, axis=0), delimiter=',',
                header=' segment_id,1st_place_mouse,2nd_place_mouse,confidence,snout_distance,p_value,time,freq,ampitude,ID_B,ID_A')

def adjust_color(img):
    img2 = copy.deepcopy(img)
    img2 = cv2.cvtColor(img2, cv2.COLOR_RGB2GRAY)
    img2 = cv2.equalizeHist(img2)
    img2 = cv2.cvtColor(img2, cv2.COLOR_GRAY2RGB)
    return img2

def calc_conf_lut(R, D, show_figs=True):

    d = D[1:]
    p = abs(R)
    c = R>0

    thr = 10**np.arange(0, -50.1, -0.1)
    err_cnt = np.zeros([thr.shape[0], d.shape[0]])
    hit_cnt = np.zeros([thr.shape[0], d.shape[0]])

    for i_thr in range(thr.shape[0]):
        for i_d in range(d.shape[0]):
            I = p[:,i_d] < thr[i_thr]
            err_cnt[i_thr, i_d] = np.sum(np.logical_and(I, c[:, i_d]==0))
            hit_cnt[i_thr, i_d] = np.sum(np.logical_and(I, c[:, i_d]==1))

    A = scipy.ndimage.gaussian_filter(hit_cnt, [10, 1]) / scipy.ndimage.gaussian_filter(hit_cnt+err_cnt, [10, 1])
    I = err_cnt+hit_cnt < 500
    A[I] = np.nan

    for i_d in range(A.shape[1]):
        for i_p in range(A.shape[0]):
            a = copy.deepcopy(A[0:(i_p+1), i_d])
            a_max = np.nanmax(a)
            A[i_p, i_d] = a_max

    A[np.isnan(A)] = 0
    conf_lut = A

    if show_figs:
        plt.pcolormesh(rad2deg(d), np.log10(thr), conf_lut, cmap='jet', vmin=0.9, vmax=1.0)
        plt.colorbar()
        plt.xlabel('distance (deg)')
        plt.ylabel('threshold p value (log10)')
        plt.show()

        conf_lut_raw = hit_cnt / (hit_cnt+err_cnt)
        I = err_cnt+hit_cnt < 500
        conf_lut_raw[I] = np.nan
        plt.pcolormesh(rad2deg(d), np.log10(thr), conf_lut_raw, cmap='jet', vmin=0.9, vmax=1.0)
        plt.colorbar()
        plt.xlabel('distance (deg)')
        plt.ylabel('threshold p value (log10)')
        plt.show()

    return thr, d, conf_lut

def calc_flattened_spec(x, fs, med=None):

    fftsize = 512
    timestep = 0.0005

    mtsp = usvseg.multitaperspec(x,fs,fftsize,timestep)

    if med is None:
        fltnd, med = usvseg.flattening(mtsp, liftercutoff=6)
    else:
        fltnd, med = usvseg.flattening(mtsp, med, liftercutoff=6)

    return fltnd, med

def calib_micpos(data_dir, SEG, P, calibfile=None, h5f_outpath=None, pos_lim=None, pos_searchlim=None, pos_init=None, vis_progress=False, ftol=1e-3, maxiter=100):

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        speedOfSound = f['/misc/speedOfSound'][()]
        pressure_calib = f['/misc/pressure_calib'][()]
        mic0pos = f['/misc/mic0pos'][()]

    n_mic = n_ch

    with open_dat(data_dir) as fp_dat:

        def get_error(dx, P, SEG, data_dir, mic0pos, speedOfSound):

            dx = np.reshape(dx, (n_mic-1,3))
            dx = np.vstack([np.array([0,0,0]), dx])
            dx = dx + mic0pos

            n_pos_check = len(SEG)
            pwr = np.zeros(n_pos_check)

            for i_seg in range(n_pos_check):

                seg2 = copy.deepcopy(SEG[i_seg])

                pt = np.array([P[i_seg]])

                tau = get_tau(data_dir, None, speedOfSound, points=pt, micpos=dx)

                S, _ = calc_seg_power(fp_dat, seg2, tau, fs, n_ch, pressure_calib, return_average=True)

                pwr[i_seg] = S

            e = -np.sum(pwr**2)

            return e

        global Nfeval
        Nfeval = 1
        def callbackF(Xi):
            global Nfeval
            e = get_error(Xi, P, SEG, data_dir, mic0pos, speedOfSound)
            print('iter:{0:4d}, f(x) = '.format(Nfeval) + str(-e))

            if vis_progress:
                xx = np.reshape(Xi, (n_mic-1,3))
                plt.plot(0, 0, 'x')
                plt.plot(xx[:,0], xx[:,1], 'o')
                plt.xlim(pos_lim[0][0], pos_lim[1][0])
                plt.ylim(pos_lim[0][1], pos_lim[1][1])
                plt.gca().invert_xaxis()
                plt.gca().invert_yaxis()
                plt.gca().set_box_aspect((pos_lim[1][1] - pos_lim[0][1]) / (pos_lim[1][0] - pos_lim[0][0]))
                plt.show()

            Nfeval += 1

        # run optimization
        if calibfile is None:
            if pos_init is None:
                dx0 = np.tile([0,0,0], (n_mic-1,1))
            else:
                dx0 = np.array(pos_init)
        else:
            with h5py.File(calibfile, mode='r') as f:
                dx0 = f['/result/micpos'][()]
            dx0 = dx0[1:,:] - dx0[0,:]

        if pos_lim is None:
            if pos_searchlim is None:
                lb = np.tile([-0.01, -0.01, -0.003], (n_mic-1,1))
                ub = np.tile([0.01, 0.01, 0.003], (n_mic-1,1))
            else:
                lb = dx0 + np.tile(pos_searchlim[0], (n_mic-1,1))
                ub = dx0 + np.tile(pos_searchlim[1], (n_mic-1,1))
        else:
            lb = np.tile(pos_lim[0], (n_mic-1,1))
            ub = np.tile(pos_lim[1], (n_mic-1,1))

        pos_lim = [np.min(lb, axis=0), np.max(ub, axis=0)]

        dx0 = dx0.flatten()
        lb = lb.flatten()
        ub = ub.flatten()

        b = scipy.optimize.Bounds(lb, ub)

        print('Running optimization...')
        R = scipy.optimize.minimize(get_error, x0=dx0, args=(P, SEG, data_dir, mic0pos, speedOfSound), method='L-BFGS-B', 
                                    bounds=b, callback=callbackF, options={'ftol': ftol, 'maxiter': maxiter})

        dx_pred = R.x
        dx_pred = np.reshape(dx_pred, (n_mic-1,3))
        dx_pred = np.vstack([np.array([0,0,0]), dx_pred])
        micpos = dx_pred + mic0pos
        print('estimated micpos:')
        print(micpos)

        if h5f_outpath is None:
            h5f_outpath = data_dir + '/micpos.h5'
        with h5py.File(h5f_outpath, mode='w') as f:
            f.create_dataset('/result/micpos', data = micpos)
            for i_seg in range(len(SEG)):
                f.create_dataset('/seg/seg{:06}'.format(i_seg), data = SEG[i_seg])
                f.create_dataset('/snout_pos/seg{:06}'.format(i_seg), data = P[i_seg])

        print('the result is saved in: ' + h5f_outpath)

def calc_seg_stft(fp_dat, seg, fs, n_ch, pressure_calib):
    
    mrgn = 0.01
    t_intv = [np.min(seg[:,0]) - mrgn, np.max(seg[:,0]) + mrgn]

    if t_intv[0] < 0:
        t_intv[0] = 0

    n_ch  = n_ch.astype(np.int64)

    simg = read_dat(fp_dat, t_intv, fs, n_ch)

    simg = simg / pressure_calib

    nfft = 192
    noverlap = 0

    f, t, X = scipy.signal.stft(simg, fs=fs, window='hamming', nperseg=nfft, noverlap=noverlap, axis=0)

    return f, t, X

def calc_seg_power(fp_dat, seg, tau, fs, n_ch, pressure_calib, return_average=True, stft_data=None):

    mrgn = 0.01
    t_intv = [np.min(seg[:,0]) - mrgn, np.max(seg[:,0]) + mrgn]
    
    if t_intv[0] < 0:
        t_intv[0] = 0

    n_ch  = n_ch.astype(np.int64)

    if stft_data is None:
        f, t, X = calc_seg_stft(fp_dat, seg, fs, n_ch, pressure_calib)
    else:
        f, t, X = stft_data

    t = t+t_intv[0]

    I_t = np.searchsorted(t, seg[:,0])
    I_f = np.searchsorted(f, seg[:,1])-1

    if return_average:
        S = np.zeros(tau.shape[0])
    else:
        S = np.zeros((seg.shape[0]*2, tau.shape[0]))

    for i_segpoint in range(seg.shape[0]):
        i_f = I_f[i_segpoint]
        i_t = I_t[i_segpoint]

        x = X[i_f,:, i_t]
        if usegpu and tau.shape[0] > 100:
            s = calc_point_power_gpu(x, f[i_f], tau, n_ch)
        else:
            s = calc_point_power(x, f[i_f], tau, n_ch)

        if return_average:
            S = S + s
        else:
            S[i_segpoint*2, :] = s

        x = X[i_f+1,:, i_t]
        if usegpu and tau.shape[0] > 100:
            s = calc_point_power_gpu(x, f[i_f+1], tau, n_ch)
        else:
            s = calc_point_power(x, f[i_f+1], tau, n_ch)

        if return_average:
            S = S + s
        else:
            S[i_segpoint*2+1, :] = s

    if return_average:
        S = S/(seg.shape[0]*2)
        return S, None
    else:
        tf_point = np.repeat(seg[:,0:2],2, axis=0)
        return S, tf_point

def calc_point_power_gpu(x, f, tau, n_ch):

    xx = cupy.array(x)
    tau2 = cupy.array(tau)

    xx = cupy.reshape(xx, [n_ch, 1])
    xx = xx/cupy.linalg.norm(xx)

    a = -2j*np.pi*f
    svec = cupy.exp(a*tau2)
    svec = svec.T
    svec = svec/cupy.sqrt(n_ch)

    s = cupy.abs(cupy.dot(cupy.conjugate(xx.T), svec))
    s = s*s

    return cupy.asnumpy(s)

def calc_point_power(x, f, tau, n_ch):

    x = np.reshape(x, [n_ch, 1])
    x = x/np.linalg.norm(x)

    a = -2j*np.pi*f
    svec = ne.evaluate('exp(a*tau)')
    svec = svec.T
    svec = svec/np.sqrt(n_ch)

    s = np.abs(np.dot(np.conjugate(x.T), svec))
    s = ne.evaluate('s*s')

    return s

def calc_sspec_all_frame(data_dir, calibfile, fpath_out, t_end=-1, d=5):

    #with open(config_path, 'r') as f:
    #    usvcam_cfg = yaml.safe_load(f)
    #speedOfSound = float(usvcam_cfg['speed_of_sound'])
    #pressure_calib = np.array(usvcam_cfg['sound_pressure_calibration'])

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        im_w = f['/camera_param/color_intrin/width'][()]
        im_h = f['/camera_param/color_intrin/height'][()]
        speedOfSound = f['/misc/speedOfSound'][()]
        pressure_calib = f['/misc/pressure_calib'][()]
    if '/camera_param/cam_delay' in f:
        cam_delay = f['/camera_param/cam_delay'][()]
    else:
        cam_delay = cam_delay_default

    tau, grid_shape = get_tau(data_dir, calibfile, speedOfSound, d=d)

    seg = load_usvsegdata_ss(data_dir)

    vid_size = [im_w, im_h]

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    if t_end > 0:
        I = T[:,1]/fs < t_end
        T = T[I,:]

    n_frame = T.shape[0]
    wsize = fs/10

    def calc_mean_S_in_intv(buf, t0, t1):
        """
        buf: deque of (tf, S)
        tf.shape = (m, 2); tf[:,0] is time  
        S.shape  = (m, X)
        """
        S_sum = None
        count = 0

        for tf, S in buf:
            I = (tf[:, 0] > t0) & (tf[:, 0] < t1)

            c = np.sum(I)
            if c == 0:
                continue

            m = np.sum(S[I, :], axis=0)

            S_sum = m if S_sum is None else (S_sum + m)
            count += c

        if count == 0:
            return None, 0
        
        return S_sum / count, count

    t_calc_done = 0.0
    dt_calc = 2.0

    buf = deque()

    with h5py.File(fpath_out, mode='w') as fp_out:
        with open_dat(data_dir) as fp_dat:
            for i_frame in tqdm(range(n_frame)):

                t_intv = np.array([T[i_frame,1]-wsize/2, T[i_frame,1]+wsize/2])/fs - cam_delay

                while t_intv[1] > t_calc_done:
                    #print('calc chunk...')
                    calc_intv = [t_calc_done, t_calc_done+dt_calc]
                    I = np.logical_and(seg[:,0] > calc_intv[0], seg[:,0] < calc_intv[1])
                    if np.sum(I) > 0:
                        seg2 = seg[I,:]
                        S, tf = calc_seg_power(fp_dat, seg2, tau, fs, n_ch, pressure_calib, return_average=False)
                        #S = np.reshape(S, (-1, grid_shape[0], grid_shape[1]))
                        buf.append((tf, S))

                    t_calc_done += dt_calc
                    #print('done')
                
                while len(buf) > 3:
                    buf.popleft()       # delete old memory

                S, n_used = calc_mean_S_in_intv(buf, t_intv[0], t_intv[1])
                if n_used < 3 * 2: # should include >= 3 points in seg
                    continue
                S = np.reshape(S, grid_shape)
                fp_out.create_dataset('/sspec/' + '/{:06}'.format(i_frame), data = S)

def calc_vm_stats(data_dir, calibfile, roi=None, iter_id=None):

    D = np.arange(0.0, 15, 0.5) /180 * np.pi
    n_trial = 7

    #with open(config_path, 'r') as f:
    #    usvcam_cfg = yaml.safe_load(f)
    #speedOfSound = float(usvcam_cfg['speed_of_sound'])
    #pressure_calib = np.array(usvcam_cfg['sound_pressure_calibration'])

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        im_w = f['/camera_param/color_intrin/width'][()]
        im_h = f['/camera_param/color_intrin/height'][()]
        speedOfSound = f['/misc/speedOfSound'][()]
        pressure_calib = f['/misc/pressure_calib'][()]
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    loc = np.genfromtxt(data_dir + '/loc.csv', delimiter=',')

    snout_pos = np.genfromtxt(data_dir + '/snout.csv', delimiter=',')
    snout_pos = snout_pos[:,:2]

    seg = load_usvsegdata_ss(data_dir)
    _, I_ss = np.unique(seg[:,3:5], axis=0, return_inverse=True)

    if iter_id is None:
        fpath_out = data_dir + '/vmstat.h5'
    else:
        fpath_out = data_dir + '/vmstat.{:d}.h5'.format(iter_id)

    def get_power_at_vm(data_dir, calibfile, seg2, fs, n_ch, pressure_calib, D, snout_pos, roi, speedOfSound, stft_data):
        angles, pts = gen_rand_points(data_dir, calibfile, D, snout_pos, roi)
        tau = get_tau(data_dir, calibfile, speedOfSound, points=pts)  
        s, _ = calc_seg_power(None, seg2, tau, fs, n_ch, pressure_calib, return_average=False, stft_data=stft_data)
        return s

    with h5py.File(fpath_out, mode='w') as fp_out:
        fp_out.create_dataset('/D', data = D)
        fp_out.create_dataset('/n_trial', data = n_trial)

        with open_dat(data_dir) as fp_dat:
            for i_ss in tqdm(range(np.max(I_ss)+1)):

                seg2 = seg[I_ss==i_ss,:]

                i_frame = time2vidframe(data_dir, (np.min(seg2[:,0])+np.max(seg2[:,0]))/2, T, fs, cam_delay)

                if np.isnan(snout_pos[i_frame, 0]):
                    continue

                peaks = np.reshape(loc[i_ss,3:], [-1, 2])
                if np.sum(np.logical_not(np.isnan(peaks[:,0]))) == 0:
                    #print('not localized')
                    continue
                
                stft_data = calc_seg_stft(fp_dat, seg2, fs, n_ch, pressure_calib)

                S = np.zeros([seg2.shape[0]*2, D.shape[0], n_trial])
                result = joblib.Parallel(n_jobs=-1)(
                            joblib.delayed(get_power_at_vm)(data_dir, calibfile, seg2, fs, n_ch, pressure_calib, D, snout_pos[i_frame,:], roi, speedOfSound, stft_data) for i_trial in range(n_trial))
                
                for i_trial in range(n_trial):
                    S[:, :, i_trial] = result[i_trial]

                fp_out.create_dataset('/rslt/seg{:06}'.format(i_ss), data = S)
                fp_out.create_dataset('/segdata/seg{:06}'.format(i_ss), data = seg2)

def clean_data_dir(data_dir, filekeep=[]):

    remove_folders = ['loc', 'seg2']
    remove_files = ['loc.csv', 'vmstat.*.h5', 'vid.loc.mp4', 
                    'assign.csv', 'result.csv', 'vid.asgn.mp4']

    for rf in remove_folders:
        if rf in filekeep:
            continue
        L = glob.glob(data_dir + '/' + rf)
        for fn in L:
            shutil.rmtree(fn)

    for rf in remove_files:
        if rf in filekeep:
            continue
        L = glob.glob(data_dir + '/' + rf)
        for fn in L:
            os.remove(fn)

def create_assignment_video(data_dir, n_mice, color_eq=False, min_dist=2.0):

    audiblewavfile = glob.glob(data_dir + '/*.audible.wav')
    if len(audiblewavfile) == 0:
        audiblewavfile = glob.glob(data_dir + '/*.audible.flac')
    audiblewavfile = audiblewavfile[0]
    
    os.makedirs('./tmp', exist_ok=True)
    tmpvidfile = './tmp/tmp.mp4'

    print('making video with assignment...')
    draw_assign_on_all_vidframe(tmpvidfile, data_dir, n_mice, color_eq=color_eq, min_dist=min_dist)

    print('combining sound and video...')
    outfile = data_dir + './vid.asgn.mp4'
    os.system('ffmpeg -y -i "' + tmpvidfile + '" -i "' + audiblewavfile + '" -c:v copy -c:a aac "' + outfile + '"')

    print('done')

def create_localization_video(data_dir, calibfile, t_end=-1, color_eq=False, only_max=None, min_peak_lev=1.6):
    
    if only_max is None:
        only_max = only_max_default

    audiblewavfile = glob.glob(data_dir + '/*.audible.wav')
    if len(audiblewavfile) == 0:
        audiblewavfile = glob.glob(data_dir + '/*.audible.flac')
    audiblewavfile = audiblewavfile[0]

    os.makedirs('./tmp', exist_ok=True)
    tmpvidfile = './tmp/tmp.mp4'
    tmpsspecfile = './tmp/tmp_sspec.h5'

    print('processing: ', data_dir)
    L = glob.glob(data_dir + '/seg/*.csv')
    if len(L) == 0:
        print('no file was found in seg folder.')
        return

    print('calculating spatial spectrum in each frame...')
    calc_sspec_all_frame(data_dir, calibfile, tmpsspecfile, t_end=t_end)

    print('overlay spatial spectrum on video...')
    draw_spect_on_all_vidframe(tmpvidfile, data_dir, tmpsspecfile, t_end=t_end, color_eq=color_eq, only_max=only_max, min_peak_lev=min_peak_lev)

    print('combining sound and video...')
    if t_end > 0:
        tmpwavfile = './tmp/tmp.wav'
        os.system('ffmpeg -y -i "' + audiblewavfile + '" -t ' + str(t_end) + ' "' + tmpwavfile + '"')
    else:
        tmpwavfile = audiblewavfile

    outfile = data_dir + './vid.loc.mp4'
    os.system('ffmpeg -y -i "' + tmpvidfile + '" -i "' + tmpwavfile + '" -c:v copy -c:a aac "' + outfile + '"')

    print('done')

def dat2wav(data_dir, i_ch, offset=0):

    fpath_dat = data_dir + '/snd.dat'
    fpath_wav = data_dir + '/snd.ch{:d}.wav'.format(i_ch)
    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]

    fsize = os.path.getsize(fpath_dat)

    sw = 2
    readsize = fs

    cnt = 0
    with open(fpath_dat, 'rb') as f:
        with wave.open(fpath_wav, 'wb') as f_out:
            f_out.setnchannels(1)
            f_out.setsampwidth(sw)
            f_out.setframerate(fs)
            f_out.setcomptype('NONE', 'not compressed')
            while True:
                a = np.fromfile(f, np.int16, n_ch*readsize)
                if len(a) == 0:
                    break
                a = a - offset
                x = a.reshape([-1, n_ch])
                xx = x[:,i_ch]
                xxx = np.zeros(xx.shape[0], np.int16)
                xxx[:] = xx
                f_out.writeframes(xxx)

def dat2flac(data_dir, i_ch, offset=0):

    fpath_dat = data_dir + '/snd.dat'
    fpath_flac = data_dir + '/snd.ch{:d}.flac'.format(i_ch)
    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]

    fsize = os.path.getsize(fpath_dat)

    a = np.fromfile(fpath_dat, np.int16)
    a = a - offset
    x = a.reshape([-1, n_ch])

    soundfile.write(fpath_flac, x[:,i_ch], fs)

def flac2flac(data_dir, i_ch, offset=0):

    fpath_flac_in = data_dir + '/snd.flac'
    fpath_flac_out = data_dir + '/snd.ch{:d}.flac'.format(i_ch)
    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]

    with soundfile.SoundFile(fpath_flac_in, "r") as f:
        with soundfile.SoundFile(fpath_flac_out, "w", samplerate=f.samplerate, channels=1, subtype="PCM_16") as g:
            for block in soundfile.blocks(fpath_flac_in, blocksize=262144, always_2d=True, dtype="int16"):
                g.write(block[:, i_ch] - offset)

def detect_bbv(data_dir, thr=1.4):

    # Note: thr = 1.4 -> 28 dB

    # load wav
    L = glob.glob(data_dir + '/snd.*.wav')
    for l in L:
        if 'audible' in l:
            continue
        wavfile = l
    fs, x = scipy.io.wavfile.read(wavfile)

    # downsample
    r = 10
    x = scipy.signal.decimate(x, r)
    fs = fs/r

    # spectrogram audible band
    f,t,s = scipy.signal.spectrogram(x, fs, nperseg=int(fs/100))

    I = np.logical_and(f>2000, f<16000)
    f = f[I]
    s = s[I,:]
    s = np.log10(s+1e-20)

    # filtering spectrogram
    med = np.median(s,axis=1)
    s = s - np.tile(med, [s.shape[1],1]).T
    s[s<0] = 0
    s = scipy.ndimage.median_filter(s, [5, 1])
    s = s - np.median(s, axis=0)

    # find interval above threshold
    intv = to_intervals(np.max(s, axis=0)>thr, t)

    # merge intervals
    intv[:,0] -= 0.025
    intv[:,1] += 0.025

    t_ori = np.arange(0, x.shape[0]/fs, 1/fs)

    tree = intervaltree.IntervalTree.from_tuples(intv)
    tree.merge_overlaps(strict=False)
    tree = sorted(tree)
    tree = list(tree)

    mrg_intv = np.zeros((len(tree), 2))
    for i in range(len(tree)):
        ii = np.array(list(tree[i]))
        mrg_intv[i,:] = ii[0:2]
    mrg_intv[:,0] += 0.02
    mrg_intv[:,1] -= 0.02

    I = mrg_intv[:,1] - mrg_intv[:,0] >= 0.01
    mrg_intv = mrg_intv[I,:]

    print(mrg_intv.shape[0], 'BBVs were found; total time is', np.sum(mrg_intv[:,1]-mrg_intv[:,0]))

    # draw fine spectrogram around intervals
    def fun(x, t_ori, intv):
        I = np.logical_and(t_ori>intv[0]-0.05, t_ori<intv[1]+0.05)
        xx = x[I]
        ff,tt,ss = scipy.signal.spectrogram(xx, fs, window=('hamming'), nperseg=int(fs/100), noverlap=int(fs/100 * 0.9))

        return np.log10(ss[20:160,:]+1e-20)

    specs = joblib.Parallel(n_jobs=-1)(joblib.delayed(fun)(x, t_ori, intv) for intv in mrg_intv)

    # save spectrogram as images
    out_dir = data_dir + '/bbv_detection'
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    for i_start in range(0, len(x), int(fs*5)):
        xx = x[i_start:(i_start+int(fs*5))]
        ff,tt,ss = scipy.signal.spectrogram(xx, fs, window=('hamming'), nperseg=int(fs/100), noverlap=int(fs/100 * 0.9))
        ss = np.log10(ss[20:160,:]+1e-20)
        ttt = i_start/fs + tt
        marker = np.ones([20, ss.shape[1]])*0.5
        for intv in mrg_intv:
            I = np.logical_and(ttt >= intv[0], ttt <= intv[1])
            marker[:,I] = 1
        marker *= np.max(ss)
        ss = np.concatenate([marker,ss],axis=0)
        img = ss
        img[img<0] = 0
        img = 1.0-img/np.max(img)
        img = np.flipud(img)
        img = cv2.resize(img,(1000,200))
        img = (img * 255).astype(np.uint8)
        img = cv2.cvtColor(img,cv2.COLOR_GRAY2BGR)
        cv2.imwrite(out_dir + '/{:04d}.jpg'.format(int(i_start/fs)), img)

    # save intervals
    mask = np.concatenate([mrg_intv, np.ones([mrg_intv.shape[0],1])], axis=1)
    np.savetxt(data_dir + '/mask.csv', mask, delimiter=',', header='onset, offset, type')

def draw_assign_on_all_vidframe(fpath_out, data_dir, n_mice, conf_thr=0.99, color_eq=False, min_dist=2.0):

    spect_height = 100

    with open(config_path, 'r') as f:
        usvcam_cfg = yaml.safe_load(f)

    ch_monitor = usvcam_cfg['ch_monitor']

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    vid_file = data_dir + '/vid.mp4'

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        n_ch = n_ch.astype(np.int64)
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default
        if '/camera_param/v_fps' in f:
            v_fps = f['/camera_param/v_fps'][()]
        else:
            v_fps = v_fps_default

    snout_pos = np.genfromtxt(data_dir + '/snout.csv', delimiter=',')

    T[:,1] = T[:,1]/fs
    t = np.arange(0, T[-10,1], 1/v_fps)
    n_frame = t.shape[0]

    vr = cv2.VideoCapture(vid_file)
    fmt = cv2.VideoWriter_fourcc('m', 'p', '4', 'v') 
    vw = cv2.VideoWriter(fpath_out, fmt, v_fps, (int(vr.get(cv2.CAP_PROP_FRAME_WIDTH)), int(vr.get(cv2.CAP_PROP_FRAME_HEIGHT)+spect_height*2)), isColor=True) 

    seg = np.genfromtxt(data_dir + '/assign.csv', delimiter=',')

    I = np.logical_or(seg[:,4] <= min_dist/180*np.pi, seg[:,3] < conf_thr)
    seg[I,1] = -1

    mask = None
    if os.path.exists(data_dir + '/mask.csv'):
        mask = np.loadtxt(data_dir + '/mask.csv', delimiter=',')
        print('mask data is found', flush=True)

    with open_dat(data_dir) as fp_dat:
        
        crnt_ts = -10
        cnt = 0  
        pre_t = -1
        med = None
        for i_frame in tqdm(range(n_frame)):

            while crnt_ts-cam_delay < t[i_frame]:
                ret, frame = vr.read()
                if color_eq:
                    frame = adjust_color(frame)
                crnt_ts = T[cnt,1]
                cnt += 1

            #### draw spectrogram with assignment
            if t[i_frame] - pre_t >= 1:
                pre_t = int(t[i_frame])
                t_intv = [pre_t, pre_t+1]

                #spec, med = get_sound_spec(fp_dat, t_intv, fs, n_ch, ch_monitor, med)
                spec, _ = get_sound_spec(fp_dat, t_intv, fs, n_ch, ch_monitor)

                I = np.logical_and(seg[:,6] >= t_intv[0], seg[:,6] <= t_intv[1])
                seg2 = seg[I, :]
                _, I_ss = np.unique(seg2[:,0], return_inverse=True)

                dt = 0.0005
                fftsize = 512
                clrs = [(255,255,0),(255,0,255),(0,255,255),(0,255,0),(0,128,255),(0,0,255),(255,0,128)]

                spec_show = spec*255
                spec_show = cv2.cvtColor(spec_show.astype(np.uint8), cv2.COLOR_GRAY2BGR)
                spec_ori_show = copy.deepcopy(spec_show)

                if I_ss.shape[0] > 0:
                    for i_ss in range(np.max(I_ss)+1):
                        seg3 = seg2[I_ss==i_ss,:]
                        i_t = np.round((seg3[:,6] - t_intv[0])/dt)
                        i_f = np.round(seg3[:,7]/(fs/2) * (fftsize/2))
                        pt = np.array([i_t, i_f]).T
                        I = np.argsort(i_t)
                        pt = pt[I,:]
                        cv2.polylines(spec_show, [pt.astype(int)], False, (0,0,0), thickness=10)
                        i_mouse = seg3[0,1]
                        if i_mouse >= 0:
                            clr = clrs[int(i_mouse)]
                            cv2.polylines(spec_show, [pt.astype(int)], False, clr, thickness=4)

                spec_show_masked = copy.deepcopy(spec_show)
                if mask is not None:
                    for intv in mask:
                        if intv[2] == 1:
                            x1 = (intv[0] - t_intv[0])/(t_intv[1]-t_intv[0]) * spec_show.shape[1]
                            x2 = (intv[1] - t_intv[0])/(t_intv[1]-t_intv[0]) * spec_show.shape[1]
                            cv2.rectangle(spec_show_masked, (int(x1), 0), (int(x2), spec_show.shape[0]), (0,0,50), thickness=-1)
                    
                    for intv in mask:
                        if intv[2] == 2:
                            x1 = (intv[0] - t_intv[0])/(t_intv[1]-t_intv[0]) * spec_show.shape[1]
                            x2 = (intv[1] - t_intv[0])/(t_intv[1]-t_intv[0]) * spec_show.shape[1]
                            cv2.rectangle(spec_show_masked, (int(x1), 0), (int(x2), spec_show.shape[0]), (0,0,0), thickness=-1)

                    
                spec_show = cv2.addWeighted(spec_show, 0.5, spec_show_masked, 0.5, 0)
                            
            #### draw snout markers
            t_intv = [t[i_frame]-(1/v_fps/2), t[i_frame]+(1/v_fps/2)]
            I = np.logical_and(seg[:,6] >= t_intv[0], seg[:,6] <= t_intv[1])
            seg2 = seg[I, :]
            if snout_pos.shape[0] >= cnt:
                sp = np.reshape(snout_pos[cnt-1, :], [-1, 2])
            else:
                sp = np.zeros([n_mice, 2])
                sp[:,:] = np.nan
                
            for i_mouse in range(n_mice):
                if np.sum(np.isnan(sp[i_mouse, :]))==0:
                    cv2.circle(frame, (int(sp[i_mouse, 0]), int(sp[i_mouse, 1])), 3, clrs[i_mouse], thickness=-1, lineType=cv2.LINE_AA)
                    if np.sum(seg2[:,1]==i_mouse) > 0:
                        cv2.circle(frame, (int(sp[i_mouse, 0]), int(sp[i_mouse, 1])), 20, clrs[i_mouse], thickness=2, lineType=cv2.LINE_AA)

            #### integrate video frame with spectrograms and output to video
            w = vr.get(cv2.CAP_PROP_FRAME_WIDTH)
            spec_show2 = cv2.resize(spec_show, (int(w), spect_height))
            spec_ori_show2 = cv2.resize(spec_ori_show, (int(w), spect_height))
            spec_show2 = np.flip(spec_show2, axis=0)
            spec_ori_show2 = np.flip(spec_ori_show2, axis=0)
            d = round(w*(t[i_frame] - pre_t))
            spec_show2[:,d,:] = 0
            spec_ori_show2[:,d,:] = 0

            frame_out = np.concatenate([frame, spec_show2, spec_ori_show2], axis=0)
            frame_out[frame.shape[0]-1,:,:] = 0
            frame_out[frame.shape[0]+spec_ori_show2.shape[0]-1,:,:] = 0

            vw.write(frame_out)

    vw.release()
    vr.release()

def draw_spect_on_all_vidframe(fpath_out, data_dir, sspecfile, t_end=-1, color_eq=False, only_max=None, min_peak_lev=1.6):
    
    if only_max is None:
        only_max = only_max_default

    spect_height = 100

    with open(config_path, 'r') as f:
        usvcam_cfg = yaml.safe_load(f)

    ch_monitor = usvcam_cfg['ch_monitor']

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    vid_file = data_dir + '/vid.mp4'
    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        n_ch = n_ch.astype(np.int64)
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default
        if '/camera_param/v_fps' in f:
            v_fps = f['/camera_param/v_fps'][()]
        else:
            v_fps = v_fps_default

    if t_end > 0:
        I = T[:,1]/fs < t_end
        T = T[I,:]

    T[:,1] = T[:,1]/fs
    t = np.arange(0, T[-10,1], 1/v_fps)
    n_frame = t.shape[0]

    vr = cv2.VideoCapture(vid_file)
    fmt = cv2.VideoWriter_fourcc('m', 'p', '4', 'v') 
    vw = cv2.VideoWriter(fpath_out, fmt, v_fps, (int(vr.get(cv2.CAP_PROP_FRAME_WIDTH)), int(vr.get(cv2.CAP_PROP_FRAME_HEIGHT)+spect_height)), isColor=True) 
    with h5py.File(sspecfile, mode='r') as f:
        with open_dat(data_dir) as fp_dat:
            frame_snd = []  
            if '/sspec' in f:
                for key in f['/sspec'].keys():
                    frame_snd.append(int(key))
            crnt_ts = -10
            cnt = 0  

            spec = np.zeros((10,10))

            pre_t = -1
            med = None
            for i_frame in tqdm(range(n_frame)):
                while crnt_ts-cam_delay < t[i_frame]:
                    ret, frame = vr.read()
                    crnt_ts = T[cnt,1]
                    cnt += 1
                    
                if frame is None:
                    continue

                if t[i_frame] - pre_t >= 1:
                    pre_t = int(t[i_frame])

                    t_intv = [pre_t, pre_t+1]
                    spec, _ = get_sound_spec(fp_dat, t_intv, fs, n_ch, ch_monitor)
                    spec = np.flipud(spec)

                w = vr.get(cv2.CAP_PROP_FRAME_WIDTH)
                spec_show = cv2.resize(spec, (int(w), spect_height))
                d = round(w*(t[i_frame] - pre_t))
                spec_show = spec_show*255
                spec_show = cv2.cvtColor(spec_show.astype(np.uint8), cv2.COLOR_GRAY2BGR)
                spec_show[:,d,:] = 0

                if cnt-1 in frame_snd:
                    sspec = f['/sspec/' + '/{:06}'.format(cnt-1)][()]
                    img = draw_spec_on_vidframe(sspec, frame, color_eq=color_eq, only_max=only_max, min_peak_lev=min_peak_lev)
                    frame_out = np.concatenate([img.astype(np.uint8), spec_show], axis=0)
                    vw.write(frame_out)
                else:
                    if color_eq:
                        frame2 = adjust_color(frame)
                    else:
                        frame2 = frame
                    frame_out = np.concatenate([frame2, spec_show], axis=0)
                    vw.write(frame_out)

    vw.release()
    vr.release()

def draw_spec_on_vidframe(sspec, frame, vid_mrgn=0, color_eq=False, only_max=None, z_range=None, min_peak_lev=1.6):
    
    if only_max is None:
        only_max = only_max_default

    sspec_z = (sspec - np.mean(sspec)) / np.std(sspec)

    if z_range is None:
        z_range = [0, 0]
        z_range[0] = max(1.0, np.max(sspec_z)*0.9)
        z_range[1] = max(4.0, np.max(sspec_z))
    
    sspec_z_disp = (sspec_z - z_range[0]) / (z_range[1] - z_range[0])
    
    sspec_z_disp[sspec_z_disp < 0.0] = 0.0
    sspec_z_disp[sspec_z_disp > 1.0] = 1.0
    sspec_z_disp *= 255
    sspec_z_disp = cv2.applyColorMap(sspec_z_disp.astype(np.uint8), cv2.COLORMAP_JET)

    if color_eq:
        img1 = adjust_color(frame)
    else:
        img1 = frame 

    s = cv2.resize(sspec_z_disp, (frame.shape[1]+vid_mrgn*2, frame.shape[0]+vid_mrgn*2))
    if vid_mrgn > 0:
        s = s[vid_mrgn:-vid_mrgn, vid_mrgn:-vid_mrgn]

    img2 = img1 * 0.5 + s * 0.5
    
    mask = cv2.resize(sspec_z, (frame.shape[1]+vid_mrgn*2, frame.shape[0]+vid_mrgn*2)) > z_range[0]
    if vid_mrgn > 0:
        mask = mask[vid_mrgn:-vid_mrgn, vid_mrgn:-vid_mrgn]

    img3 = np.zeros(frame.shape)
    for i in range(3):  
        img3[:,:,i] = img2[:,:,i]*mask + img1[:,:,i]*np.logical_not(mask)
    
    img3 = img3.astype(np.uint8)

    # draw peak pos

    peaks = get_sspec_peaks(sspec, (frame.shape[1], frame.shape[0]), vid_mrgn=vid_mrgn, only_max=only_max, min_peak_lev=min_peak_lev)
    
    if peaks is not None:
        if peaks.shape[0] > 0:
            for i in range(peaks.shape[0]):
                cv2.drawMarker(img3, (int(peaks[i,0]), int(peaks[i,1])), (255, 255, 255), markerSize=10, thickness=2)
    
    return img3

class ErrDensityModel:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    @staticmethod
    def func(x, a, b):
        return a * np.exp(-x**2 * (0.5 / b**2))
    
    def __call__(self, x, scale=True):
        return self.func(x, self.a, self.b)

def estimate_assign_param_simple(data_dirs, calib_files, outfile, d_max=10, bin_count=200, show_figs=True):
    """
    This algorithm is considered beta.
    It is grounded in existing literature, but the present method and its
    implementation have not yet been subject to peer review.
    """

    ####### get list of [snout pos, sspec peak pos] (X)
    X = []

    for i_data, data_dir in enumerate(data_dirs): 

        snout_pos = np.genfromtxt(data_dir + '/snout.csv', delimiter=',')
        loc_data = np.genfromtxt(data_dir + '/loc.csv', delimiter=',', skip_header=1)

        calib_file = calib_files[i_data]

        for ld in loc_data:

            i_frame = int(ld[2])

            p = np.reshape(ld[3:], [-1, 2])

            p_s = snout_pos[i_frame,:]

            if np.isnan(p_s[0]):
                continue

            if np.isnan(p[0,0]):
                continue

            p_s = np.reshape(p_s, [1,2])

            p_az, p_el = point2angle(data_dir, calib_file, p)
            p_s_az, p_s_el = point2angle(data_dir, calib_file, p_s)

            p_angle = np.array([p_az, p_el]).T
            p_s_angle = np.array([p_s_az, p_s_el]).T

            d = p_angle - p_s_angle
            d = np.sqrt(np.sum(d**2, axis=1))
            i_peak = np.nanargmin(d)
            x = [rad2deg(p_az[i_peak]), rad2deg(p_el[i_peak]), rad2deg(p_s_az[0]), rad2deg(p_s_el[0])]
            X.append(x)

    X = np.squeeze(np.array(X))

    ##### calculate relative position (error) and estimate the error distribution 
    RP = []
    for i_data in range(X.shape[0]):
        rp = X[i_data, 0:2] - X[i_data, 2:4]
        RP.append(rp)
    RP = np.array(RP)

    D = np.linalg.norm(RP, axis=1)

    I = D < d_max
    
    D = D[I]
    RP = RP[I,:]

    dd = d_max/20

    # Define bin edges
    D_sorted = np.sort(D)
    cut_idx = np.arange(bin_count, D_sorted.size, bin_count)
    eds = np.r_[D_sorted[0], D_sorted[cut_idx], D_sorted[-1]]
    i = np.searchsorted(eds, d_max, side="right")
    eds = eds[:i+1]

    # Histogram counts
    C, _ = np.histogram(D, bins=eds)

    # Compute areas
    A = np.pi * (eds ** 2)
    A = np.diff(A)

    # Density calculation
    density = C / A

    # Midpoints of bins
    x = (eds[:-1] + eds[1:]) / 2

    y = density / np.sum(density)

    # Fit model definition
    def model_func(x, a, b):
        return a * np.exp(-x**2 * (0.5 / b**2))

    # Fit the data with constraints

    bounds=([0, 0], [np.inf, d_max]) 

    popt, pcov = scipy.optimize.curve_fit(
        ErrDensityModel.func,
        x,
        y,
        bounds=bounds,
        method='trf'  # Trust Region Reflective algorithm suitable for bounded problems
    )

    a, b = popt
    fitted_model = ErrDensityModel(a,b)

    if show_figs:
        # Create finer x for smooth plotting
        xx = np.arange(0, d_max, 0.1)
        yy = fitted_model(xx)

        # Plot original density points
        plt.figure(figsize=(8, 6))
        plt.plot(x, density, 'o', label='Data')

        # Plot the fitted curve scaled back
        plt.plot(xx, yy * np.sum(density), label='Fit')
        plt.xlabel('Distance (deg)')
        plt.ylabel('Density (1/deg^2)')
        plt.legend()
        plt.grid(True)
        plt.show()

    print('median error = {:f} deg, N of segments = {:d}'.format(np.median(D), D.shape[0]))

    with h5py.File(outfile, mode='w') as f_out:
        f_out.create_dataset('/model', data = 'simple')
        f_out.create_dataset('/a', data = a)
        f_out.create_dataset('/b', data = b)
    
    print(a, b)
    print('done')
    
def estimate_assign_param(data_dirs, calibfiles, outfile, n_trial=7, iter_ids=None, show_figs=True):

    with h5py.File(outfile, mode='w') as f_out:
        
        ##### distance error distribution for screening  ######
        print('calculating distance error distribution...')

        Da = []
        Dp = []
        for i_data, data_dir in enumerate(data_dirs):

            loc = np.genfromtxt(data_dir + '/loc.csv', delimiter=',')
            snout_pos = np.genfromtxt(data_dir + '/snout.csv', delimiter=',')
            snout_pos = snout_pos[:,:2]

            snout_pos = snout_pos[loc[:,2].astype(int), :]

            az, el = point2angle(data_dir, calibfiles[i_data], snout_pos)
            snout_pos_angle = np.array([az, el]).T

            for i_seg in tqdm(range(loc.shape[0])):
                peaks = np.reshape(loc[i_seg,3:], [-1, 2])
                if np.isnan(peaks[0,0]):
                    continue
                if np.isnan(snout_pos[i_seg,0]):
                    continue
                az, el = point2angle(data_dir, calibfiles[i_data], peaks)
                peaks_angle = np.array([az, el]).T
                d = peaks_angle - snout_pos_angle[i_seg,:]
                d = np.sqrt(np.sum(d**2, axis=1))
                d = np.nanmin(d)
                Da.append(d)
                d = peaks - snout_pos[i_seg,:]
                d = np.sqrt(np.sum(d**2, axis=1))
                d = np.nanmin(d)
                Dp.append(d)

        Da = np.array(Da)
        Dp = np.array(Dp)

        d95 = np.percentile(Da, 95)
        d99 = np.percentile(Da, 99)

        f_out.create_dataset('/dist_err_dist_in_rad', data = Da)
        f_out.create_dataset('/dist_err_dist_in_px', data = Dp)

        if show_figs:
            print('distance between snout and peak, 75, 95 and 99 percentile: {:.2f} deg, {:.2f} deg, {:.2f} deg'.format(rad2deg(np.percentile(Da, 75)), rad2deg(d95), rad2deg(d99)))
            print('median = {:.2f} deg'.format(rad2deg(np.percentile(Da, 50))))
            plt.hist(rad2deg(Da), bins=500, range=[0, 40])
            plt.xlabel('error (deg)')
            plt.show()

            print('in px, 75, 95 and 99 percentile: {:.2f} px, {:.2f} px, {:.2f} px'.format(np.percentile(Dp, 75), np.percentile(Dp, 95), np.percentile(Dp, 99)))
            print('median = {:.2f} px'.format(np.percentile(Dp, 50)))
            plt.hist(Dp, bins=500, range=[0, 500])
            plt.xlabel('error (px)')
            plt.show()

        print('done')
        
        ##### LUT of confidence #############
        print('calculating the lookup table for confidence estimation...')

        if iter_ids is None:
            iter_ids = [0]
            flag_no_itr = True
        else:
            flag_no_itr = False

        for n_other_mouse in range(1,n_trial+1):

            print('checking n = {:d}...'.format(n_other_mouse+1), flush=True)

            R = []
            for data_dir in data_dirs:

                for iter_id in iter_ids:
                    if flag_no_itr:
                        vmstatfile = data_dir + '/vmstat.h5'
                    else:
                        vmstatfile = data_dir + '/vmstat.{:d}.h5'.format(iter_id)
                    with h5py.File(vmstatfile, mode='r') as f:
                        D = f['/D'][()]
                        keys = f['/segdata'].keys()
                        for i_k, k in enumerate(tqdm(keys)):
                            sd = f['/segdata/' + k][()]

                            r = f['/rslt/' + k][()]

                            # add real mouse to the virtual mouse data
                            r0 = np.reshape(r[:,0,0], [-1, 1, 1])
                            r0 = np.tile(r0, [1, r.shape[1], 1])
                            r2 = np.concatenate([r0, r], axis=2)

                            # select mice according to N of mice
                            r2 = r2[:,:,0:n_other_mouse+1]

                            # find best 2 from the group
                            r_mean = np.median(r2, axis=0)
                            I = np.argsort(r_mean, axis=1)
                            I = I[:,-2:]

                            # do stat test between the best 2, p become minus when real mouse is not the best
                            rr = np.zeros(D.shape[0]-1)
                            for i_d in range(D.shape[0]-1):
                                _, p = scipy.stats.wilcoxon(r2[:,i_d+1, I[i_d+1, 0]], r2[:,i_d+1, I[i_d+1, 1]])
                                if I[i_d+1, 1] == 0:    # real mouse is the best
                                    rr[i_d] = p
                                else:
                                    rr[i_d] = -p

                            R.append(rr)

            R = np.array(R)

            R = R[np.logical_not(np.isnan(R[:,0])),:]

            p, d, conf_lut = calc_conf_lut(R, D, show_figs)

            f_out.create_dataset('/conf/N_other_{:d}/lut'.format(n_other_mouse), data = conf_lut)
            f_out.create_dataset('/conf/N_other_{:d}/d'.format(n_other_mouse), data = d)
            f_out.create_dataset('/conf/N_other_{:d}/p'.format(n_other_mouse), data = p)

            print('done')

def find_nearest_pos(candidate_pos, ref_pos):
    if len(candidate_pos.shape) == 2:
        d = np.sqrt(np.sum((candidate_pos-ref_pos)**2, axis=1))
        p = candidate_pos[np.argmin(d),:]
        d = np.min(d)
    else:
        d = np.sqrt(np.sum((candidate_pos-ref_pos)**2))
        p = candidate_pos
    return p, d

def gen_rand_points(data_dir, calibfile, D, p0, roi):

    az0, el0 = point2angle(data_dir, calibfile, np.array([p0]))
    
    A = []
    P = []
    for i_d in range(D.shape[0]):
        while True:
            theta = np.random.rand()*2*np.pi
            a = [D[i_d]*np.cos(theta)+az0[0], D[i_d]*np.sin(theta)+el0[0]]
            p = angle2point(data_dir,calibfile, np.array([a[0]]), np.array([a[1]]))
            p = p.ravel()

            if roi is None:
                A.append(a)
                P.append(p)
                break
            elif all([p[0] >= roi[0], p[0] <= roi[2], p[1] >= roi[1], p[1] <= roi[3]]):
                A.append(a)
                P.append(p)
                break

    A = np.array(A)
    P = np.array(P)
    return A, P

def get_sound_spec(fp_dat, t_intv, fs, n_ch, tgt_ch, med=None):

    if t_intv[0] < 0:
        t_intv[0] = 0

    simg = read_dat(fp_dat, t_intv, fs, n_ch)
    
    x = simg[:,tgt_ch]

    b, a = scipy.signal.butter(1, 1000/fs/2, btype='high')
    x = scipy.signal.filtfilt(b, a, x)

    fltnd, med = calc_flattened_spec(x, fs, med)
    fltnd[fltnd>30] = 30
    fltnd[fltnd<0] = 0
    fltnd /= 30
    spec = 1-fltnd

    return spec, med

def get_sspec_peaks(sspec, vid_size, vid_mrgn=0, roi=None, only_max=None, min_peak_lev=1.6):   # z=1.6 means around 95% in cum dist
        
    if only_max is None:
        only_max = only_max_default

    peaks = None

    if only_max:

        if (np.max(sspec)-np.mean(sspec))/np.std(sspec) > min_peak_lev:
            I = np.unravel_index(np.argmax(sspec), sspec.shape)
            peaks = np.array([I[1], I[0]], dtype=float)
            peaks = np.reshape(peaks, [1,2])
            
    else:
        sspec_f = scipy.ndimage.maximum_filter(sspec, size=5)
        I = np.where(np.logical_and(sspec_f==sspec, (sspec-np.mean(sspec))/np.std(sspec)>min_peak_lev))
        peaks = np.array([I[1], I[0]], dtype=float).T

    if peaks is not None:
        r_x = (vid_size[0]+vid_mrgn*2)/sspec.shape[1]
        r_y = (vid_size[1]+vid_mrgn*2)/sspec.shape[0]

        peaks[:,0] *= r_x
        peaks[:,1] *= r_y
        peaks -= vid_mrgn

    return peaks

def get_tau(data_dir, calibfile, speedOfSound, points=None, d=5, micpos=None, vid_mrgn=0):
    if micpos is None:
        with h5py.File(calibfile, mode='r') as f:
            micpos = f['/result/micpos'][()]
    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        r = f['/camera_param/depth_to_color_extrin/rotation'][()]
        t = f['/camera_param/depth_to_color_extrin/translation'][()]
        ppx = f['/camera_param/color_intrin/ppx'][()]
        ppy = f['/camera_param/color_intrin/ppy'][()]
        fx = f['/camera_param/color_intrin/fx'][()]
        fy = f['/camera_param/color_intrin/fy'][()]
        coeff = f['/camera_param/color_intrin/coeffs'][()]
        im_w = f['/camera_param/color_intrin/width'][()]
        im_h = f['/camera_param/color_intrin/height'][()]
        fs = f['/daq_param/fs'][()]
        camera_height = f['/camera_param/camera_height'][()]
        #if '/camera_param/camera_height' in f.keys():
        #    camera_height = f['/camera_param/camera_height'][()]
        #else:
        #    with open(config_path, 'r') as f_cfg:
        #        usvcam_cfg = yaml.safe_load(f_cfg)
        #        camera_height = usvcam_cfg['camera_height']

    coeff = coeff.ravel()

    r = np.reshape(r, [3, 3]).T
    cm = np.zeros([4,4])
    cm[:3, :3] = r
    cm[:3, 3] = t
    cm[3,3] = 1

    if points is None:
        im_u = np.arange(d/2 - vid_mrgn, im_w-1 + vid_mrgn, d)
        im_v = np.arange(d/2 - vid_mrgn, im_h-1 + vid_mrgn, d)
        U, V = np.meshgrid(im_u, im_v)
        U = np.ravel(U)
        V = np.ravel(V)
        grid_shape = (im_v.shape[0], im_u.shape[0])
    else:
        U = points[:,0]
        V = points[:,1]

    # undistort points
    mtx = [[fx, 0, ppx],
           [0, fy, ppy],
           [0, 0, 1]]
    mtx = np.array(mtx)
    dist = coeff[np.newaxis, :]
    p = np.array([U,V]).T
    p_undist = cv2.undistortPoints(p, mtx, dist)
    p_undist = p_undist[:,0,:]
    U = p_undist[:,0]
    V = p_undist[:,1]

    P = np.zeros([U.shape[0], 3])

    for i in range(U.shape[0]):
        u = U[i]
        v = V[i]
        A = np.vstack([u*cm[2,:] - cm[0, :], v*cm[2,:] - cm[1, :]])
        b = -(A[:,2]*camera_height + A[:,3])
        a = A[:, 0:2]
        x = np.matmul(np.linalg.inv(a), b)
        P[i,:] = np.array([x[0], x[1], camera_height])

    nMic = micpos.shape[0]

    TAU = np.zeros([P.shape[0], nMic])
    for i in range(nMic):
        mp0 = np.reshape(micpos[0,:], [1, 3])
        mp = np.reshape(micpos[i,:], [1, 3])
        D = np.sqrt(np.sum((P-mp)**2, axis=1)) - np.sqrt(np.sum((P-mp0)**2, axis=1))
        tau = D/speedOfSound
        TAU[:,i] = tau

    if points is None:
        return TAU, grid_shape
    else:
        return TAU

def is_sspec_localized(sspec, thr):
    return np.max(sspec) > np.mean(sspec) + thr * np.std(sspec)

def load_usvsegdata_ss(data_dir):
    L = glob.glob(data_dir + '/seg/*.ss.csv')
    seg = np.zeros((0,5))
    for i_file, fn in enumerate(L):
        if os.path.getsize(fn) == 0:
            continue
        a = np.genfromtxt(fn, delimiter=',')
        if a.ndim == 1:
            a = np.reshape(a, [1,-1])
        I = a[:,3]>0
        a = a[I,:]

        b = np.zeros((a.shape[0], 5))
        b[:,0:4] = a[:,0:4]
        b[:,4] = i_file+1

        seg = np.vstack((seg, b))
    
    return seg

def locate_all_segs(data_dir, calibfile, vid_mrgn=100, roi=None, out_sspec=False, color_eq=False, only_max=None, loc_thr=2.3, min_peak_lev=1.6): # z=2.3 means around 99% in cum dist
        
    if only_max is None:
        only_max = only_max_default

    #with open(config_path, 'r') as f:
    #    usvcam_cfg = yaml.safe_load(f)
    #speedOfSound = float(usvcam_cfg['speed_of_sound'])
    #pressure_calib = np.array(usvcam_cfg['sound_pressure_calibration'])

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        im_w = f['/camera_param/color_intrin/width'][()]
        im_h = f['/camera_param/color_intrin/height'][()]
        speedOfSound = f['/misc/speedOfSound'][()]
        pressure_calib = f['/misc/pressure_calib'][()]
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default

    vid_size = [im_w, im_h]

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    if out_sspec:
        os.makedirs(data_dir + '/loc', exist_ok=True)
        vid_file = data_dir + '/vid.mp4'
        vr = cv2.VideoCapture(vid_file)

    tau, grid_shape = get_tau(data_dir, calibfile, speedOfSound, d=5, vid_mrgn=vid_mrgn)

    seg = load_usvsegdata_ss(data_dir)
    if seg.shape[0] == 0:
        return

    _, I_ss = np.unique(seg[:,3:5], axis=0, return_inverse=True)

    rslt = np.zeros((np.max(I_ss)+1, 100+3))
    rslt[:,:] = np.nan
    with open_dat(data_dir) as fp_dat:
        for i_ss in tqdm(range(np.max(I_ss)+1)):

            seg2 = seg[I_ss==i_ss,:]

            i_frame = time2vidframe(data_dir, (np.min(seg2[:,0])+np.max(seg2[:,0]))/2, T, fs, cam_delay)
            
            S, peaks = locate_seg(fp_dat, seg2, tau, grid_shape, fs, n_ch, pressure_calib, vid_size, vid_mrgn=vid_mrgn, roi=roi, only_max=only_max, min_peak_lev=min_peak_lev)
            
            if out_sspec:
                vr.set(cv2.CAP_PROP_POS_FRAMES,i_frame)
                ret, frame = vr.read()
                img = draw_spec_on_vidframe(S, frame, vid_mrgn=vid_mrgn, color_eq=color_eq, only_max=only_max, min_peak_lev=min_peak_lev)
                if is_sspec_localized(S, loc_thr):
                    cv2.imwrite(data_dir + '/loc/{:06}-{:03}.jpg'.format(int(seg2[0,4]), int(seg2[0,3])), img)
                else:
                    cv2.imwrite(data_dir + '/loc/_{:06}-{:03}.jpg'.format(int(seg2[0,4]), int(seg2[0,3])), img)

            rslt[i_ss, 0] = seg2[0,4]
            rslt[i_ss, 1] = seg2[0,3]
            rslt[i_ss, 2] = i_frame
            
            if peaks is None:
                continue
            if peaks.shape[0] > 0:
                if is_sspec_localized(S, loc_thr):
                    p = peaks.ravel()
                    rslt[i_ss, 3:(3+p.shape[0])] = p

    np.savetxt(data_dir + '/loc.csv', rslt, delimiter=',',
                header=' ID_A,ID_B,video_frame,peak_locations (x1-y1-x2-y2-...),')

def locate_seg(fp_dat, seg, tau, grid_shape, fs, n_ch, pressure_calib, vid_size, vid_mrgn=0, roi=None, only_max=None, min_peak_lev=1.6):
        
    if only_max is None:
        only_max = only_max_default
        
    S, _ = calc_seg_power(fp_dat, seg, tau, fs, n_ch, pressure_calib, return_average=True)

    S = np.reshape(S, grid_shape)

    peaks = get_sspec_peaks(S, vid_size, vid_mrgn, roi, only_max=only_max, min_peak_lev=min_peak_lev)

    return S, peaks
     
def merge_assigned_segs(data_dir, n_mice, gap_min=0.03, conf_thr=0.99, min_dist=2.0):

    os.makedirs(data_dir + '/seg2', exist_ok=True)

    A = np.genfromtxt(data_dir + '/assign.csv', delimiter=',')

    I = np.logical_or(A[:,4] <= min_dist/180*np.pi, A[:,3] < conf_thr)
    A[I,1] = -1

    n_seg = int(np.max(A[:,0])+1)

    seg_id = np.zeros(n_seg)
    mouse_id = np.zeros(n_seg)
    seg_intv = np.zeros([n_seg,2])
    for i in range(n_seg):
        I = A[:,0] == i
        t = A[I,6]
        a = A[I,1]
        a = a[0]
        seg_id[i] = i
        mouse_id[i] = a
        seg_intv[i,0] = np.min(t)
        seg_intv[i,1] = np.max(t)

    data_summary = []

    cnt_seg = 0
    for i_mouse in range(-1,n_mice):
        if i_mouse < 0:
            I = np.where(mouse_id >= -10) # all
        else:
            I = np.where(np.logical_or(mouse_id==i_mouse, mouse_id<0))

        I = I[0]
        
        intv = seg_intv[I,:]
        intv[:,0] -= gap_min/2
        intv[:,1] += gap_min/2

        tree = intervaltree.IntervalTree.from_tuples(intv)
        tree.merge_overlaps(strict=False)
        tree = sorted(tree)
        tree = list(tree)

        mrg_intv = np.zeros((len(tree), 2))
        for i in range(len(tree)):
            ii = np.array(list(tree[i]))
            mrg_intv[i,:] = ii[0:2]

        grp = []
        for i in range(mrg_intv.shape[0]):
            grp.append([])
        for i in range(intv.shape[0]):
            c = np.mean(intv[i,:])
            ii = np.where(np.logical_and(mrg_intv[:,0] < c, mrg_intv[:,1] > c))
            ii = ii[0][0]
            grp[ii].append(I[i])

        for i in range(len(grp)):
            seg = []
            for j in range(len(grp[i])):
                seg.append(A[A[:, 0]==grp[i][j] , :])

            seg = np.concatenate(seg, axis=0)

            if np.max(seg[:,1]) <= -5:
                continue

            if i_mouse == -1:
                assign_rate = np.sum(seg[:,1]<0)/seg.shape[0]
            else:
                assign_rate = np.sum(seg[:,1]==i_mouse)/seg.shape[0]

            if assign_rate > 0 and i_mouse >= 0:
                mid = i_mouse+1
            elif assign_rate == 1 and i_mouse == -1:
                mid = -1
            else:
                continue

            np.savetxt(data_dir + '/seg2/seg{:06d}.csv'.format(cnt_seg), seg, delimiter=',',
                        header=' segment_id,assigned_mouse,--,confidence,snout_distance,p_value,time,freq,ampitude,ID_B,ID_A')

            t_start = np.min(seg[:,6])
            t_end = np.max(seg[:,6])
            dur = (t_end - t_start)*1000
            I_maxamp = np.argmax(seg[:,8])
            f_max = (seg[I_maxamp, 7])/1000
            a_max = seg[I_maxamp, 8]
            f_mean = (np.mean(seg[:, 7]))/1000
            f_cv = np.std(seg[:, 7])/(f_mean*1000)

            data_summary.append([cnt_seg, mid, assign_rate, t_start, t_end, dur, f_max, a_max, f_mean, f_cv])

            cnt_seg += 1

    data_summary = np.array(data_summary)
    if data_summary.shape[0] > 0:
        np.savetxt(data_dir + '/result.csv', data_summary, 
                    header='#,mouseID,assign_rate,start,end,duration,maxfreq,maxamp,meanfreq,cvfreq',
                    fmt='%d, %d, %.3f, %f, %f, %f, %f, %f, %f, %f', comments='')
    else:
        np.savetxt(data_dir + '/result.csv', data_summary, 
                    header='#,mouseID,assign_rate,start,end,duration,maxfreq,maxamp,meanfreq,cvfreq')

def pick_seg_for_calib(data_dir, n_pos_check=20):

    print('[How to use the GUI]')
    print('pick clear (noiseless, strong) vocal segments - space key, pick; other keys, skip...')

    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default

    T = np.genfromtxt(data_dir + '/sync.csv', delimiter=',')

    seg = load_usvsegdata_ss(data_dir)
    _, I_ss = np.unique(seg[:,4], axis=0, return_inverse=True)
    n_ss = np.max(I_ss)+1

    snout_pos = np.genfromtxt(data_dir + '/snout.csv', delimiter=',')
    snout_pos = snout_pos[:,:2]

    # snout pos of each call
    snout_pos_voc = np.zeros([n_ss, 2])
    for i_ss in range(n_ss):
        seg2 = seg[I_ss==i_ss,:]
        mrgn = 0.01
        t_intv = [np.min(seg2[:,0]) - mrgn, np.max(seg2[:,0]) + mrgn]
        i_frame = time2vidframe(data_dir, (np.min(seg2[:,0])+np.max(seg2[:,0]))/2, T, fs, cam_delay)
        snout_pos_voc[i_ss,:] = snout_pos[i_frame, :]

    I_notnan = np.logical_not(np.isnan(snout_pos_voc[:, 0]))

    kmeans_model = sklearn.cluster.KMeans(n_clusters=n_pos_check, random_state=10).fit(snout_pos_voc[I_notnan,:])
    labels = kmeans_model.labels_

    labels2 = -1*np.ones(n_ss)
    labels2[I_notnan] = labels
    labels = labels2

    L = glob.glob(data_dir + '/*.usvseg_dat.csv')
    wav_name = os.path.splitext(os.path.splitext(os.path.basename(L[0]))[0])[0]

    P = []
    SEG = []
    for i in range(n_pos_check):
        I = np.where(labels==i)
        I = I[0]

        for j in range(I.shape[0]):
            ii = I[j]
            seg2 = seg[I_ss==ii,:]
            imgfile = data_dir + '/seg/' + wav_name + '_{:04}.jpg'.format(int(seg2[0,4]))
            img = cv2.imread(imgfile)
            cv2.putText(img, 'pos{:02}'.format(i), (5, 20), cv2.FONT_HERSHEY_PLAIN, 1, (0, 0, 255), 1, cv2.LINE_AA)
            cv2.imshow('press_space_when_ok', img)
            k = cv2.waitKey()
            if k == 32: #Space key
                break
        
        P.append(snout_pos_voc[ii,:])
        SEG.append(seg2)

    cv2.destroyAllWindows()

    return SEG, P

def pick_seg_for_calib_manual(data_dir):
    print('[How to use the GUI]')
    print('CLICK=select sound source; SPACE(non-ESC key)=go to next sound; ESC=Quit')

    w_max = 600

    paramfile_path = data_dir + '/param.h5'
    vidfile_path = data_dir + '/vid.mp4'
    syncfile_path = data_dir + '/sync.csv'

    vr = cv2.VideoCapture(vidfile_path)

    with h5py.File(paramfile_path, mode='r') as f:
        fs = f['/daq_param/fs'][()]
        n_ch = f['/daq_param/n_ch'][()]
        if '/camera_param/cam_delay' in f:
            cam_delay = f['/camera_param/cam_delay'][()]
        else:
            cam_delay = cam_delay_default

    T = np.genfromtxt(syncfile_path, delimiter=',')

    seg = load_usvsegdata_ss(data_dir)
    _, I_ss = np.unique(seg[:, 4], axis=0, return_inverse=True)
    n_ss = np.max(I_ss) + 1

    # snout pos of each call
    P = []
    SEG = []
    L = glob.glob(data_dir + '/*.usvseg_dat.csv')
    wav_name = os.path.splitext(os.path.splitext(os.path.basename(L[0]))[0])[0]

    def update_disp():
        disp_img2 = copy.deepcopy(disp_img)
        if not np.isnan(p_crnt[0]):
            cv2.circle(disp_img2, (int(p_crnt[0]), int(p_crnt[1])), color=(0, 0, 255), thickness=2, radius=5)
        for p in P:
            cv2.drawMarker(disp_img2, (int(p[0]*r_scale), int(p[1]*r_scale)), color=(0, 255, 255), thickness=2, markerSize=5)

        cv2.imshow(wname, disp_img2)

    def onMouse(event, x, y, flag, params):
        if event == cv2.EVENT_LBUTTONDOWN:
            p_crnt[0] = x
            p_crnt[1] = y
            update_disp()

    wname = 'click sound source'
    cv2.namedWindow(wname)
    cv2.setMouseCallback(wname, onMouse, [])
    disp_img = np.zeros([int(vr.get(cv2.CAP_PROP_FRAME_HEIGHT)), int(vr.get(cv2.CAP_PROP_FRAME_WIDTH)), 3], np.uint8)
    
    if w_max < vr.get(cv2.CAP_PROP_FRAME_WIDTH):
        r_scale = w_max / vr.get(cv2.CAP_PROP_FRAME_WIDTH)
    else:
        r_scale = None
    
    p_crnt = np.zeros([2], float)
    p_crnt[:] = np.nan
    for i_ss in range(n_ss):
        seg2 = seg[I_ss == i_ss, :]
        i_frame = time2vidframe(data_dir, (np.min(seg2[:, 0]) + np.max(seg2[:, 0])) / 2, T, fs, cam_delay)
        vr.set(cv2.CAP_PROP_POS_FRAMES, i_frame)
        ret, frame = vr.read()

        imgfile = data_dir + '/seg/' + wav_name + '_{:04}.jpg'.format(int(seg2[0, 4]))
        img = cv2.imread(imgfile)

        disp_img = frame
        if r_scale is not None:
            disp_img = cv2.resize(disp_img, (int(disp_img.shape[1]*r_scale), int(disp_img.shape[0]*r_scale)))
        
        r = 100
        a = cv2.resize(img, (r, r))
        disp_img[-r:, -r:] = a
        update_disp()

        k = cv2.waitKey()
        if k == 27:  # ESC
            break

        if not np.isnan(p_crnt[0]):
            P.append(copy.deepcopy(p_crnt/r_scale))
            SEG.append(seg2)

        p_crnt[:] = np.nan

    cv2.destroyAllWindows()

    return SEG, P

def point2angle(data_dir, calibfile, points):

    with h5py.File(calibfile, mode='r') as f:
        micpos = f['/result/micpos'][()]
    paramfile = data_dir + '/param.h5'
    with h5py.File(paramfile, mode='r') as f:    
        r = f['/camera_param/depth_to_color_extrin/rotation'][()]
        t = f['/camera_param/depth_to_color_extrin/translation'][()]
        ppx = f['/camera_param/color_intrin/ppx'][()]
        ppy = f['/camera_param/color_intrin/ppy'][()]
        fx = f['/camera_param/color_intrin/fx'][()]
        fy = f['/camera_param/color_intrin/fy'][()]
        coeff = f['/camera_param/color_intrin/coeffs'][()]
        im_w = f['/camera_param/color_intrin/width'][()]
        im_h = f['/camera_param/color_intrin/height'][()]
        fs = f['/daq_param/fs'][()]
        camera_height = f['/camera_param/camera_height'][()]
        #if '/camera_param/camera_height' in f.keys():
        #    camera_height = f['/camera_param/camera_height'][()]
        #else:
        #    with open(config_path, 'r') as f_cfg:
        #        usvcam_cfg = yaml.safe_load(f_cfg)
        #        camera_height = usvcam_cfg['camera_height']

    r = np.reshape(r, [3, 3]).T
    cm = np.zeros([4,4])
    cm[:3, :3] = r
    cm[:3, 3] = t
    cm[3,3] = 1

    U = points[:,0]
    V = points[:,1]

    # undistort points
    mtx = [[fx, 0, ppx],
           [0, fy, ppy],
           [0, 0, 1]]
    mtx = np.array(mtx)
    dist = coeff[np.newaxis, :]
    p = np.array([U,V]).T
    p_undist = cv2.undistortPoints(p, mtx, dist)
    p_undist = p_undist[:,0,:]
    U = p_undist[:,0]
    V = p_undist[:,1]

    P = np.zeros([U.shape[0], 3])

    for i in range(U.shape[0]):
        u = U[i]
        v = V[i]
        A = np.vstack([u*cm[2,:] - cm[0, :], v*cm[2,:] - cm[1, :]])
        b = -(A[:,2]*camera_height + A[:,3])
        a = A[:, 0:2]
        x = np.matmul(np.linalg.inv(a), b)
        P[i,:] = np.array([x[0], x[1], camera_height])

    mp0 = np.reshape(micpos[0,:], [1, 3])
    P = P-mp0

    x = P[:,2]
    y = P[:,0]
    z = P[:,1]

    azimuth = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))

    return azimuth, elevation

def rad2deg(th):
    return th/np.pi*180

def time2vidframe(data_dir, t, T, fs, cam_delay):
    # T: sync.csv

    TT = T[:,1]/fs
    i_frame = np.argmin(abs(TT-(t+cam_delay)))

    return i_frame

def to_intervals(f,t):
    if f[-1] == True:
        f = np.concatenate([f, np.zeros([1], dtype=bool)])
    f = np.concatenate([np.zeros([1], dtype=bool), f])
    d = np.diff(f.astype(int))
    intv = np.array([t[np.where(d==1)[0]], t[np.where(d==-1)[0]-1]])
    return intv.T

def save_intrinsics(f, group_name, intrin):
    f.create_dataset(group_name + '/coeffs', data = intrin.coeffs)
    f.create_dataset(group_name + '/fx', data = intrin.fx)
    f.create_dataset(group_name + '/fy', data = intrin.fy)
    f.create_dataset(group_name + '/width', data = intrin.width)
    f.create_dataset(group_name + '/height', data = intrin.height)
    #f.create_dataset(group_name + '/model', data = intrin.model)
    f.create_dataset(group_name + '/ppx', data = intrin.ppx)
    f.create_dataset(group_name + '/ppy', data = intrin.ppy)

def save_extrinsics(f, group_name, extrin):
    f.create_dataset(group_name + '/rotation', data = extrin.rotation)
    f.create_dataset(group_name + '/translation', data = extrin.translation)

def update_paramfile(data_dir):
    paramfile = data_dir + '/param.h5'

    with open(config_path, 'r') as f:
        usvcam_cfg = yaml.safe_load(f)
    speedOfSound = float(usvcam_cfg['speed_of_sound'])
    pressure_calib = np.array(usvcam_cfg['sound_pressure_calibration'])
    mic0pos = np.asarray(usvcam_cfg['mic0pos'])/1000.0

    with h5py.File(paramfile, mode='r+') as f:  
        if '/misc/speedOfSound' in f.keys():
            del f['/misc/speedOfSound']
            
        if '/misc/pressure_calib' in f.keys():
            del f['/misc/pressure_calib']

        if '/misc/mic0pos' in f.keys():
            del f['/misc/mic0pos']

        f.create_dataset('/misc/speedOfSound', data=speedOfSound)
        f.create_dataset('/misc/pressure_calib', data=pressure_calib)
        f.create_dataset('/misc/mic0pos', data=mic0pos)


