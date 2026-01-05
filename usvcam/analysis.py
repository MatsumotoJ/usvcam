import usvseg
import glob
import os
import yaml
import shutil

from . import tool

def enable_gpu():
    tool.enable_gpu()

def disable_gpu():
    tool.disable_gpu()

def set_default_analysis_params_for_device(dev_name):
    tool.set_default_analysis_params_for_device(dev_name)

def assign_vocalizations(data_dir, calibfile, assignfile, n_mice, conf_thr=0.99, gap_min=0.03, only_max=None, loc_thr=2.3, min_peak_lev=1.6, min_dist=2.0):
    
    print("cleaning data directory...")
    tool.clean_data_dir(data_dir, filekeep=['vid.loc.mp4', 'vmstat.*.h5'])

    if len(glob.glob(data_dir + '/seg/*.ss.csv')) == 0:
        print('--NO USV SEGMENT IN THE DATA')
        return

    print('localize segments...')
    tool.locate_all_segs(data_dir, calibfile, vid_mrgn=100, only_max=only_max, loc_thr=loc_thr, min_peak_lev=min_peak_lev)

    if not os.path.exists(data_dir + '/loc.csv'):
        print('--NO LOCALIZED USV SEGMENT ')
        return

    print('assign segments...')
    tool.assign_all_segs(data_dir, calibfile, assignfile, n_mice, conf_thr=conf_thr)

    print('merge assigned segments...')
    tool.merge_assigned_segs(data_dir, n_mice, gap_min=gap_min, min_dist=min_dist)
    
    print('done.')

def calib_micpos(data_dir, mannual_mode=False, n_sample=20, outpath=None, pos_lim=None, pos_searchlim=None, pos_init=None, vis_progress=False, calibfile=None, ftol=1e-3, maxiter=100):
    if mannual_mode:
        SEG, P = tool.pick_seg_for_calib_manual(data_dir)
    else:
        SEG, P = tool.pick_seg_for_calib(data_dir, n_sample)
    tool.calib_micpos(data_dir, SEG, P, calibfile=calibfile, h5f_outpath=outpath, pos_lim=pos_lim, pos_searchlim=pos_searchlim, pos_init=pos_init, vis_progress=vis_progress, ftol=ftol, maxiter=maxiter)
    print('done.')

def estimate_assign_param(data_dirs, calibfiles, assignfile, n_iter=8, n_trial=7, show_figs=False, only_max=None, loc_thr=2.3, min_peak_lev=1.6):
    
    n_data = len(data_dirs)
    
    print("cleaning data directories...")
    for i_data in range(n_data):
        tool.clean_data_dir(data_dirs[i_data])

    print("localizing semgents...")
    for i_data in range(n_data):
        tool.locate_all_segs(data_dirs[i_data], calibfiles[i_data], only_max=only_max, loc_thr=loc_thr, min_peak_lev=min_peak_lev)

    print("calculate stats with virtual mice... (This may take hours)")
    for i_data in range(n_data):
        print('- processing {}/{}...'.format(i_data+1, n_data), flush=True)
        for iter_id in range(n_iter):
            tool.calc_vm_stats(data_dirs[i_data], calibfiles[i_data], roi=None, iter_id=iter_id)

    print("estimating parameters for assignment...")
    tool.estimate_assign_param(data_dirs, calibfiles, assignfile, n_trial=n_trial, iter_ids=list(range(n_iter)), show_figs=show_figs)

def create_assignment_video(data_dir, n_mice, color_eq=False, min_dist=2.0):
    if not os.path.exists(data_dir + '/assign.csv'):
        print('--NO ASSIGNED USV SEGMENT ')
        return
    tool.create_assignment_video(data_dir, n_mice, color_eq=color_eq, min_dist=min_dist)

def create_localization_video(data_dir, calibfile, t_end=-1, color_eq=False, only_max=None, min_peak_lev=1.6):
    tool.create_localization_video(data_dir, calibfile, t_end=t_end, color_eq=color_eq, only_max=only_max, min_peak_lev=min_peak_lev)

def dat2wav(data_dir, i_ch, offset=0):
    tool.dat2wav(data_dir, i_ch, offset)

def extract_ch(data_dir, i_ch, offset=0):

    if os.path.isfile(data_dir + '/snd.dat'):
        tool.dat2flac(data_dir, i_ch, offset)
    elif os.path.isfile(data_dir + '/snd.flac'):
        tool.flac2flac(data_dir, i_ch, offset)

def run_usvseg(data_dir, usvseg_prm_file, i_ch=0, offset=0):

    outp = data_dir + '/seg'

    # clean files to avoid contamination with the previous calculation results
    if os.path.exists(outp):
        shutil.rmtree(outp)
    L = glob.glob(data_dir + '/*.audible.*')
    for l in L:
        os.remove(l)
    L = glob.glob(data_dir + '/*.ch*.wav')
    for l in L:
        os.remove(l)
    L = glob.glob(data_dir + '/*.ch*.flac')
    for l in L:
        os.remove(l)
    L = glob.glob(data_dir + '/*.usvseg_dat.csv')
    for l in L:
        os.remove(l)

    # extrac ch for usvseg
    if i_ch == 0 and os.path.isfile(data_dir + '/snd.flac'):
        fp = data_dir + '/snd.flac'
    else: 
        extract_ch(data_dir, i_ch, offset)
        L = glob.glob(data_dir + '/snd.ch*.flac')
        fp = L[0]

    # run usvseg
    with open(usvseg_prm_file, 'r') as f:
        params = yaml.load(f, Loader=yaml.SafeLoader)
    savefp = os.path.splitext(fp)[0] + '.usvseg_dat.csv'
    fname_audiblewav = os.path.splitext(fp)[0] + '.audible.flac'

    usvseg.proc_wavfile(params, fp, savefp, outp, fname_audiblewav=fname_audiblewav, usvcamflg=True)