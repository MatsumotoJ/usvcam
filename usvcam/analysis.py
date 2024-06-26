import usvseg
import glob
import os
import yaml

from . import tool

def assign_vocalizations(data_dir, calibfile, assignfile, n_mice, conf_thr=0.99, gap_min=0.03):
    
    print("cleaning data directory...")
    tool.clean_data_dir(data_dir, filekeep=['vid.loc.mp4', 'vmstat.*.h5'])

    if len(glob.glob(data_dir + '/seg/*.ss.csv')) == 0:
        print('--NO USV SEGMENT IN THE DATA')
        return

    print('localize segments...')
    tool.locate_all_segs(data_dir, calibfile, vid_mrgn=100)

    if not os.path.exists(data_dir + '/loc.csv'):
        print('--NO LOCALIZED USV SEGMENT ')
        return

    print('assign segments...')
    tool.assign_all_segs(data_dir, calibfile, assignfile, n_mice, conf_thr=conf_thr)

    print('merge assigned segments...')
    tool.merge_assigned_segs(data_dir, n_mice, gap_min=gap_min)
    
    print('done.')

def calib_with_voc(data_dir, mannual_mode=False, n_sample=20, outpath=None, pos_lim=None, pos_init=None, vis_progress=False):
    if mannual_mode:
        SEG, P = tool.pick_seg_for_calib_manual(data_dir)
    else:
        SEG, P = tool.pick_seg_for_calib(data_dir, n_sample)
    tool.calc_micpos_with_voc(data_dir, SEG, P, h5f_outpath=outpath, pos_lim=pos_lim, pos_init=pos_init, vis_progress=vis_progress)
    print('done.')

def estimate_assign_param(data_dirs, calibfiles, assignfile, n_iter=8, n_trial=7, show_figs=False):
    
    n_data = len(data_dirs)
    
    print("cleaning data directories...")
    for i_data in range(n_data):
        tool.clean_data_dir(data_dirs[i_data])

    print("localizing semgents...")
    for i_data in range(n_data):
        tool.locate_all_segs(data_dirs[i_data], calibfiles[i_data])

    print("calculate stats with virtual mice... (This may take hours)")
    for i_data in range(n_data):
        print('- processing {}/{}...'.format(i_data+1, n_data), flush=True)
        for iter_id in range(n_iter):
            tool.calc_vm_stats(data_dirs[i_data], calibfiles[i_data], roi=None, iter_id=iter_id)

    print("estimating parameters for assignment...")
    tool.estimate_assign_param(data_dirs, calibfiles, assignfile, n_trial=n_trial, iter_ids=list(range(n_iter)), show_figs=show_figs)

def create_assignment_video(data_dir, n_mice, color_eq=False):
    if not os.path.exists(data_dir + '/assign.csv'):
        print('--NO ASSIGNED USV SEGMENT ')
        return
    tool.create_assignment_video(data_dir, n_mice, color_eq=color_eq)

def create_localization_video(data_dir, calibfile, t_end=-1, color_eq=False):
    tool.create_localization_video(data_dir, calibfile, t_end=t_end, color_eq=color_eq)

def dat2wav(data_dir, i_ch, offset=0):
    tool.dat2wav(data_dir, i_ch, offset)

def run_usvseg(data_dir, usvseg_prm_file):

    fp = glob.glob(data_dir + '/*.wav')
    fp = fp[0]

    with open(usvseg_prm_file, 'r') as f:
        params = yaml.load(f, Loader=yaml.SafeLoader)

    savefp = os.path.splitext(fp)[0] + '.usvseg_dat.csv'
    outp = data_dir + '/seg'
    fname_audiblewav = os.path.splitext(fp)[0] + '.audible.wav'

    usvseg.proc_wavfile(params, fp, savefp, outp, fname_audiblewav=fname_audiblewav, usvcamflg=True)