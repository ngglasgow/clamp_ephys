import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
#%matplotlib
import matplotlib.pyplot as plt
import numpy as np

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
paths = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs')
tables = paths.tables
figures = paths.figures

'''####################### SET PROPER PATH_TO_DATA_NOTES ########################## '''
data_path = r"C:\Users\jhuang\Documents\phd_projects\Injected_GC_data\VC_pairs\tables\\p2_data_notes.csv"

''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
lowpass_freq = 500      # Hz
stim_time = 500         # ms
post_stim = 250         # ms, amount of time after stimulus to look for max value
baseline_start = 3000   # ms, time in the sweep to start taking the baseline
baseline_end = 6000     # ms, time in the sweep at which baseline ends
tp_start = 5            # ms, time of start of test pulse
vm_jump = 10            # mV, test pulse voltage jump
pre_tp = 3              # ms, amount of time before test pulse start to get baseline
unit_scaler = -12       # unitless, scaler to get back to A, from pA
amp_factor = 1          # scaler for making plots in pA
fs = 25                 # kHz, the sampling frequency

'''#################### THIS LOOP RUNS THE SINGLE CELL ANALYSIS #################### '''
#cell_path = os.path.join(r'C:\Users\jhuang\Documents\phd_projects\Injected_GC_data\VC_pairs\data\p2', 'JH200304_c2_light100.ibw')

for path in paths.p2_paths:

    data = clamp_ephys.workflows.cell(path, fs=fs, path_to_data_notes=data_path, timepoint='p2', amp_factor=amp_factor)

    data.get_raw_peaks(stim_time, post_stim)
    data.filter_traces(lowpass_freq)
    data.get_filtered_peaks(stim_time, post_stim)

    '''===================================================================================
    testing for kinetics
    - notes on find_peaks, peak_prominence, and peak_widths
        - use std of entire signal as a thresholding - std of baseline seems noisier
        - should std_baseline allow pre_stim as well? don't remember why it doesn't...
        - heights are being calculated as simply above 0, so just x, not that useful
        - prominence is acconting for lowest contour, so inherently corrects for baseline (very immediate baseline)
        - any peaks MUST be positive, so need to invert for sure
        - 
    '''
    peaks = data.peaks_filtered_indices
    subtracted_data = data.traces_filtered - data.baseline_filtered

    window_start = (stim_time + 20) * fs
    baseline_start = 3000 * fs

    # this below would be outside of the per-cell loop

    for sweep in range(len(subtracted_data.columns)):
        
        # using non-subtracted data so I can see which sweeps to drop
        # window omits TP and pre-stimulus time
        x = data.traces_filtered.iloc[window_start:, sweep].values
        baseline = data.traces_filtered.iloc[baseline_start:, sweep].values
        thresh = 3 * baseline.std()
        sweep_length = len(x)
        sweep_time = np.arange(0, sweep_length/fs, 1/fs)

        # finding all peaks
        peaks, properties = scipy.signal.find_peaks(x * -1, prominence=thresh)
        prominence_data = tuple(properties.values())

        # correct peaks time for fs
        peaks_corr = peaks/fs

        fig = plt.figure()
        fig.suptitle('Sweep {}'.format(sweep))
        plt.plot(sweep_time, x)
        plt.plot(peaks_corr, x[peaks], 'x')

        filename = '{}_sweep_{}.png'.format(file_id, sweep)
        fig_path = os.path.join(r'C:\Users\jhuang\Documents\phd_projects\Injected_GC_data\VC_pairs\figures\p2\sweeps_incl_depol', path, filename)
        metadata.check_create_dirs(fig_path)
        fig.savefig(fig_path, dpi=300, format='png')

