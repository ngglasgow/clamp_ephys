import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
#%matplotlib
import matplotlib.pyplot as plt
import numpy as np

''' This file loops through each .ibw file in the p2 directory and plots all individual 
sweeps for each cell.
'''

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
paths = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs')
tables = paths.tables
figures = paths.figures

'''####################### SET PROPER PATH_TO_DATA_NOTES ########################## '''
data_path = '/home/jhuang/Documents/phd_projects/Injected_GC_data/VC_pairs/tables/p2_data_notes.csv'

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

path_to_figures = os.path.join(figures, 'p2', 'sweeps_wdepol')

for path in paths.p2_paths:

    data = clamp_ephys.workflows.cell(path, fs=fs, path_to_data_notes=data_path, timepoint='p2', amp_factor=amp_factor)

    data.get_raw_peaks(stim_time, post_stim)
    data.filter_traces(lowpass_freq)
    data.get_filtered_peaks(stim_time, post_stim)

    peaks = data.peaks_filtered_indices
    
    # plotting filtered, unsubtracted traces
    sweeps = data.traces_filtered

    for sweep in range(len(sweeps.columns)):

        data.plot_sweeps(sweep, stim_time, baseline_start, save_fig=True, path_to_figures=path_to_figures)