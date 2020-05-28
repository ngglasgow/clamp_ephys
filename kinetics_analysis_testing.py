import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import numpy as np

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
width = 3    

# paths to clamp_ephys/test_data/
project_path = os.path.join(os.path.expanduser('~'), 'Documents', 'phd_projects', 'Injected_GC_data')
notes_path = os.path.join(project_path, 'VC_pairs', 'tables', 'p2_data_notes.csv')
cell_path = os.path.join(project_path, 'VC_pairs', 'data', 'p2', 'JH200313_c3_light100.ibw')

data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=notes_path, timepoint='p2', amp_factor=amp_factor)

data.get_raw_peaks(stim_time, post_stim, polarity='-', baseline_start=baseline_start, baseline_end=baseline_end)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)
data.sweeps = data.traces - data.new_baseline_raw

data.get_peaks_widths(stim_time, fs, width)
   
data.get_peaks_kinetics(stim_time, fs)

