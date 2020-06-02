import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
#%matplotlib
import matplotlib.pyplot as plt
import numpy as np


# test file paths that are part of the repo with relative paths
test_data_path = 'test_data'
test_data_files = [os.path.join(test_data_path, file) for file in os.listdir(test_data_path) if file.endswith('ibw')]
test_notes = os.path.join(test_data_path, 'p14_data_notes.csv')

''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
lowpass_freq = 500      # Hz
stim_time = 500         # ms
post_stim = 250         # ms, amount of time after stimulus to look for max value
baseline_start = 3000   # ms, time in the sweep to start taking the baseline
baseline_end = 6000     # ms, time in the sweep at which baseline ends
tp_start = 5            # ms, time of start of test pulse
vm_jump = -5            # mV, test pulse voltage jump
pre_tp = 11              # ms, amount of time before test pulse start to get baseline
unit_scaler = -12       # unitless, scaler to get back to A, from pA
amp_factor = 1000          # scaler for making plots in pA
fs = 10                 # kHz, the sampling frequency
width = 3               # ms, the required width for classifying an event as a peak
timepoint = 'p14'

test_cell = 'test_data/JH190828_c6_light100_1.ibw'
#test_cell = test_data_files[0]
data = clamp_ephys.workflows.cell(test_cell, fs=fs, path_to_data_notes=test_notes, timepoint=timepoint, amp_factor=amp_factor)

data.get_raw_peaks(stim_time, post_stim, polarity='-', baseline_start=baseline_start, baseline_end=baseline_end)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)
data.get_series_resistance(tp_start, vm_jump, pre_tp, unit_scaler)
data.get_sweep_data()
data.get_responses(threshold=5)
data.get_sweepavg_summary()

data.get_peaks_widths(stim_time, width)
data.get_peaks_kinetics(stim_time)
