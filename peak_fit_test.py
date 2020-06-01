import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
#%matplotlib
import matplotlib.pyplot as plt
import numpy as np

lowpass_freq = 500  # Hz
stim_time = 500     # ms
post_stim = 250     # ms, amount of time after stimulus to look for max value
baseline_start = 3000   # ms, time in the sweep to start taking the baseline
baseline_end = 6000     # ms, time in the sweep at which baseline ends
tp_start = 50        # ms, time of start of test pulse
vm_jump = -5        # mV, test pulse voltage jump
pre_tp = 11          # ms, amount of time before test pulse start to get baseline
unit_scaler = -12   # unitless, scaler to get back to A, from pA
amp_factor = 1000      # scaler for making plots in pA
fs = 10             # kHz, the sampling frequency
width = 3               # ms, the required width for classifying an event as a peak
timepoint = 'p14'

# paths to clamp_ephys/test_data/
project_path = os.path.join(os.path.expanduser('~'), 'Documents', 'phd_projects', 'Injected_GC_data')
notes_path = os.path.join(project_path, 'VC_pairs', 'tables', 'p14_data_notes.csv')
cell_path = os.path.join(project_path, 'VC_pairs', 'data', 'p14', 'JH190828_c6_light100_1.ibw')

data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=notes_path, timepoint=timepoint, amp_factor=amp_factor)

data.get_raw_peaks(stim_time, post_stim, polarity='-', baseline_start=baseline_start, baseline_end=baseline_end)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)

data.get_peaks_widths(stim_time, width)   

sweeps = data.traces_filtered - data.new_baseline_raw
window_start = (stim_time + 20) * fs
all_peaks_kinetics_df = pd.DataFrame()
all_peaks_kinetics_avg_df = pd.DataFrame()
first3_kinetics_avg_df = pd.DataFrame()



#%%
sweep = 3
peak_number = 11     
trace = sweeps.iloc[window_start:, sweep].values
thresh = 2.5 * trace.std()

peaks, properties = scipy.signal.find_peaks(trace * -1, distance=0.5*fs, prominence=thresh, width=width*fs)
prominence_data = list(properties.values())[0:3]

fig = plt.figure()
fig.suptitle('Sweep {}'.format(sweep))
plt.plot(trace)
plt.plot(peaks, trace[peaks], 'x')

if len(peaks) > 0:
    ninety_left = data.all_widths_df.loc[(sweep), 'ninety_left']
    ten_left = data.all_widths_df.loc[(sweep), 'ten_left']
    ten_to_ninety = (ninety_left - ten_left) / fs
    peak_time = data.all_widths_df.loc[(sweep), 'peaks_index'] / fs
    hw_time = data.all_widths_df.loc[(sweep), 'half_widths'] / fs

    if data.mean_peak_filtered > 0:
        invert = 1

    else:
        invert = -1

    peak = peaks[peak_number]

    tau = data.get_tau(trace, sweep, peak_number)

else:
    print('No peaks in sweep {}'.format(sweep))


# %%
