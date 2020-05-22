import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
#%matplotlib
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
project_path = os.path.join(os.path.expanduser('~'), 'Documents', 'phd_projects', 'git_repos','clamp_ephys')
notes_path = os.path.join(project_path, 'test_data', 'p2_data_notes.csv')
cell_path = os.path.join(project_path, 'test_data', 'JH200303_c1_light100.ibw')

data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=notes_path, timepoint='p2', amp_factor=amp_factor)

data.get_raw_peaks(stim_time, post_stim, polarity='-', baseline_start=baseline_start, baseline_end=baseline_end)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)

sweeps = data.traces - data.new_baseline_raw
window_start = (stim_time + 20) * fs

trace = sweeps.iloc[window_start:, 0].values
thresh = 2.5 * trace.std()

peaks, properties = scipy.signal.find_peaks(trace * -1, distance=0.5*fs, prominence=thresh, width=width*fs)
prominence_data = list(properties.values())[0:3]
# fig = plt.figure()
# fig.suptitle('Sweep {}'.format(sweep))
# plt.plot(trace)
# plt.plot(peaks, trace[peaks], 'x')

# calculate 10 to 90% and FWHM
ten_widths, ten_height, ten_left, ten_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=0.9, prominence_data=prominence_data)
ninety_widths, ninety_height, ninety_left, ninety_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=0.1, prominence_data=prominence_data)
half_widths, hw_height, hw_left, hw_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=0.5, prominence_data=prominence_data)
full_widths, fw_height, fw_left, fw_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=1, prominence_data=prominence_data)

hw_height = hw_height * -1
fw_height = fw_height * -1
ten_height = ten_height * -1
ninety_height = ninety_height * -1

ten_to_ninety = (ninety_left - ten_left) / fs
prominences = properties['prominences']
hw_time = half_widths / fs

peak = peaks[1]

# make sure to pull this out of the loops
def decay_func(time, current_peak, tau, offset):
    return current_peak * np.exp(-time/tau) + offset

decay_end = fw_right[1].astype(int) # indexing needs to be a whole number
peak_trace = trace[peak:decay_end + 1] * -1
xtime_adj = np.arange(0, len(peak_trace) / fs, 1 / fs)

decay_end_ninety = ten_right[1].astype(int) # indexing needs to be a whole number
peak_trace_ninety = trace[peak:decay_end_ninety + 1] * -1
xtime_adj_ninety = np.arange(0, len(peak_trace_ninety) / fs, 1 / fs)

# if you actually need guesses, i'm taking a guess at tau
# take first index value that goes below 1*tau value, then multiply by 2 to account for noise
# and divide by sampling time to get ms time scale
# these should be the same regardless of where stopping
guess_tau_time  = np.where(peak_trace < (peak_trace[1] * 0.37))[0][0] * 2 / fs
starting_params = [peak_trace[1], guess_tau_time, 0]

# fits
popt, pcov = scipy.optimize.curve_fit(f=decay_func, xdata=xtime_adj, ydata=peak_trace)
current_peak, tau, offset = popt

popt_ninety, pcov_ninety = scipy.optimize.curve_fit(f=decay_func, xdata=xtime_adj_ninety, ydata=peak_trace_ninety)
current_peak_ninety, tau_ninety, offset_ninety = popt_ninety

plt.figure()
plt.plot(xtime_adj, peak_trace, color='k', label='data')
plt.plot(xtime_adj, decay_func(xtime_adj, *popt), color='r', label=f'fit on FW; tau (ms): {round(tau, 2)}')
plt.plot(xtime_adj_ninety, decay_func(xtime_adj_ninety, *popt_ninety), color='b', label=f'fit on 90%; tau (ms): {round(tau_ninety, 2)}')
plt.legend()

