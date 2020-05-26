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
sweep = 18

data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=notes_path, timepoint='p2', amp_factor=amp_factor)

data.get_raw_peaks(stim_time, post_stim, polarity='-', baseline_start=baseline_start, baseline_end=baseline_end)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)

sweeps = data.traces - data.new_baseline_raw
window_start = (stim_time + 20) * fs

trace = sweeps.iloc[window_start:, sweep].values
thresh = 2.5 * trace.std()

peaks, properties = scipy.signal.find_peaks(trace * -1, distance=0.5*fs, prominence=thresh, width=width*fs)
prominence_data = list(properties.values())[0:3]
# fig = plt.figure()
# fig.suptitle('Sweep {}'.format(sweep))
# plt.plot(trace)
# plt.plot(peaks, trace[peaks], 'x')

if len(peaks) > 0:

    # calculate 10 to 90% and FWHM
    ten_widths, ten_height, ten_left, ten_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=0.9, prominence_data=prominence_data)
    ninety_widths, ninety_height, ninety_left, ninety_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=0.1, prominence_data=prominence_data)
    half_widths, hw_height, hw_left, hw_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=0.5, prominence_data=prominence_data)
    full_widths, fw_height, fw_left, fw_right = scipy.signal.peak_widths(trace * -1, peaks, rel_height=1, prominence_data=prominence_data)
    results_full = scipy.signal.peak_widths(trace * -1, peaks, rel_height=1, prominence_data=prominence_data)

    hw_height = hw_height * -1
    fw_height = fw_height * -1
    ten_height = ten_height * -1
    ninety_height = ninety_height * -1


    plt.plot(trace)
    plt.plot(peaks, trace[peaks], "x")
    
    plt.hlines(ten_height, ten_left, ten_right, color="C3")
    plt.show()

    ten_to_ninety = (ninety_left - ten_left) / fs
    prominences = properties['prominences']
    hw_time = half_widths / fs
    fw_time = full_widths / fs

    # using 10% to 10% as trace for now
    trace_start = ten_left[8].astype(int)
    trace_end = ten_right[8].astype(int) # indexing needs to be a whole number
    whole_trace = trace[trace_start:trace_end+1] * -1
    xtime_adj = np.arange(0, len(whole_trace) / fs, 1 / fs)


    # take baseline as 10 ms before event onset
    baseline_window_start = trace_start - 10 * fs
    event_baseline = np.mean(trace[baseline_window_start:trace_start]) * -1
    whole_trace_sub = whole_trace - event_baseline

    # # test on noise
    # noise = np.random.normal(0,1,100000000)
    # test_baseline = np.mean(noise)
    # noise_sub = noise - test_baseline
    # test_x = np.arange(0, len(noise))
    # test_integral = scipy.integrate.trapz(noise_sub, test_x)

    
    charge = scipy.integrate.trapz(whole_trace_sub, xtime_adj)

else:
    print('No peaks in sweep {}'.format(sweep))

