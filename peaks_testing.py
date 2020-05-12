import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import numpy as np

'''####################### SET PROPER PATH_TO_DATA_NOTES ########################## '''
data_path = os.path.join(os.getcwd(), 'test_data', 'p2_data_notes.csv')

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
cell_path = os.path.join(os.getcwd(), 'test_data', 'JH200311_c1_light100.ibw')
data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=data_path, timepoint='p2', amp_factor=amp_factor)

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

# finds all peaks, and prominences for an entire trace
test_x = subtracted_data.iloc[1000:, 0].values
non_sub_x = data.traces_filtered.iloc[1000:, 0].values
thresh = test_x.std()
thresh_non_sub = non_sub_x.std()

peaks, properties = scipy.signal.find_peaks(non_sub_x * -1, prominence=thresh)
prominence_data = tuple(properties.values())
plt.figure()
plt.plot(non_sub_x)
plt.plot(peaks, non_sub_x[peaks], 'x')

half_widths, hw_height, hw_left, hw_right = scipy.signal.peak_widths(non_sub_x * -1, peaks, rel_height=0.5, prominence_data=prominence_data)
full_widths, fw_height, fw_left, fw_right = scipy.signal.peak_widths(non_sub_x * -1, peaks, rel_height=1, prominence_data=prominence_data)
hw_height = hw_height * -1
fw_height = fw_height * -1

plt.figure()
plt.plot(non_sub_x)
plt.plot(peaks, non_sub_x[peaks], 'x')
plt.hlines(hw_height, xmin=hw_left, xmax=hw_right, color='C2')
plt.hlines(fw_height, xmin=fw_left, xmax=fw_right, color='C3')

# find half width of just a defined peak
test_peak = [peaks[2]]
test_hw, test_hw_height, test_hw_left, test_hw_right = scipy.signal.peak_widths(non_sub_x * -1, test_peak, rel_height=0.5)
test_hw_height = test_hw_height * -1

plt.figure()
plt.plot()
plt.plot(non_sub_x)
plt.plot(test_peak, non_sub_x[test_peak], 'x')
plt.hlines(test_hw_height, test_hw_left, test_hw_right, color='C2')

# find the half width of a barely respond defined peak in window
x = data.traces_filtered[0].values
peak = [data.peaks_filtered_indices[0]]

half_widths = scipy.signal.peak_widths(x * -1, peak, rel_height=0.5)
plt.figure()
plt.plot()
plt.plot(x * -1)
plt.plot(peak, x[peak] * -1, 'x')
plt.hlines(*half_widths[1:], color='C2')
hw = half_widths[0]

hw_df = pd.DataFrame()
for i in range(len(data.traces_filtered.columns)):
    x = data.traces_filtered[i].values * -1
    peak = [data.peaks_filtered_indices[i]]
    hw, hw_height, hw_left, hw_right = scipy.signal.peak_widths(x, peak, rel_height=0.5)
    hw_time = hw / fs
    hw_df = pd.concat([hw_df, pd.DataFrame(hw_time)], ignore_index=True)
hw_df.columns = ['Max peak half-width (ms)']

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# test new class function
data.get_max_peak_half_width()


# drop first 520 s from sweep to account for TP and only look at time from stimulus onset to end of sweep
window_start = (stim_time + 20) * fs

# this below would be outside of the per-cell loop
p2_kinetics_summary = pd.DataFrame()

columns_index = ['peak_index', 'prominence', '10 to 90% RT (ms)', 'half-width (ms)']
first_peaks_properties_df = pd.DataFrame()

for sweep in range(len(subtracted_data.columns)):
    
    x = subtracted_data.iloc[window_start:, sweep].values
    thresh = 3 * x.std()

    # finding all peaks
    peaks, properties = scipy.signal.find_peaks(x * -1, prominence=thresh)
    prominence_data = tuple(properties.values())
    plt.figure()
    plt.plot(x)
    plt.plot(peaks, x[peaks], 'x')

    '''
    ten_widths, ten_height, ten_left, ten_right = scipy.signal.peak_widths(x * -1, peaks, rel_height=0.9, prominence_data=prominence_data)
    ninety_widths, ninety_height, ninety_left, ninety_right = scipy.signal.peak_widths(x * -1, peaks, rel_height=0.1, prominence_data=prominence_data)
    half_widths, hw_height, hw_left, hw_right = scipy.signal.peak_widths(x * -1, peaks, rel_height=0.5, prominence_data=prominence_data)
    full_widths, fw_height, fw_left, fw_right = scipy.signal.peak_widths(x * -1, peaks, rel_height=1, prominence_data=prominence_data)

    hw_height = hw_height * -1
    fw_height = fw_height * -1
    ten_height = ten_height * -1
    ninety_height = ninety_height * -1

    ten_to_ninety = (ninety_left - ten_left) / fs
    prominences = properties['prominences']
    hw_time = half_widths / fs

    peaks_array = np.array((peaks, prominences, ten_to_ninety, hw_time)).T
    peaks_data = pd.DataFrame(peaks_array, columns=columns_index)

    # extracting kinetics of the first response/peak
    first_peak_data = peaks_data.iloc[0]
    first_peaks_properties_df = pd.concat([first_peaks_properties_df, first_peak_data], axis=1, ignore_index=True).T
    '''
    plt.plot(x)

first_peaks_properties_avg = first_peaks_properties_df.mean
p2_kinetics_summary = pd.concat([p2_kinetics_summary, first_peaks_properties_avg], ignore_index=True)

'''
# finding latency to response, i.e. the time at which current exceeds 3x std of baseline
baseline_start = baseline_start * fs
baseline_end = baseline_end * fs
baseline = new_mean_baseline(data, fs, baseline_start, baseline_end)
baseline_std = new_std_baseline(data, fs, baseline_start, baseline_end)
'''

# find 10-90% rise time
# asking to eval at rel_height equivalent to evaluating at height of Peak(height) - prominence * rel_height
# so evaluating at rel_height = 0.1 does not equal evaluating at 10% of Peak(height)
test_hw, test_hw_height, test_hw_left, test_hw_right = scipy.signal.peak_widths(test_x * -1, first_event_index, rel_height=0.1)
test_hw_height = test_hw_height * -1

plt.figure()
plt.plot()
plt.plot(test_x)
plt.plot(first_event_index, test_x[first_event_index], 'x')
plt.hlines(test_hw_height, test_hw_left, test_hw_right, color='C2')

test_hw, test_hw_height, test_hw_left, test_hw_right = scipy.signal.peak_widths(test_x * -1, first_event_index, rel_height=0.9)
test_hw_height = test_hw_height * -1

plt.figure()
plt.plot()
plt.plot(test_x)
plt.plot(first_event_index, test_x[first_event_index], 'x')
plt.hlines(test_hw_height, test_hw_left, test_hw_right, color='C2')


