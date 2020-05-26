import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import numpy as np

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
width = 3               # ms, the width necessary for an event to be identified as synaptic

'''#################### THIS LOOP RUNS THE SINGLE CELL ANALYSIS #################### '''
cell_path = '/home/jhuang/Documents/phd_projects/Injected_GC_data/VC_pairs/data/p2/JH200303_c8_light100.ibw'
data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=data_path, timepoint='p2', amp_factor=amp_factor)

data.get_raw_peaks(stim_time, post_stim, polarity='-', baseline_start=baseline_start, baseline_end=baseline_end)
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
# peaks = data.peaks_filtered_indices
# subtracted_data = data.traces_filtered - data.baseline_filtered

# finding events using subtracted, unfiltered data
sweeps = data.traces - data.new_baseline_raw

# drop first 520 s from sweep to account for TP and only look at time from stimulus onset to end of sweep
window_start = (stim_time + 20) * fs

# this below would be outside of the per-cell loop
p2_kinetics_summary = pd.DataFrame()

columns_index = ['peak_index', 'prominence', '10 to 90% RT (ms)', 'tau', 'half-width (ms)']
first_peaks_properties_df = pd.DataFrame()

for sweep in range(len(sweeps.columns)):
    
    trace = sweeps.iloc[window_start:, sweep].values
    thresh = 2.5 * trace.std()

    # finding and plotting all peaks
    peaks, properties = scipy.signal.find_peaks(trace * -1, distance=0.5*fs, prominence=thresh, width=width*fs)

    if len(peaks) == 0:
        print('No peaks in sweep {}'.format(sweep))
    else:
        prominence_data = list(properties.values())[0:3]
        fig = plt.figure()
        fig.suptitle('Sweep {}'.format(sweep))
        plt.plot(trace)
        plt.plot(peaks, trace[peaks], 'x')

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

        # fit exponential curve to decay phase of each peak to find tau
        # decay phase is defined by time of peak to the right index at full width (fw_right)
        peak_index = 0
        tau_list = []

        # for peak in peaks:
            
        #     def decay_func(time, current_peak, tau):
	    #         return current_peak * np.exp(-time/tau)

        #     decay_end = fw_right[peak_index].astype(int) # indexing needs to be a whole number
        #     time_xdata = np.arange(peak, decay_end+1)
        #     current_ydata = sweeps.iloc[peak:decay_end+1, sweep].values * -1
            
        #     current_peak0 = sweeps.iloc[peak, sweep] * -1
        #     tau0 = 2*fs
        #     starting_params = [current_peak0, tau0]
            
        #     popt, pcov = scipy.optimize.curve_fit(f=decay_func, xdata=time_xdata, ydata=current_ydata, p0=starting_params)
        #     current_peak, tau = popt

        #     tau_list.append(tau)
            
        #     peak_index += 1

        peaks_array = np.array((peaks, prominences, ten_to_ninety, hw_time)).T
        peaks_data = pd.DataFrame(peaks_array, columns=columns_index)

        # extracting kinetics of the first three responses/peaks
        first_peaks_data = peaks_data.iloc[:3] #this slice indexing works even if len(peaks) < 3
        
        # extract delay to response (index to first peak)
        delay_to_response = peaks[0] / fs
            
        avg_first_peaks = pd.DataFrame(first_peaks_data.mean(axis=0)).T
        avg_first_peaks.insert(0, 'delay_to_response (ms)', delay_to_response)    
        first_peaks_properties_df = pd.concat([first_peaks_properties_df, avg_first_peaks], axis=0, ignore_index=True)


first_peaks_properties_avg = pd.DataFrame(first_peaks_properties_df.mean(axis=0)).T
p2_kinetics_summary = pd.concat([p2_kinetics_summary, first_peaks_properties_avg], axis=0, ignore_index=True)

# '''
# # finding latency to response, i.e. the time at which current exceeds 3x std of baseline
# baseline_start = baseline_start * fs
# baseline_end = baseline_end * fs
# baseline = new_mean_baseline(data, fs, baseline_start, baseline_end)
# baseline_std = new_std_baseline(data, fs, baseline_start, baseline_end)
# '''

# # find 10-90% rise time
# # asking to eval at rel_height equivalent to evaluating at height of Peak(height) - prominence * rel_height
# # so evaluating at rel_height = 0.1 does not equal evaluating at 10% of Peak(height)
# test_hw, test_hw_height, test_hw_left, test_hw_right = scipy.signal.peak_widths(test_x * -1, first_event_index, rel_height=0.1)
# test_hw_height = test_hw_height * -1

# plt.figure()
# plt.plot()
# plt.plot(test_x)
# plt.plot(first_event_index, test_x[first_event_index], 'x')
# plt.hlines(test_hw_height, test_hw_left, test_hw_right, color='C2')

# test_hw, test_hw_height, test_hw_left, test_hw_right = scipy.signal.peak_widths(test_x * -1, first_event_index, rel_height=0.9)
# test_hw_height = test_hw_height * -1

# plt.figure()
# plt.plot()
# plt.plot(test_x)
# plt.plot(first_event_index, test_x[first_event_index], 'x')
# plt.hlines(test_hw_height, test_hw_left, test_hw_right, color='C2')


