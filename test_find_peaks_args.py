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
lowpass_freq = 500  # Hz
stim_time = 500     # ms
post_stim = 250     # ms, amount of time after stimulus to look for max value
tp_start = 5        # ms, time of start of test pulse
vm_jump = 10        # mV, test pulse voltage jump
pre_tp = 3          # ms, amount of time before test pulse start to get baseline
unit_scaler = -12   # unitless, scaler to get back to A, from pA
amp_factor = 1      # scaler for making plots in pA
fs = 25             # kHz, the sampling frequency

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

non_sub_x = data.traces_filtered.iloc[1000:, 21].values
thresh_non_sub = non_sub_x.std()
x = non_sub_x * -1
# parameters for finding peaks and widths
titles =  ['peak_index', 'prominence', '10 to 90% RT (ms)', 'half-width (ms)']
prominence = thresh_non_sub*3

wlens = [None, 100, 250, 500, 1000, 1250, 2000, 5000, 10000, 25000]

peaks_df = pd.DataFrame()

for wlen in wlens:
    peaks, properties = scipy.signal.find_peaks(x, prominence=prominence, wlen=wlen)

    prominence_data = tuple(properties.values())

    '''plt.figure()
    plt.plot(non_sub_x)
    plt.plot(peaks, non_sub_x[peaks], 'x')
    plt.hlines(base_values, properties['left_bases'], peaks)
    plt.vlines(peaks, non_sub_x[peaks], base_values, color='r')
    '''
    ten_widths, ten_height, ten_left, ten_right = scipy.signal.peak_widths(x, peaks, rel_height=0.9, prominence_data=prominence_data)
    ninety_widths, ninety_height, ninety_left, ninety_right = scipy.signal.peak_widths(non_sub_x * -1, peaks, rel_height=0.1, prominence_data=prominence_data)
    half_widths, hw_height, hw_left, hw_right = scipy.signal.peak_widths(x, peaks, rel_height=0.5, prominence_data=prominence_data)
    full_widths, fw_height, fw_left, fw_right = scipy.signal.peak_widths(x, peaks, rel_height=1, prominence_data=prominence_data)

    hw_height = hw_height * -1
    fw_height = fw_height * -1
    ten_height = ten_height * -1
    ninety_height = ninety_height * -1

    ten_to_ninety = (ninety_left - ten_left) / 25
    prominences = properties['prominences']
    hw_time = half_widths / 25

    titles_iter = [[str(wlen)], titles]
    column_index = pd.MultiIndex.from_product(titles_iter, names=['wlen', 'properties'])
    peaks_array = np.array((peaks, prominences, ten_to_ninety, hw_time)).T
    peaks_data = pd.DataFrame(peaks_array, columns=column_index)
    peaks_df = pd.concat([peaks_df, peaks_data], axis=1)

    plt.figure()
    plt.plot(non_sub_x)
    plt.plot(peaks, non_sub_x[peaks], 'x')
    plt.hlines(hw_height, xmin=hw_left, xmax=hw_right, color='r')
    plt.hlines(fw_height, xmin=fw_left, xmax=fw_right, color='k')
    plt.hlines(ten_height, xmin=ten_left, xmax=ten_right, color='green')
    plt.hlines(ninety_height, xmin=ninety_left, xmax=ninety_right, color='orange')
    plt.title(f'peaks (n={len(peaks)}) and widths with prominence: {prominence}; and wlen: {wlen}')

peak_idx = peaks_df.loc[:, (slice(None), 'peak_index')]
prominence = peaks_df.loc[:, (slice(None), 'prominence')]
rt = peaks_df.loc[:, (slice(None), '10 to 90% RT (ms)')]
hw = peaks_df.loc[:, (slice(None), 'half-width (ms)')]


 '''
 Notes on loop for changing wlen
 - idea of wlen is to limit the signal area from a peak to determine the prominence
    - the problem with ephys data is that the wlen centers on the peak, and goes in both directions, whereas
      ephys data is never (probably?) symmetrical and thus, the winwdows have to be much longer than expected
      based on e.g. rise time.
    - overall i think it probably isn't worth giving it a cut off value
    - need to test with actually evoked events which will likely be much longer
 - explore data with the above variables: peak_idx, prominence, rt, hw
 - making the wlen < 5000 tends to cut off prominence drastically, due to the longer absolute decay time
 - only save a few pA on prominence or < 1 ms on half widths when having 5000 wlen vs no wlen
 - wlen < 2000 miss some peaks
 - wlen < 5000 cuts off a significatn amount of peak, probably too much (based on vis)
 - 10-90 could have some longer (not real) times when there are two close events where the 10% of event 2 is
   clearly before the event starts 
    - solution here would be to make a projection to baseline, which seems complicated, must be a better way
'''
