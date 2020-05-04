import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib.pyplot as plt

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
paths = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs')
tables = paths.tables
figures = paths.figures

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
cell_path = os.path.join(os.getcwd(), 'test_data', 'JH200303_c1_light100.ibw')
data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=data_path, timepoint='p2', amp_factor=amp_factor)

data.get_raw_peaks(stim_time, post_stim)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)

def plot_half_width(trace, data):
    '''plots a trace with the peak identified and the half width drawn
    trace: int
        index of the trace you want to see
    data: data object
        needs to be a data object
    '''
    x = data.traces_filtered[trace]
    peak = data.peaks_filtered_indices[trace]
    hw_data = data.max_peak_half_widths.iloc[trace, 1:].values

    fig, axs = plt.subplots()
    axs.plot(x)
    axs.plot(peak, x[peak], 'x')
    axs.hlines(hw_data[0] * -1, hw_data[1], hw_data[2], color='r')

    return fig
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# test new class function
data.get_max_peak_half_width()

# pick an example trace to plot the actual half width of a given peak 
%matplotlib

fig = plot_half_width(1, data)

