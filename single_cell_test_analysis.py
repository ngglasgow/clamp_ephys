import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
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
data = clamp_ephys.workflows.cell(cell_path, fs=fs, path_to_data_notes=data_path, timepoint='p2', amp_factor=amp_factor, drop_sweeps=True)

data.get_raw_peaks(stim_time, post_stim)
data.filter_traces(lowpass_freq)

# collect peak amplitudes and time to peak of all the max peaks
peaks_max = data.get_filtered_peaks(stim_time, post_stim)
time_to_peak_max = (data.peaks_filtered_indices / fs) - stim_time

def plot_half_width(trace, data):
    '''plots a trace with the peak identified and the half width drawn
    trace: int
        index of the trace you want to see
    data: data object
        needs to be a data object
    '''
    x = data.traces_filtered[trace]
    peak = data.peaks_filtered_indices[trace]
    hw_data = data.get_fwhm_peak_max().iloc[trace, 1:].values

    fig, axs = plt.subplots()
    axs.plot(x)
    axs.plot(peak, x[peak], 'x')
    axs.hlines(hw_data[0] * -1, hw_data[1], hw_data[2], color='r')

    return fig
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# test new class function; collect FWHMs of all the peak IPSCs
halfwidths_peak_max = data.get_fwhm_peak_max()['Max peak half-width (ms)']

# pick an example trace to plot the actual half width of a given peak 
%matplotlib widget
fig = plot_half_width(1, data)

import numpy as np # we will use this later, so import it now

from bokeh.io import output_notebook, show
from bokeh.plotting import figure
output_notebook()


trace = data.traces.iloc[:, 0].values
sd = trace[3400 * data.fs:3900 * data.fs].std()
time = np.arange(0, len(trace) / data.fs, 1 / data.fs)
peaks, properties = scipy.signal.find_peaks(trace * -1, prominence=30)
time_peaks = (peaks / data.fs)

wavelet = scipy.signal.ricker(40, 4)
ctime = np.arange(0, len(cwt)/ data.fs, 1 / data.fs)
cwt = scipy.signal.convolve(trace, wavelet)

p = figure(plot_width=800, plot_height=400)

p.line(time, trace, line_width=2)
p.circle(time_peaks, trace[peaks], size=5, line_color="navy", fill_color="orange", fill_alpha=0.5)

p.line(ctime, cwt, color='red')
show(p)

drop_path = '/home/nate/urban/clamp_ephys/test_data/dropped_sweeps.csv'

dropped_sweeps = pd.read_csv(drop_path, index_col=[0])

if data.filename in dropped_sweeps.index:

strsweeps = dropped_sweeps.loc[data.filename].values[0][1:-1].split(', ')
drop_sweeps = [int(sweep) for sweep in strsweeps]
drop_sweeps


show(p)

scipy.signal.convolve()