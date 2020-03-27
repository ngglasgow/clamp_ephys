import clamp_ephys
import pandas as pd
import os

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
paths = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs')
tables = paths.tables
figures = paths.figures

'''####################### SET PROPER PATH_TO_DATA_NOTES ########################## '''
data_path = r"C:\Users\jhuang\Documents\phd_projects\Injected_GC_data\VC_pairs\tables\\p14_data_notes.csv"


''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
lowpass_freq = 500  # Hz
stim_time = 500     # ms
post_stim = 250     # ms, amount of time after stimulus to look for max value
tp_start = 50        # ms, time of start of test pulse
vm_jump = -5        # mV, test pulse voltage jump
pre_tp = 11          # ms, amount of time before test pulse start to get baseline
unit_scaler = -12   # unitless, scaler to get back to A, from pA
amp_factor = 1000      # scaler for making plots in pA
fs = 10             # kHz, the sampling frequency

'''#################### THIS LOOP RUNS THE SINGLE CELL ANALYSIS #################### '''
# p14_summary will hold the single line summary for each cell
p14_summary = pd.DataFrame()

for path in paths.p14_paths:
    data = clamp_ephys.workflows.cell(path, fs=fs, path_to_data_notes=data_path, timepoint='p14', amp_factor=amp_factor)

    data.get_raw_peaks(stim_time, post_stim)
    data.filter_traces(lowpass_freq)
    data.get_filtered_peaks(stim_time, post_stim)
    data.get_series_resistance(tp_start, vm_jump, pre_tp, unit_scaler)
    data.get_sweep_data()
    data.get_responses(threshold=5)
    data.get_sweepavg_summary()

    p14_summary = pd.concat([p14_summary, data.sweepavg_summary], ignore_index=True)
    
    fig = data.plot_peaks_rs(amp_factor, save_fig=True, path_to_figures=figures)

    data.save_metadata(tables)
    data.save_sweepavg_summary(tables)
    data.save_mean_filtered_trace(tables)
    data.save_mean_subtracted_trace(tables)

summary_path = os.path.join(tables, 'p14_summary.csv')
p14_summary.to_csv(summary_path, float_format='%8.4f', index=False)