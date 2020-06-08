import clamp_ephys
import pandas as pd
import os

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
project_path = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs')
tables = project_path.tables
figures = project_path.figures

'''####################### SET PROPER PATH_TO_DATA_NOTES ########################## '''
notes_path = os.path.join(tables, 'p14_data_notes.csv')


''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
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

'''#################### THIS LOOP RUNS THE SINGLE CELL ANALYSIS #################### '''
# p14_summary will hold the single line summary for each cell
p14_summary = pd.DataFrame()
p14_kinetics_summary = pd.DataFrame()
p14_meansweep_kinetics_summary = pd.DataFrame()

ncells = len(project_path.p14_paths)

for file_number, cell in enumerate(project_path.p14_paths, 1):
    data = clamp_ephys.workflows.cell(cell, fs=fs, path_to_data_notes=notes_path, timepoint=timepoint, amp_factor=amp_factor)

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
    # data.save_mean_peak_time(tables)

    # kinetics analysis for all sweeps
    data.get_peaks_widths(stim_time, width)
    data.get_peaks_kinetics(stim_time)

    data.save_all_peaks_kinetics(tables)
    data.save_all_peaks_kinetics_avg(tables)
    data.save_first3_kinetics_avg(tables)

    p14_kinetics_summary = pd.concat([p14_kinetics_summary, data.avg_first3_kinetics_avg_df])

    # kinetics analysis for mean sweep
    data.get_peaks_widths(stim_time, width, mean=True)
    data.get_peaks_kinetics(stim_time, mean=True)

    data.save_all_peaks_kinetics(tables, mean=True)
    data.save_all_peaks_kinetics_avg(tables, mean=True)
    data.save_first3_kinetics_avg(tables, mean=True)

    p14_meansweep_kinetics_summary = pd.concat([p14_meansweep_kinetics_summary, data.avg_first3_kinetics_avg_df])

    print(f'Finished analysis for {file_number} of {ncells} cells')

summary_path = os.path.join(tables, 'p14_summary.csv')
kinetics_summary_path = os.path.join(tables, timepoint, 'p14_first3_kinetics_summary.csv')
meansweep_kinetics_summary_path = os.path.join(tables, timepoint, 'p14_meansweep_first3_kinetics_summary.csv')

p14_summary.to_csv(summary_path, float_format='%8.4f', index=False)
p14_kinetics_summary.to_csv(kinetics_summary_path, float_format='%8.4f', index=False)
p14_meansweep_kinetics_summary.to_csv(meansweep_kinetics_summary_path, float_format='%8.4f', index=False)