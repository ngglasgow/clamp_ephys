import clamp_ephys
import pandas as pd
import os

timepoint = 'p2'

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
project_path = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs', timepoint)
tables = project_path.tables
figures = project_path.figures

if timepoint == 'p2':
    notes_path = os.path.join(tables, 'p2_data_notes.csv')
else:
    notes_path = os.path.join(tables, 'p14_data_notes.csv')

''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
lowpass_freq = 500      # Hz
stim_time = 500         # ms
post_stim = 250         # ms, amount of time after stimulus to look for max value
baseline_start = 3000   # ms, time in the sweep to start taking the baseline
baseline_end = 6000     # ms, time in the sweep at which baseline ends
unit_scaler = -12       # unitless, scaler to get back to A, from pA
width = 3               # ms, the required width for classifying an event as a peak


# p14_summary will hold the single line summary for each cell
summary = pd.DataFrame()
kinetics_summary = pd.DataFrame()
meansweep_kinetics_summary = pd.DataFrame()

ncells = len(project_path.paths)

for file_number, cell in enumerate(project_path.paths, 1):
    data = clamp_ephys.workflows.cell(cell, path_to_data_notes=notes_path, timepoint=timepoint)

    data.get_raw_peaks(stim_time, post_stim)
    data.filter_traces(lowpass_freq)
    data.get_filtered_peaks(stim_time, post_stim)
    data.get_series_resistance(unit_scaler)
    data.get_sweep_data()
    data.get_responses(threshold=5)
    data.get_sweepavg_summary()

    summary = pd.concat([summary, data.sweepavg_summary], ignore_index=True)
    
    fig = data.plot_peaks_rs(save_fig=True, path_to_figures=figures)

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

    kinetics_summary = pd.concat([kinetics_summary, data.avg_first3_kinetics_avg_df])

    # kinetics analysis for mean sweep for light conditions only
    if 'light' in cell:
        data.get_peaks_widths(stim_time, width, mean=True)

        if len(data.all_widths_df) == 0:    # skips kinetics analysis for mean traces without peaks
            print('There are no peaks to analyze in the mean trace of {}'.format(data.file_id))
        else:
            data.get_peaks_kinetics(stim_time, mean=True)

            data.save_all_peaks_kinetics(tables, mean=True)
            data.save_all_peaks_kinetics_avg(tables, mean=True)
            data.save_first3_kinetics_avg(tables, mean=True)

            meansweep_kinetics_summary = pd.concat([meansweep_kinetics_summary, data.avg_first3_kinetics_avg_df])
    else:
        pass    

    print(f'Finished analysis for {file_number} of {ncells} cells')

summary_path = os.path.join(tables, f'{timepoint}_summary.csv')
kinetics_summary_path = os.path.join(tables, timepoint, f'{timepoint}_first3_kinetics_summary.csv')
meansweep_kinetics_summary_path = os.path.join(tables, timepoint, f'{timepoint}_meansweep_first3_kinetics_summary.csv')

summary.to_csv(summary_path, float_format='%8.4f', index=False)
kinetics_summary.to_csv(kinetics_summary_path, float_format='%8.4f', index=False)
meansweep_kinetics_summary.to_csv(meansweep_kinetics_summary_path, float_format='%8.4f', index=False)