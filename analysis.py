import clamp_ephys
import pandas as pd
import os

timepoint = 'p14'

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
project_path = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs', timepoint)
tables = project_path.tables
figures = project_path.figures

if timepoint == 'p2':
    notes_path = os.path.join(tables, 'p2_data_notes.csv')
else:
    notes_path = os.path.join(tables, 'p14_data_notes.csv')


# p14_summary will hold the single line summary for each cell
summary = pd.DataFrame()
kinetics_summary = pd.DataFrame()
meansweep_kinetics_summary = pd.DataFrame()

ncells = len(project_path.paths)

for file_number, cell in enumerate(project_path.paths, 1):
    data = clamp_ephys.workflows.cell(cell, path_to_data_notes=notes_path, timepoint=timepoint)

    data.get_raw_peaks()
    data.filter_traces()
    data.get_filtered_peaks()
    data.get_series_resistance()
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
    data.get_peaks_widths()
    data.get_peaks_kinetics()

    data.save_all_peaks_kinetics(tables)
    data.save_all_peaks_kinetics_avg(tables)
    data.save_first3_kinetics_avg(tables)

    kinetics_summary = pd.concat([kinetics_summary, data.avg_first3_kinetics_avg_df])

    # kinetics analysis for mean sweep
    data.get_peaks_widths(mean=True)
    data.get_peaks_kinetics(mean=True)

    data.save_all_peaks_kinetics(tables, mean=True)
    data.save_all_peaks_kinetics_avg(tables, mean=True)
    data.save_first3_kinetics_avg(tables, mean=True)

    meansweep_kinetics_summary = pd.concat([meansweep_kinetics_summary, data.avg_first3_kinetics_avg_df])

    print(f'Finished analysis for {file_number} of {ncells} cells')

summary_path = os.path.join(tables, '{timepoint}_summary.csv')
kinetics_summary_path = os.path.join(tables, timepoint, '{timepoint}_first3_kinetics_summary.csv')
meansweep_kinetics_summary_path = os.path.join(tables, timepoint, '{timepoint}_meansweep_first3_kinetics_summary.csv')

summary.to_csv(summary_path, float_format='%8.4f', index=False)
kinetics_summary.to_csv(kinetics_summary_path, float_format='%8.4f', index=False)
meansweep_kinetics_summary.to_csv(meansweep_kinetics_summary_path, float_format='%8.4f', index=False)