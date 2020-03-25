from . import clamp
from . import metadata
from . import responses
import elephant
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import os
import platform

class cell:
    def __init__(self, path_to_file, fs, path_to_data_notes, timepoint):
        self.filepath = path_to_file
        self.filename = self.filepath.split('/')[-1]
        self.fs = fs
        self.notes_path = path_to_data_notes

        self.traces = clamp.igor_to_pandas(self.filepath)
        self.metadata = metadata.get_metadata(self.filename, self.notes_path)
        self.cell_id = self.metadata['Cell ID'][0]
        self.condition = self.metadata['Condition'][0]
        self.cell_type = self.metadata['Cell Type'][0]
        self.timepoint = timepoint


    def filter_traces(self, lowpass_freq):
        '''
        add filtered traces attrbute to data object
        lowpass_freq: frequency in Hz to pass to elephant filter
        '''
        traces_filtered = elephant.signal_processing.butter(self.traces.T, lowpass_freq=lowpass_freq, fs=self.fs * 1000)
        self.traces_filtered = pd.DataFrame(traces_filtered).T


    def get_raw_peaks(self, stim_time, post_stim, polarity='-', pre_stim=100):
        '''
        Finds the baseline and peaks of the raw traces based on passthrough arguments to clamp.
        adds peaks_raw attribute to data object: pandas.Series
        '''
        self.baseline_raw = clamp.mean_baseline(self.traces, self.fs, stim_time, pre_stim)
        self.peaks_raw = clamp.epsc_peak(self.traces, self.baseline_raw, self.fs, stim_time, post_stim, polarity)


    def get_filtered_peaks(self, stim_time, post_stim, polarity='-', pre_stim=100):
        '''
        Finds the baseline and peaks of the filtered traces thorugh passthrough arguments to clamp.
        adds peaks_filtered attribute to data object: pandas.Series
        '''
        self.baseline_filtered = clamp.mean_baseline(self.traces_filtered, self.fs, stim_time, pre_stim)
        self.peaks_filtered = clamp.epsc_peak(self.traces_filtered, self.baseline_filtered, self.fs, stim_time, post_stim, polarity)

    
    def get_series_resistance(self, tp_start, vm_jump, pre_tp, unit_scaler):
        '''
        Finds the series resistance of raw traces with passthrough arguments to clamp.
        adds rs attribute to data object: pandas.Series of float in MOhms
        '''
        self.rs = clamp.series_resistance(self.traces, self.fs, tp_start, vm_jump, pre_tp, unit_scaler)


    def plot_peaks_rs(self, amp_factor):
        '''
        Takes the data traces and plots the current summary of peaks plot
        Parameters
        ----------
        amp_factor: int
            is for scaling current values to = pA
        timepoint: str
            for labeling graph, what injection timepoint p2 or p14
        
        Returns
        -------
        fig: matplotlib.pyplot fig
            the figure object created
        '''
        # set up auto y max for peak plots (min since negative)
        y_min = self.peaks_filtered.min()
        y_min_lim = y_min * 1.15 * amp_factor

        # set up logic for Rs y scaling: if < 20 MOhms, don't scale, if > scale
        if self.rs.max() <= 20:
            rs_y_min = 0
            rs_y_max = 20
        else:
            rs_y_min = self.rs.min() * 0.5
            rs_y_max = self.rs.max() * 1.2
        
        # make a figure with 2 plots
        fig, axs = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
        fig.suptitle('Summary for {} {} {} {}'.format(self.timepoint, self.cell_type, self.cell_id, self.condition))
        # optional for plotting unfiltered on same graph for comparison
        axs[0, 0].plot(self.peaks_raw*amp_factor, marker='.', color='darkgray', linestyle='', label='raw')

        # plot the filterd peak currents NOTE: convert peak values to pA
        axs[0, 0].plot(self.peaks_filtered*amp_factor, color='k', marker='.', linestyle='', label='filtered')
        axs[0, 0].set_xlabel('Stimulus Number')
        axs[0, 0].set_ylabel('EPSC Peak (pA)')
        axs[0, 0].set_ylim(0, y_min_lim)
        axs[0, 0].legend()

        # plot the series resistance values
        axs[0, 1].plot(self.rs, marker='.', color='k', linestyle='')
        axs[0, 1].set_xlabel('Stimulus Number')
        axs[0, 1].set_ylabel('Rs (MOhm)')
        axs[0, 1].set_ylim(rs_y_min, rs_y_max)

        ''' Plot averaged EPSC trace overlaying all the individual traces '''
        # calculate the mean and the SEM of the entire time series
        filt_subtracted = self.traces_filtered - self.baseline_filtered
        filt_data_mean = filt_subtracted.mean(axis=1)
        filt_data_std = filt_subtracted.std(axis=1)

        # calculate auto y min limit for mean + std
        mean_std = (filt_data_mean - filt_data_std)
        y_min_mean_std = mean_std[5000:].min()
        y_min_mean_lim = y_min_mean_std * 1.1 * amp_factor

        # set up time value for length of traces and window of what to plot
        sweep_length = len(self.traces)                  # allow for different sweep length
        sweep_time = np.arange(0, sweep_length/self.fs, 1/self.fs)     # time of sweeps in ms

        # set up length of line for light stimulation
        blue_start = 500    # ms, time blue light comes on
        blue_stop = 550     # ms, time blue light turns off

        # plot mean data trace with all traces in gray behind
        axs[1, 0].plot(sweep_time, filt_subtracted*amp_factor, color='darkgray', linewidth=0.5)
        axs[1, 0].plot(sweep_time, filt_data_mean*amp_factor, color='k')
        axs[1, 0].hlines(75, blue_start, blue_stop, color='deepskyblue')
        axs[1, 0].set_xlabel('Time (ms)')
        axs[1, 0].set_ylabel('Current (pA)')
        axs[1, 0].set_xlim(450, 1000)
        axs[1, 0].set_ylim(y_min_lim, 100)

        # plot mean data trace with shaded SEM gray behind
        axs[1, 1].plot(sweep_time, filt_data_mean*amp_factor, color='k', label='mean')
        axs[1, 1].fill_between(sweep_time,
                            (filt_data_mean - filt_data_std) * amp_factor,
                            (filt_data_mean + filt_data_std) * amp_factor,
                            color='darkgray',
                            label='st. dev.')
        axs[1, 1].hlines(75, blue_start, blue_stop, color='deepskyblue')
        axs[1, 1].set_xlabel('Time (ms)')
        axs[1, 1].set_ylabel('Current (pA)')
        axs[1, 1].set_xlim(450, 1000)
        axs[1, 1].set_ylim(y_min_mean_lim, 100)
        axs[1, 1].legend(loc=1)

        return fig


    def save_fig(self, path_to_figures, figure):
        '''
        Saves the figure object to the path_to_figures
        Parameters
        ----------
        path_to_figures: str
            path to the figure directory
        figure: plt.pyplot fig
            figure object
        '''
        filename = '{}_{}_{}_{}_summary.png'.format(self.cell_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_figures, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        figure.savefig(path, dpi=300, format='png')


    def save_metadata(self, path_to_tables):
        '''
        Takes the metadata, appends the sweep information and saves it
        Parameters:
        -----------
        path_to_tables: str
            path to the directory for tables
        '''
        data_dict = OrderedDict()
        data_dict['Raw Peaks (nA)'] = self.peaks_raw
        data_dict['Filtered Peaks (nA)'] = self.peaks_filtered
        data_dict['Rs (MOhms)'] = self.rs
        self.sweep_data = pd.DataFrame(data_dict)

        # fill in n=sweeps of metadata_df so that can join with peaks for clean df
        metadata_df = pd.DataFrame()
        for i in range(len(self.peaks_raw)):
            metadata_df = pd.concat([metadata_df, metadata], ignore_index=True)

        # join summary data with metadata
        sweep_meta_data = metadata_df.join(self.sweep_data)

        filename = '{}_{}_{}_{}_all_sweeps_data.csv'.format(self.cell_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        sweep_meta_data.to_csv(path, float_format='%8.4f', index=False, header=True)


    def save_mean_data(self, path_to_tables):
        '''
        Takes the metadata and sweep data finds the means and appends to save summary to tables
        Parameters
        ----------
        path_to_tables: str
            path to the tables directory
        '''
        # make metadata into a series for easy appending
        metaseries = self.metadata.loc[0]

        # find mean, st dev, and sem of all sweeps for raw, filt, and rs
        mean = metaseries.append(self.sweep_data.mean())
        std = metaseries.append(self.sweep_data.std())
        sem = metaseries.append(self.sweep_data.sem())

        # combine into dataframe, add measure type string and # of sweeps
        summary_data = pd.DataFrame([mean, std, sem])
        measures = pd.DataFrame([['mean'],['st. dev.'],['sem']], columns=['Measure'])
        n_sweeps = pd.DataFrame(len(self.sweep_data), index=range(3), columns=['# Sweeps'])

        summary_data = pd.concat([summary_data, measures, n_sweeps], axis=1)

        # define path for saving file and save it
        filename = '{}_{}_{}_{}_summary_data.csv'.format(self.cell_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        summary_data.to_csv(path, float_format='%8.4f', index=False)


    def save_mean_filtered_trace(self, path_to_tables):
        '''
        Saves the mean trace from the self.traces_filtered time series
        Parameters
        ----------
        path_to_tables: str
            path to the tables directory
        '''
        mean_filtered_trace = self.traces_filtered.mean(axis=1)
        filename = '{}_{}_{}_{}_mean_timeseries.csv'.format(self.cell_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        mean_filtered_trace.to_csv(path, float_format='%8.4f', index=False, header=False)


    def __repr__(self):
        return 'Data object for a single cell {}'.format(self.filename)
        

class file_structure:
    def __init__(self, location, project_path):
        '''
        Creates an object with paths as attributes:
        location:   str value 'local' or 'server' only, refers to where you are
                    doing the actual work, 'local' by default.
        project_path:   str of the root project path from wherever your home dir is
        '''
        machine = platform.uname()[0]

        if location == 'local':
            if machine == 'Darwin':
                home_dir = '/Volumes/Urban'

            elif machine == 'Linux':
                home_dir = os.path.join(os.path.expanduser('~'), 'urban/neurobio/Huang')

            elif machine == 'Windows':
                home_dir = r"C:\Users\jhuang\Documents\phd_projects"

            else:
                print("OS not recognized. \nPlease see Nate for correction.")

        elif location == 'server':
            if machine == 'Darwin':
                home_dir = '/Volumes/Urban'

            elif machine == 'Linux':
                home_dir = os.path.join(os.path.expanduser('~'), 'urban/neurobio/Huang')

            elif machine == 'Windows':
                home_dir = r"N:\Huang"

            else:
                print("OS not recognized. \nPlease see Nate for correction.")

        self.project = os.path.join(home_dir, project_path)
        self.figures = os.path.join(self.project, 'figures')
        self.tables = os.path.join(self.project, 'tables')
        self.p2 = os.path.join(self.project, 'data/p2')
        self.p2_files = os.listdir(self.p2)
        self.p14 = os.path.join(self.project, 'data/p14')
        self.p14_files = os.listdir(self.p14)

    def __repr__(self):
        return 'Project file structure and file lists for {}'.format(self.project)
