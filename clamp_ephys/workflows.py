from . import clamp
from . import metadata
import elephant
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import os
import platform
import scipy

class cell:
    def __init__(self, path_to_file, fs, path_to_data_notes, timepoint, amp_factor, drop_sweeps=False):
        self.filepath = path_to_file

        machine = platform.uname()[0]
        if machine == 'Windows':
            self.filename = self.filepath.split('\\')[-1]
        else:
            self.filename = self.filepath.split('/')[-1]

        self.fs = fs
        self.notes_path = path_to_data_notes
        self.file_id = self.filename.split('.')[0]
        self.traces = clamp.igor_to_pandas(self.filepath) * amp_factor

        if drop_sweeps is True:
            filename_length = len(self.filename) + 1
            dropped_path = os.path.join(self.filepath[:-filename_length], 'dropped_sweeps.csv')
            dropped_sweeps = pd.read_csv(dropped_path, index_col=[0])
            strsweeps = dropped_sweeps.loc[self.filename].values[0][1:-1].split(', ')
            if '' in strsweeps:
                pass
            else:
                sweeps_to_drop = [int(sweep) for sweep in strsweeps]
                self.traces.drop(columns=sweeps_to_drop, inplace=True)

        self.mean_traces = self.traces.mean(axis=1)
        self.time = np.arange(0, len(self.traces), 1000 / self.fs)
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
        self.mean_traces_filtered = self.traces_filtered.mean(axis=1)

    def get_first_raw_peak(self, stim_time, post_stim, polarity='-', pre_stim=100):
        '''
        Finds the baseline and the first peak, defined as the first event to exceed
        threshold greater than 3 std above the baseline (defined as last 3s of the sweep)
        '''
        



    def get_raw_peaks(self, stim_time, post_stim, polarity='-', pre_stim=100, baseline_start=3000, baseline_end=6000):
        '''
        Finds the baseline and peaks of the raw traces based on passthrough arguments to clamp.
        adds peaks_raw attribute to data object: pandas.Series
        '''
        self.baseline_raw = clamp.mean_baseline(self.traces, self.fs, stim_time, pre_stim)
        self.new_baseline_raw = clamp.new_mean_baseline(self.traces, self.fs, baseline_start, baseline_end)
        self.peaks_raw = clamp.epsc_peak(self.traces, self.baseline_raw, self.fs, stim_time, post_stim, polarity)
        self.mean_baseline_raw = clamp.mean_baseline(self.mean_traces, self.fs, stim_time, pre_stim)
        self.mean_peak_raw, peak_index = clamp.epsc_peak(self.mean_traces, self.mean_baseline_raw, self.fs, stim_time, post_stim, polarity, index=True)
        self.mean_peak_raw_std = self.traces.std(axis=1)[peak_index]
        self.mean_peak_raw_sem = self.traces.sem(axis=1)[peak_index]


    def get_filtered_peaks(self, stim_time, post_stim, polarity='-', pre_stim=100):
        '''
        Finds the baseline and peaks of the filtered traces thorugh passthrough arguments to clamp.
        adds peaks_filtered attribute to data object: pandas.Series
        '''
        self.baseline_filtered = clamp.mean_baseline(self.traces_filtered, self.fs, stim_time, pre_stim)
        self.baseline_filtered_std = clamp.std_baseline(self.traces_filtered, self.fs, stim_time)
        self.peaks_filtered, self.peaks_filtered_indices = clamp.epsc_peak(self.traces_filtered, self.baseline_filtered, self.fs, stim_time, post_stim, polarity, index=True)
        self.mean_baseline_filtered = clamp.mean_baseline(self.mean_traces_filtered, self.fs, stim_time, pre_stim)
        self.mean_baseline_std_filtered = clamp.std_baseline(self.mean_traces_filtered, self.fs, stim_time) 
        self.mean_peak_filtered, self.mean_peak_index = clamp.epsc_peak(self.mean_traces_filtered, self.mean_baseline_filtered, self.fs, stim_time, post_stim, polarity, index=True)
        self.mean_peak_filtered_std = self.traces_filtered.std(axis=1)[self.mean_peak_index]
        self.mean_peak_filtered_sem = self.traces_filtered.sem(axis=1)[self.mean_peak_index]
        self.mean_peak_filtered_time = self.mean_peak_index / self.fs


    def get_fwhm_peak_max(self):
        n = len(self.traces_filtered.columns)
        hw_df = pd.DataFrame()

        if self.mean_peak_filtered > 0:
            invert = 1
        
        else:
            invert = -1

        for i in range(n):
            x = self.traces_filtered[i].values * invert
            peak = [self.peaks_filtered_indices[i]]
            hw, hw_height, hw_left, hw_right = scipy.signal.peak_widths(x, peak, rel_height=0.5)
            hw_time = hw / self.fs
            hw_peak = pd.DataFrame({'Max peak half-width (ms)': hw_time, 'HW height': hw_height, 'HW left index': hw_left, 'HW right index': hw_right}, index=range(1))
            hw_df = pd.concat([hw_df, hw_peak], ignore_index=True)

        self.max_peak_half_widths = hw_df

        return hw_df


    def get_filtered_peaks_kinetics(self, fs, stim_time, post_stim, polarity='-'):
        '''
        Finds the peak in each filtered trace and calculates latency to start, full width at half max,
        total charge transferred, time to peak, duration        
        '''

        self.individual_peaks_filtered, self.individual_peaks_index = self.traces_filtered.apply(clamp.epsc_peak(self.traces_filtered, self.baseline_filtered, self.fs, stim_time, post_stim, polarity, index=True))
        self.individual_subtracted_traces = self.traces_filtered - self.baseline_filtered

        # define window to perform analyses on        
        start = stim_time * fs
        end = (stim_time + post_stim) * fs
        self.traces_filtered_windowed = self.individual_subtracted_traces.iloc[start:end]
        
        # this is latency to start (time to get to 20% of peak IPSC)
        self.current_start = self.traces_filtered_windowed.argmax(self.traces_windowed < (self.individual_peaks_filtered * 0.2))

        # this is finding full width at half max
        self.fwhm = scipy.signal.peak_widths(self.traces_filtered_windowed, self.individual_peaks_index, rel_height=0.5)
        
        # plt.plot(self.traces_filtered_windowed)
        # plt.plot(self.individual_peaks_index, self.traces_filtered_windowed[self.individual_peaks_index], "x")
        # plt.hlines(*self.fwhm[1:], color="C2")
        # plt.show()
      


    def get_series_resistance(self, tp_start, vm_jump, pre_tp):
        '''
        Finds the series resistance of raw traces with passthrough arguments to clamp.
        adds rs attribute to data object: pandas.Series of float in MOhms
        '''
        unit_scaler = -12
        self.rs = clamp.series_resistance(self.traces, self.fs, tp_start, vm_jump, pre_tp, unit_scaler)


    def get_sweep_data(self):
        '''
        Takes all the data and returns a DataFrame with sweep_data.
        '''
        data_dict = OrderedDict()
        data_dict['Raw Peaks (pA)'] = self.peaks_raw
        data_dict['Filtered Peaks (pA)'] = self.peaks_filtered
        data_dict['Rs (MOhms)'] = self.rs
        self.sweep_data = pd.DataFrame(data_dict)


    def get_responses(self, threshold=None):
        '''
        Decides on whether there is a response above 2x, 3x above the baseline std,
        or a user-selectable cutoff.
        Parameters
        ----------
        threshold: int, float (optional)
            If supplied, will provide another threshold in addition to the 2x and 3x
            above the baseline std to threshold the response checker.
        Returns
        -------
        self.responses: pd.DataFrame(bool)
            A DataFrame with bool for responses above the threshold in the column header.
        '''

        baseline_std = self.mean_baseline_std_filtered
        peak_mean = self.mean_peak_filtered.mean()

        response_2x = abs(peak_mean) > baseline_std * 2
        response_3x = abs(peak_mean) > baseline_std * 3

        if threshold is None:
            self.responses = pd.DataFrame({'Response 2x STD': response_2x,
                                           'Response 3x STD': response_3x},
                                           index=range(1))
        else:
            response_threshold = abs(peak_mean) > baseline_std * threshold
            response_string = 'Response {}x STD'.format(threshold)

            self.responses = pd.DataFrame({'Response 2x STD': response_2x,
                                           'Response 3x STD': response_3x,
                                           response_string: response_threshold},
                                           index=range(1))


    def get_sweepavg_summary(self):
        '''
        Accumulates all data and reports metadata, responses, means, stds, and sems.
        '''
        summary_data = pd.DataFrame({'Raw SweepAvg Peak (pA)': self.mean_peak_raw,
                                     'Raw SweepAvg Peak STD (pA)': self.mean_peak_raw_std,
                                     'Raw SweepAvg Peak SEM (pA)': self.mean_peak_raw_sem,
                                     'Filtered SweepAvg Peak (pA)': self.mean_peak_filtered,
                                     'Filtered SweepAvg Peak STD (pA)': self.mean_peak_filtered_std,
                                     'Filtered SweepAvg Peak SEM (pA)': self.mean_peak_filtered_sem,
                                     }, index=range(1))

        self.sweepavg_summary = pd.concat([self.metadata, self.responses, summary_data], axis=1)


    def get_summary_data(self):
        '''
        Accumulates all data and reports metadata, responses, means, stds, and sems.
        '''
        mean = self.sweep_data.mean()
        mean.index = mean.index + ' mean'
        std = self.sweep_data.std()
        std.index = std.index + ' std'
        sem = self.sweep_data.sem()
        sem.index = sem.index + ' sem'

        # put data together in single row of a pd.DataFrame
        summary_data = pd.DataFrame(pd.concat([mean, std, sem])).T

        # add in metadata and responses
        self.summary_data = pd.concat([self.metadata, self.responses, summary_data], axis=1)


    def plot_sweeps(self, sweep, stim_time, baseline_start, save_fig=False, path_to_figures=None):
        '''
        Plots all the sweeps in every cell.
        Parameters
        ----------
        sweep: int
            sweep number being plotted, 0-indexed
        stim_time: int
            time of stimulus onset, ms
        baseline_start: int
            time to start measuring the baseline, ms
        save_fig: bool 
            tells function to either save and close the plot (true) or display the plot (false)
        path_to_figures: str
            path to figures IF save_fig=True
        Returns
        -------
        fig: matplotlib.pyplot fig
            the figure object created
        '''

        window_start = (stim_time + 20) * self.fs
        baseline_window_start = baseline_start * self.fs

        # using filtered, non-subtracted data so I can see which sweeps to drop
        # window omits TP and pre-stimulus time
        x = self.traces_filtered.iloc[window_start:, sweep].values
        baseline = self.traces_filtered.iloc[baseline_window_start:, sweep].values
        thresh = 3 * baseline.std()
        sweep_length = len(x)
        sweep_time = np.arange(0, sweep_length/self.fs, 1/self.fs)

        # finding all peaks
        # peaks, properties = scipy.signal.find_peaks(x * -1, prominence=thresh)
        
        # correct peaks time for fs
        # peaks_corrected = peaks/self.fs

        fig = plt.figure()
        fig.suptitle('Sweep {}'.format(sweep))
        plt.plot(sweep_time, x)
        # plt.plot(peaks_corrected, x[peaks], 'x')

        filename = '{}_sweep_{}.png'.format(self.file_id, sweep)
        base_path = os.path.join(path_to_figures, self.file_id)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        fig.savefig(path, dpi=300, format='png')
        plt.close(fig)



    def plot_peaks_rs(self, amp_factor, save_fig=False, path_to_figures=None):
        '''
        Takes the data traces and plots the current summary of peaks plot
        Parameters
        ----------
        amp_factor: int
            is for scaling current values to = pA
        timepoint: str
            for labeling graph, what injection timepoint p2 or p14
        save_fig: bool 
            tells function to either save and close the plot (true) or display the plot (false)
        path_to_figures: str
            path to figures IF save_fig=True
        Returns
        -------
        fig: matplotlib.pyplot fig
            the figure object created
        '''
        # set up auto y max for peak plots (min since negative)
        y_min = self.peaks_filtered.min()
        y_min_lim = y_min * 1.15

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
        axs[0, 0].plot(self.peaks_raw, marker='.', color='darkgray', linestyle='', label='raw')

        # plot the filterd peak currents NOTE: convert peak values to pA
        axs[0, 0].plot(self.peaks_filtered, color='k', marker='.', linestyle='', label='filtered')
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
        y_min_mean_lim = y_min_mean_std * 1.1

        # set up time value for length of traces and window of what to plot
        sweep_length = len(self.traces)                  # allow for different sweep length
        sweep_time = np.arange(0, sweep_length/self.fs, 1/self.fs)     # time of sweeps in ms

        # set up length of line for light stimulation
        blue_start = 500    # ms, time blue light comes on
        blue_stop = 550     # ms, time blue light turns off

        # plot mean data trace with all traces in gray behind
        axs[1, 0].plot(sweep_time, filt_subtracted, color='darkgray', linewidth=0.5)
        axs[1, 0].plot(sweep_time, filt_data_mean, color='k')
        axs[1, 0].hlines(75, blue_start, blue_stop, color='deepskyblue')
        axs[1, 0].set_xlabel('Time (ms)')
        axs[1, 0].set_ylabel('Current (pA)')
        axs[1, 0].set_xlim(450, 1000)
        axs[1, 0].set_ylim(y_min_lim, 100)

        # plot mean data trace with shaded SEM gray behind
        axs[1, 1].plot(sweep_time, filt_data_mean, color='k', label='mean')
        axs[1, 1].fill_between(sweep_time,
                            (filt_data_mean - filt_data_std),
                            (filt_data_mean + filt_data_std),
                            color='darkgray',
                            label='st. dev.')
        axs[1, 1].hlines(75, blue_start, blue_stop, color='deepskyblue')
        axs[1, 1].set_xlabel('Time (ms)')
        axs[1, 1].set_ylabel('Current (pA)')
        axs[1, 1].set_xlim(450, 1000)
        axs[1, 1].set_ylim(y_min_mean_lim, 100)
        axs[1, 1].legend(loc=1)

        if save_fig is False:
            return fig

        elif save_fig is True:
            
            filename = '{}_{}_{}_{}_summary.png'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
            base_path = os.path.join(path_to_figures, self.timepoint, self.cell_type, self.condition)
            metadata.check_create_dirs(base_path)

            path = os.path.join(base_path, filename)
            fig.savefig(path, dpi=300, format='png')
            plt.close()


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
        filename = '{}_{}_{}_{}_summary.png'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
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
        # join summary data with metadata
        sweep_meta_data = self.metadata.join(self.sweep_data, how='right')
        sweep_meta_data.fillna(method='ffill', inplace=True)

        filename = '{}_{}_{}_{}_all_sweeps_data.csv'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        sweep_meta_data.to_csv(path, float_format='%8.4f', index=False, header=True)


    def save_summary_data(self, path_to_tables):
        '''
        Takes the metadata and sweep data finds the means and appends to save summary to tables
        Parameters
        ----------
        path_to_tables: str
            path to the tables directory
        '''
        # define path for saving file and save it
        filename = '{}_{}_{}_{}_summary_data.csv'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        self.summary_data.to_csv(path, float_format='%8.4f', index=False)


    def save_sweepavg_summary(self, path_to_tables):
        '''
        Takes the metadata and sweep data finds the means and appends to save summary to tables
        Parameters
        ----------
        path_to_tables: str
            path to the tables directory
        '''
        # define path for saving file and save it
        filename = '{}_{}_{}_{}_sweepavg_summary.csv'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        self.sweepavg_summary.to_csv(path, float_format='%8.4f', index=False)

    def save_mean_filtered_trace(self, path_to_tables):
        '''
        Saves the mean trace from the self.traces_filtered time series
        Parameters
        ----------
        path_to_tables: str
            path to the tables directory
        '''
        filename = '{}_{}_{}_{}_mean_timeseries.csv'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        self.mean_traces_filtered.to_csv(path, float_format='%8.4f', index=False, header=False)

    def save_mean_subtracted_trace(self, path_to_tables):
        '''
        Saves the mean trace from the self.mean_traces_filtered time series
        Parameters
        ----------
        path_to_tables: str
            path to the tables directory
        '''
        subtracted_trace = self.mean_traces_filtered - self.mean_baseline_filtered

        filename = '{}_{}_{}_{}_mean_subtracted_timeseries.csv'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        subtracted_trace.to_csv(path, float_format='%8.4f', index=False, header=False)

    def save_mean_peak_time(self, path_to_tables):
        '''
        Saves the mean trace from the self.mean_traces_filtered time series
        Parameters
        ----------
        path_to_tables: str
            path to the tables directory
        '''
        filename = '{}_{}_{}_{}_mean_peak_time.csv'.format(self.file_id, self.timepoint, self.cell_type, self.condition)
        base_path = os.path.join(path_to_tables, self.timepoint, self.cell_type, self.condition)
        metadata.check_create_dirs(base_path)

        path = os.path.join(base_path, filename)
        peak_time = pd.DataFrame(self.mean_peak_filtered_time)
        peak_time.to_csv(path, float_format='%8.4f', index=False, header=False)
 
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
                home_dir = '/home/jhuang/Documents/phd_projects'

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
        self.p2 = os.path.join(self.project, 'data', 'p2')
        self.p2_paths = [os.path.join(self.p2, file) for file in os.listdir(self.p2)]
        self.p14 = os.path.join(self.project, 'data', 'p14')
        self.p14_paths = [os.path.join(self.p14, file) for file in os.listdir(self.p14)]

    def __repr__(self):
        return 'Project file structure and file lists for {}'.format(self.project)
