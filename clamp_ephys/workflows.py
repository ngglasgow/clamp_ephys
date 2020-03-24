from . import clamp
from . import metadata
from . import responses
import elephant
import pandas as pd

class cell:
    def __init__(self, path_to_file, fs, path_to_data_notes):
        self.filepath = path_to_file
        self.filename = self.filepath.split('/')[-1]
        self.fs = fs
        self.notes_path = path_to_data_notes

        self.traces = clamp.igor_to_pandas(self.filepath)
        self.metadata = metadata.get_metadata(self.filename, self.notes_path)


    def filter_traces(self, lowpass_freq):
        '''
        add filtered traces attrbute to data object
        lowpass_freq: frequency in Hz to pass to elephant filter
        '''
        filtered_traces = elephant.signal_processing.butter(self.traces.T, lowpass_freq=lowpass_freq, fs=self.fs * 1000)
        self.filtered_traces = pd.DataFrame(filtered_traces).T


    def get_raw_peaks(self, stim_time, post_stim, polarity='-', pre_stim=100):
        '''
        Finds the baseline and peaks of the raw traces based on passthrough arguments to clamp.
        adds peaks_raw attribute to data object: pandas.Series
        '''
        baseline = clamp.mean_baseline(self.traces, self.fs, stim_time, pre_stim)
        self.peaks_raw = clamp.epsc_peak(self.traces, baseline, self.fs, stim_time, post_stim, polarity)


    def get_filtered_peaks(self, stim_time, post_stim, polarity='-', pre_stim=100):
        '''
        Finds the baseline and peaks of the filtered traces thorugh passthrough arguments to clamp.
        adds peaks_filtered attribute to data object: pandas.Series
        '''
        baseline = clamp.mean_baseline(self.filtered_traces, self.fs, stim_time, pre_stim)
        self.peaks_filtered = clamp.epsc_peak(self.filtered_traces, baseline, self.fs, stim_time, post_stim, polarity)

    
    def get_series_resistance(self, tp_start, vm_jump, pre_tp, unit_scaler):
        '''
        Finds the series resistance of raw traces with passthrough arguments to clamp.
        adds rs attribute to data object: pandas.Series of float in MOhms
        '''
        self.rs = clamp.series_resistance(self.traces, self.fs, tp_start, vm_jump, pre_tp, unit_scaler)


    def __repr__(self):
        return 'Data object for a single cell {}'.format(self.filename)
        

def indiv_cell_analysis(timepoint, file, data_dir, sf=25, amp_factor=1, peak_factor=-12,
                         tp_start=5, vm_jump=10, pre_tp=3):
    '''
    Make and save plots and tables for an individual file

    Parameters
    ----------
    timepoint: str
        Name of the injection timepoint used in the analysis
    file: str
        The file being analyzed
    data_dir: str
        Path to the data folder holding ibw files to be analyzed for that time point
    sf: int or float
        The sampling frequency in kHz. Default is 10 kHz.
    amp_factor: int
        Factor by which to multiply current values by to get pA. Default is 1.
    peak_factor: int or float
        Factor to multiply rs_peak by to get pA. Default is -12.
    tp_start: int or float
        Time in ms when test pulse begins.
    vm_jump: int or float
        Amplitude of whatever windows needs here the test pulse voltage command in mV.
        This is 10 mV in MIES and -5 in Nathan's Igor program.
    pre_tp: int or float
        Time in ms before start of test pulse by which to measure the baseline.

    '''
       
    # gather metadata and set some key parameters for use later on in loop
    metadata = get_metadata(file, data_notes)
    file_path = metadata['Cell Path'][0]
    cell_id = metadata['Cell ID'][0]
    cell_type = metadata['Cell Type'][0]
    condition = metadata['Condition'][0]

   
    ''' Plot EPSC peaks and Rs over time of experiemnt '''
    # set up index markers for data | drug line and drug stimuli
    # pull out number of sweeps for both conditions and all
    n_control_sweeps = len(peaks)

    # set up auto y max for peak plots (min since negative)
    y_min = peaks.min()
    y_min_lim = y_min * 1.15 * amp_factor

    # set up logic for Rs y scaling: if < 20 MOhms, don't scale, if > scale
    if rs.max() <= 20:
        rs_y_min = 0
        rs_y_max = 20
    else:
        rs_y_min = rs.min() * 0.5
        rs_y_max = rs.max() * 1.2

    # make a figure with 2 plots
    fig, axs = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    fig.suptitle('Summary for ' + timepoint + ' ' + cell_type + ' ' + cell_id + ' ' + condition)

    # optional for plotting unfiltered on same graph for comparison
    axs[0, 0].plot(peaks*amp_factor, marker='.', color='darkgray', linestyle='', label='raw')

    # plot the filterd peak currents NOTE: convert peak values to pA
    axs[0, 0].plot(filt_peaks*amp_factor, color='k', marker='.', linestyle='', label='filtered')
    axs[0, 0].set_xlabel('Stimulus Number')
    axs[0, 0].set_ylabel('EPSC Peak (pA)')
    axs[0, 0].set_ylim(0, y_min_lim)
    axs[0, 0].legend()

    # plot the series resistance values
    axs[0, 1].plot(rs, marker='.', color='k', linestyle='')
    axs[0, 1].set_xlabel('Stimulus Number')
    axs[0, 1].set_ylabel('Rs (MOhm)')
    axs[0, 1].set_ylim(rs_y_min, rs_y_max)

    ''' Plot averaged EPSC trace overlaying all the individual traces '''
    # calculate the mean and the SEM of the entire time series
    filt_subtracted = filt_data - filt_baseline
    filt_data_mean = filt_subtracted.mean(axis=1)
    filt_data_sem = filt_subtracted.sem(axis=1)
    filt_data_std = filt_subtracted.std(axis=1)

    # calculate auto y min limit for mean + std
    mean_std = (filt_data_mean - filt_data_std)
    y_min_mean_std = mean_std[5000:].min()
    y_min_mean_lim = y_min_mean_std * 1.1 * amp_factor

    # set up time value for length of traces and window of what to plot
    sweep_length = len(data)                  # allow for different sweep length
    sweep_time = np.arange(0, sweep_length/sf, 1/sf)     # time of sweeps in ms

    # set up length of line for light stimulation
    blue_light = np.arange(500, 550, 1/sf)

    # plot mean data trace with all traces in gray behind
    axs[1, 0].plot(sweep_time, filt_subtracted*amp_factor, color='darkgray', linewidth=0.5)
    axs[1, 0].plot(sweep_time, filt_data_mean*amp_factor, color='k')
    axs[1, 0].hlines(75, 500, 550, color='deepskyblue')
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
    axs[1, 1].hlines(75, 500, 550, color='deepskyblue')
    axs[1, 1].set_xlabel('Time (ms)')
    axs[1, 1].set_ylabel('Current (pA)')
    axs[1, 1].set_xlim(450, 1000)
    axs[1, 1].set_ylim(y_min_mean_lim, 100)
    axs[1, 1].legend(loc=1)

    # fig

    # save figure to file

    fig_save_path = os.path.join(paths.figures, timepoint, cell_type, condition, cell_id)
    fig.savefig(fig_save_path +  '_' + timepoint + '_' + cell_type + '_' + condition + '_summary.png', 
        dpi=300, format='png')
    
    plt.close()

    ''' Save all sweeps metadata for raw, filtered and rs to a csv file '''
    # save each sweep raw peak, filtered peak, and Rs with metadata to summary file
    data_dict = OrderedDict()
    data_dict['Raw Peaks (nA)'] = peaks
    data_dict['Filtered Peaks (nA)'] = filt_peaks
    data_dict['Rs (MOhms)'] = rs
    sweep_data = pd.DataFrame(data_dict)

    # fill in n=sweeps of metadata_df so that can join with peaks for clean df
    metadata_df = pd.DataFrame()
    for i in range(len(peaks)):
        metadata_df = pd.concat([metadata_df, metadata], ignore_index=True)

    # join summary data with metadata
    sweep_meta_data = metadata_df.join(sweep_data)

    # if we decide later on to combine all sweeps and mean, stdev, sem into one
    # file then the following code would be used to make an initial all sweeps
    # df that we can concat the summary measures at the bottom, would need to
    # change it a bit too. I don't think it's necessary now, OK to have 2 files
    # per cell and just open each individually
    # sweep_number = pd.DataFrame(range(len(sweep_data)), columns=['Sweep #'])
    # sweep_number

    # save summary_data to file
    sweep_meta_path = os.path.join(paths.tables, timepoint, cell_type, condition, cell_id)
    sweep_meta_data.to_csv(sweep_meta_path + '_' + timepoint + '_' + cell_type + '_' + condition + 
        '_all_sweeps_data.csv', float_format='%8.4f', index=False, header=True)


    ''' Find the mean peak epscs for raw, filtered, and rs, and save '''
    # make metadata into a series for easy appending
    metaseries = metadata.loc[0]

    # find mean, st dev, and sem of all sweeps for raw, filt, and rs
    mean = metaseries.append(sweep_data.mean())
    std = metaseries.append(sweep_data.std())
    sem = metaseries.append(sweep_data.sem())

    # combine into dataframe, add measure type string and # of sweeps
    summary_data = pd.DataFrame([mean, std, sem])
    measures = pd.DataFrame([['mean'],['st. dev.'],['sem']], columns=['Measure'])
    n_sweeps = pd.DataFrame(len(sweep_data), index=range(3), columns=['# Sweeps'])

    summary_data = pd.concat([summary_data, measures, n_sweeps], axis=1)

    # define path for saving file and save it
    summary_path = os.path.join(paths.tables, timepoint, cell_type, condition, cell_id)
    summary_data.to_csv(summary_path + '_' + timepoint + '_' + cell_type + '_' + condition + 
        '_summary_data.csv', float_format='%8.4f', index=False)


    '''###### save filtered subtracted mean time series data to file #######'''
    mean_timeseries_path = os.path.join(paths.tables, timepoint, cell_type, condition, cell_id)
    filt_data_mean.to_csv(mean_timeseries_path + '_' + timepoint + '_' + cell_type + '_' + condition + 
        '_mean_timeseries.csv', float_format='%8.4f', index=False, header=False)

