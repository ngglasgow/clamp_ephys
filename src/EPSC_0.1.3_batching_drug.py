# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

"""
# from __main__ import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import elephant
from neo.io import IgorIO
import os
from collections import OrderedDict
import math


def get_metadata(file, data_notes):
    '''Takes a filename and parses it for metadata, and returns metadata in an
    orderedDict as a pandas DataFrame for saving later
    Also takes information from the cell spreadsheet in data_notes'''

    # pull out cell id, cell number, date and condition
    file_split = file.split('_')
    cell_id = file_split[0]+'_'+file_split[1]
    cell_num = cell_id[-1:]
    date = '20'+cell_id[2:4]+'-'+cell_id[4:6]+'-'+cell_id[6:8]

    if 'drug' in file:
        condition = 'TTX+4-AP'
    else:
        condition = 'control'

    # grab metadata from data notes spreadsheet
    file_data = data_notes[data_notes['Cell name'] == file]
    cell_path = file_data['File Path'].tolist()[0]
    genotype = file_data['Genotype'].tolist()[0]
    cell_type = file_data['Cell type'].tolist()[0]
    depol_sweep_start = file_data['Depol sweeps start'].tolist()[0]
    depol_sweep_stop = file_data['Depol sweeps stop'].tolist()[0]

    # save metadate into orderedDict pandas DataFrame
    dict = OrderedDict()
    dict['Date'] = date
    dict['Cell ID'] = cell_id
    dict['Cell Number'] = cell_num
    dict['Cell Path'] = cell_path
    # dict['Condition'] = condition
    dict['Genotype'] = genotype
    dict['Cell Type'] = cell_type
    dict['Exclude Sweep Start'] = depol_sweep_start
    dict['Exclude Sweep Stop'] = depol_sweep_stop
    metadata = pd.DataFrame(dict, index=range(1))

    return metadata


def igor_to_pandas(file, data_dir):
    '''This function opens an igor binary file (.ibw), extracts the time
    series data, and returns a pandas DataFrame'''

    file_path = os.path.join(data_dir, file)
    data_raw = IgorIO(filename=file_path)
    data_neo = data_raw.read_block()
    data_neo_array = data_neo.segments[0].analogsignals[0]
    data_df = pd.DataFrame(data_neo_array.as_array())

    return data_df


def mean_baseline(data, stim_time, pre_stim=100, sf=10):
    '''
    Find the mean baseline in a given time series
    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        The time series data for which you want a baseline.
    stim_time: int or float
        The time in ms when stimulus is triggered.
    pre_stim: int or float
        Time in ms before the stimulus trigger over which baseline is measured.
    sf: int or float
        The sampling frequency in kHz.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = (stim_time - pre_stim) * sf
    stop = (stim_time - 1) * sf
    window = data.iloc[start:stop]
    baseline = window.mean()

    return baseline


def epsc_peak(data, baseline, stim_time, polarity='-', post_stim=100, sf=10):
    '''
    Find the peak EPSC value for a pandas.Series or for each sweep (column) of
    a pandas.DataFrame. This finds the absolute peak value of mean baseline
    subtracted data.

    Parameters:
    -----------
    data: pandas.Series or pandas.DataFrame
        Time series data with stimulated synaptic response triggered at the
        same time for each sweep.
    baseline: scalar or pandas.Series
        Mean baseline values used to subtract for each sweep.
    stim_time: int or float
        Time in ms at which stimulus is triggered each sweep.
    polarity: str
        The expected polarity of the EPSC; negative: '-'; postitive: '+'.
        Default is '-'.
    post_stim: int or float
        Time in ms that marks the end of the sampling window post stimulus.
        Default is 100 ms.
    sf: int or float
        The sampling frequency in kHz. Default is 10 kHz.

    Returns
    -------
    epsc_peaks: pandas.Series
        The absolute peak of mean baseline subtracted time series data.
    '''

    subtracted_data = data - baseline
    start = stim_time * sf
    end = (stim_time + post_stim) * sf
    peak_window = subtracted_data.iloc[start:end]

    if polarity == '-':
        epsc_peaks = peak_window.min()
    elif polarity == '+':
        epsc_peaks =peak_window.max()
    else:
        raise ValueError(
            "polarity must either be + or -"
        )

    return epsc_peaks


def series_resistance(data, tp_start, vm_jump, sf=10):
    '''
    Calculate the approximate series resistance (Rs) from a test pulse (tp).

    Parameters
    ----------
    data: pandas.Series or pandas.DataFrame
        Raw time series daata of the v-clamp recording in nA.
    tp_start: int or float
        Time in ms when test pulse begins.
    vm_jump: int or float
        Amplitude ofwhatever windows needs here the test pulse voltage command in mV..
    sf: int or float
        Sampling frequency in kHz. Default is 10 kHz.

    Returns:
    rs: pandas.Series of float
        The series resistance for each sweep in MOhms.
    '''
    # find the baseline 10 ms pre test pulse and subtract from raw data
    rs_baseline = mean_baseline(data, stim_time=tp_start, pre_stim=11)
    rs_subtracted = data - rs_baseline

    # set up indices for starting and ending peak window
    start = tp_start * sf
    end = (tp_start + 2) * sf
    rs_window = rs_subtracted.iloc[start:end]

    if vm_jump > 0:
        rs_peak = rs_window.max()
    else:
        rs_peak = rs_window.min()

    # calculate Rs via V=IR -> Rs = V/I
    rs = ((vm_jump * 10**-3) / (rs_peak * 10**-9)) * 10**-6

    return rs


''' *********************************************************************** '''

''' ################## Define file structure on server #################### '''
# home_dir will depend on the OS, but the rest will not
# query machine identity and set home_dir from there
machine = os.uname()[0]
if machine == 'Darwin':
    home_dir = '/Volumes/Urban'

elif machine == 'Linux':
    home_dir = '/run/user/1000/gvfs/smb-share:server=130.49.237.41,share=urban'

else:
    home_dir = os.path.join('Z:', os.sep)

project_dir = os.path.join(home_dir, 'Huang', 'OSN_OMPvGg8_MTC')
figure_dir = os.path.join(project_dir, 'figures')
table_dir = os.path.join(project_dir, 'tables')
data_dir = os.path.join(project_dir, 'data')

''' ## Open the notes spreadsheet and parse for what we want to analyze ## '''
# open metadata file
data_notes = pd.read_csv(os.path.join(table_dir, 'OSN_Gg8vOMP.csv'))

# pull out cell_id for directory, file name, and make the full path
file_name_list = data_notes['Cell name'].tolist()
cell_id_list = []

for file in file_name_list:
    file_split = file.split('_')
    cell_id = file_split[0]+'_'+file_split[1]
    cell_id_list.append(cell_id)

file_path_list = []

for cell, file in zip(cell_id_list, file_name_list):
    file_path = os.path.join(cell, file + '.ibw')
    file_path_list.append(file_path)

data_notes = pd.concat([pd.DataFrame({'File Path': file_path_list}), data_notes], axis=1)

# drop cells that didn't save to igor
noigor_list = np.array(data_notes[data_notes['Igor saved?'] == 'No'].index)
data_notes = data_notes.drop(index=noigor_list)

# drop cells that don't have any # of drug sweeps
nodrug_list = np.array(data_notes[data_notes['# of drug sweeps'].isnull() == True].index)
data_notes = data_notes.drop(index=nodrug_list)

# update file name list to have only files you want to analyze after logic
file_name_list = data_notes['Cell name'].tolist()

''' ##########################################################################
This is all the analysis, figures, saving
Read in file metadata, open file from igor, convert to pandas
##############################################################################
'''
# loop through all the files in file_name_list for plots and saving
for file_name in file_name_list:
    # set file name from list
    file = file_name

    # gather metadata and set some key parameters for use later on in loop
    metadata = get_metadata(file, data_notes)
    file_path = metadata['Cell Path'][0]
    cell_id = metadata['Cell ID'][0]
    genotype = metadata['Genotype'][0]
    exclude_start = metadata['Exclude Sweep Start'][0]
    exclude_stop = metadata['Exclude Sweep Stop'][0]

    # open igor file and convert to pandas
    data = igor_to_pandas(file_path, data_dir)

    # process logic and build exclude sweeps list from metadata, and exclude sweeps
    if math.isnan(exclude_start) is False:
        # need to pull out the end of the excluded sweeps
        # if all sweeps after start are excluded
        if math.isnan(exclude_stop) is True:
            data = data.iloc[:, :int(exclude_start)]

        # else only exclude traces in between start and stop
        else:
            begin = data.iloc[:, :int(exclude_start)]
            end = data.iloc[:, int(exclude_stop):]
            data = pd.concat([begin, end], axis=1)

    else:
        pass

    '''
        Pull out EPSC peak from unfiltered signals
        Baseline 100 ms preceding blue light
        Peak within 250 ms of blue light
    '''

    baseline = mean_baseline(data, 500)
    peaks = epsc_peak(data, baseline, 500)

    '''
        Pull out EPSC peaks from filtered signals
        Baseline 100 ms preceding blue light
        Peak within 250 ms of blue light
    '''

    # filter signal with butterworth filter at 500 Hz for data
    filt_data = elephant.signal_processing.butter(data.T,
                                                  lowpass_freq=500.0,
                                                  fs=10000.0)
    filt_data = pd.DataFrame(filt_data).T

    filt_baseline = mean_baseline(filt_data, 500)
    filt_peaks = epsc_peak(filt_data, filt_baseline, 500)

    ''' Calculating Series Resistance (rs) from test pulse (tp) '''
    rs = series_resistance(data, 50, -5)

    ''' Plot EPSC peaks and Rs over time of experiemnt '''
    # set up index markers for data | drug line and drug stimuli
    # pull out number of sweeps for both conditions and all
    n_control_sweeps = len(peaks)

    # set up auto y max for peak plots (min since negative)
    y_min = peaks.min()
    y_min_lim = y_min * 1.15 * 1000

    # set up logic for Rs y scaling: if < 20 MOhms, don't scale, if > scale
    if rs.max() <= 20:
        rs_y_min = 0
        rs_y_max = 20
    else:
        rs_y_min = rs.min() * 0.5
        rs_y_max = rs.max() * 1.2

    # make a figure with 2 plots
    fig, axs = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    fig.suptitle('Summary for ' + genotype + ' ' + cell_id)

    # optional for plotting unfiltered on same graph for comparison
    axs[0, 0].plot(peaks*1000, marker='.', color='darkgray', linestyle='', label='raw')

    # plot the filterd peak currents NOTE: convert peak values to pA
    axs[0, 0].plot(filt_peaks*1000, color='k', marker='.', linestyle='', label='filtered')
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
    y_min_mean_lim = y_min_mean_std * 1.1 * 1000

    # set up time value for length of traces and window of what to plot
    sweep_length = len(data)                  # allow for different sweep length
    sweep_time = np.arange(0, sweep_length/10, 0.1)     # time of sweeps in ms

    # set up length of line for light stimulation
    blue_light = np.arange(500, 550, 0.1)

    # plot mean data trace with all traces in gray behind
    axs[1, 0].plot(sweep_time, filt_subtracted*1000, color='darkgray', linewidth=0.5)
    axs[1, 0].plot(sweep_time, filt_data_mean*1000, color='k')
    axs[1, 0].hlines(75, 500, 550, color='deepskyblue')
    axs[1, 0].set_xlabel('Time (ms)')
    axs[1, 0].set_ylabel('Current (pA)')
    axs[1, 0].set_xlim(450, 800)
    axs[1, 0].set_ylim(y_min_lim, 100)

    # plot mean data trace with shaded SEM gray behind
    axs[1, 1].plot(sweep_time, filt_data_mean*1000, color='k', label='mean')
    axs[1, 1].fill_between(sweep_time,
                           (filt_data_mean - filt_data_std) * 1000,
                           (filt_data_mean + filt_data_std) * 1000,
                           color='darkgray',
                           label='st. dev.')
    axs[1, 1].hlines(75, 500, 550, color='deepskyblue')
    axs[1, 1].set_xlabel('Time (ms)')
    axs[1, 1].set_ylabel('Current (pA)')
    axs[1, 1].set_xlim(450, 800)
    axs[1, 1].set_ylim(y_min_mean_lim, 100)
    axs[1, 1].legend(loc=1)

    # fig

    # save figure to file
    fig_save_path = os.path.join(figure_dir, genotype + '_' + cell_id)
    fig.savefig(fig_save_path + '_summary.png', dpi=300, format='png')

    ''' Save all sweeps data for raw, filtered and rs to a csv file '''
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
    sweep_meta_path = os.path.join(table_dir, genotype + '_' + cell_id)
    sweep_meta_data.to_csv(sweep_meta_path + '_all_sweeps_data.csv', float_format='%8.4f', index=False, header=True)


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
    summary_path = os.path.join(table_dir, genotype + '_' + cell_id)
    summary_data.to_csv(summary_path + '_summary_data.csv', float_format='%8.4f', index=False)


'''take the means of raw, filt, and rs and save to file with single row'''
''' Make into autogenerate a jupyter notebook? '''
