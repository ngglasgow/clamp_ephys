# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

"""
# from __main__ import *
#
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import elephant
from neo.io import IgorIO
import scipy
import os
from collections import OrderedDict
% matplotlib


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

    # take information from data_notes


    # save metadate into orderedDict pandas DataFrame
    dict = OrderedDict()
    dict['Date'] = date
    dict['Cell ID'] = cell_id
    dict['Cell Number'] = cell_num
    dict['Condition'] = condition
    metadata = pd.DataFrame(dict, index=range(1))


    return metadata


def igor_to_pandas(file):
    '''This function opens an igor binary file (.ibw), extracts the time
    series data, and returns a pandas DataFrame'''

    data_raw = IgorIO(filename=data_dir + file)
    data_neo = data_raw.read_block()
    data_neo_array = data_neo.segments[0].analogsignals[0]
    data_df = pd.DataFrame(data_neo_array.as_array())

    return data_df


def remove_depolarized_sweeps(data, depolarized_sweeps):
    '''
    This function replaces depolarized sweeps with np.nan value to remove
    depolarized sweeps from the record for further analysis.
        data: a pandas DataFrame of current data with each sweep being a column

        depolarized_sweeps: a list of depolarized sweep indices to be removed

    Returns a new cleaned pandas DataFrame with np.nan in columns of sweeps to
    be removed, while maintaining the initial unaltered data
    '''

    cleaned_data = data.copy()
    for sweep in depolarized_sweeps:
        cleaned_data[sweep] = np.nan

    return cleaned_data


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


def epsc_peak(data, baseline, stim_time, polarity='-', post_stim=250, sf=10):
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
        Default is 250 ms.
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
        epsc_peak = peak_window.max()
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
        Amplitude of the test pulse voltage command in mV..
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

# define directories
home_dir = os.path.expanduser('~')
project_dir = home_dir + '/ngglasgow@gmail.com/Data_Urban/jane/OSN Stim Analysis/'
data_dir = '/Volumes/Urban/Huang/'
figure_dir = project_dir + 'analysis/figures/'
table_dir = project_dir + 'analysis/tables/'

''' directory and file parsing '''
# set logic for choosing control and drug files
data_notes = pd.read_csv(table_dir + 'OSN_Gg8vOMP.csv')

# pull out cell_id for directory, file name, and make the full path
file_name_list = data_notes['Cell name'].tolist()
cell_id_list = []

for file in file_name_list:
    metadata = get_metadata(file)
    cell_id = metadata['Cell ID'][0]
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

# %% User Inputs
control_file = 'JH190131_c1_ctrl_1.ibw'
drug_file = 'JH190131_c1_drug_1.ibw'
depolarized_sweeps = []

'''
    Open up control and drug data from each directory
'''

control = igor_to_pandas(control_file)
control_metadata = get_metadata(control_file)
cell_id = control_metadata['Cell ID'][0]

drug = igor_to_pandas(drug_file)
drug_metadata = get_metadata(drug_file)

'''
    Pull out EPSC peak from unfiltered signals
    Baseline 100 ms preceding blue light
    Peak within 250 ms of blue light
'''

control_baseline = mean_baseline(control, 500)
control_peaks = epsc_peak(control, control_baseline, 500)

drug_baseline = mean_baseline(drug, 500)
drug_peaks = epsc_peak(drug, drug_baseline, 500)

all_peaks = pd.concat([control_peaks, drug_peaks], ignore_index=True)

'''
    Pull out EPSC peaks from filtered signals
    Baseline 100 ms preceding blue light
    Peak within 250 ms of blue light
'''

# filter signal with butterworth filter at 1 kHz for control
filt_control = elephant.signal_processing.butter(control,
                                                 lowpass_freq=1000.0,
                                                 fs=10000.0)
filt_control = pd.DataFrame(filt_control)

filt_control_baseline = mean_baseline(filt_control, 500)
filt_control_peaks = epsc_peak(filt_control, filt_control_baseline, 500)

# filter signal with butterworth filter at 1 kHz for drug
filt_drug = elephant.signal_processing.butter(drug,
                                              lowpass_freq=1000.0,
                                              fs=10000.0)
filt_drug = pd.DataFrame(filt_drug)

filt_drug_baseline = mean_baseline(filt_drug, 500)
filt_drug_peaks = epsc_peak(filt_drug, filt_drug_baseline, 500)

filt_all_peaks = pd.concat([filt_control_peaks, filt_drug_peaks], ignore_index=True)

''' Calculating Series Resistance (rs) from test pulse (tp) '''
rs_control = series_resistance(control, 50, -5)
rs_drug = series_resistance(drug, 50, -5)

all_rs = pd.concat([rs_control, rs_drug], ignore_index=True)


''' Plot EPSC peaks and Rs over time of experiemnt '''
# set up index markers for control | drug line and drug stimuli
# pull out number of sweeps for both conditions and all
n_control_sweeps = len(control_peaks)
n_drug_sweeps = len(drug_peaks)
n_all_sweeps = n_control_sweeps + n_drug_sweeps

# set up index for vline separating control and drug
condition_line = n_control_sweeps - 0.5

# set new index for drug_peaks, filtered drug peaks, and rs
drug_peaks_plot = drug_peaks.copy()
drug_peaks_plot.index = np.arange(n_control_sweeps, n_all_sweeps, 1)

filt_drug_peaks_plot = filt_drug_peaks.copy()
filt_drug_peaks_plot.index = np.arange(n_control_sweeps, n_all_sweeps, 1)

rs_drug_plot = rs_drug.copy()
rs_drug_plot.index = np.arange(n_control_sweeps, n_all_sweeps, 1)

# set up auto y max
y_min = all_peaks.min()
y_min_lim = y_min * 1.15 * 1000
% matplotlib inline
# make a figure with 2 plots
fig, axs = plt.subplots(3, 2, figsize=(6, 8), constrained_layout=True)
fig.suptitle('Summary for ' + cell_id)

# optional for plotting unfiltered on same graph for comparison
axs[0, 0].plot(all_peaks*1000, marker='.', color='darkgray', linestyle='', label='raw')

# plot the filterd peak currents NOTE: convert peak values to pA
axs[0, 0].plot(filt_control_peaks*1000, color='k', marker='.', linestyle='', label='Control')
axs[0, 0].plot(filt_drug_peaks_plot*1000, color='r', marker='.', linestyle='', label='TTX + 4-AP')
axs[0, 0].set_xlabel('Stimulus Number')
axs[0, 0].set_ylabel('EPSC Peak (pA)')
axs[0, 0].set_ylim(0, y_min_lim)
axs[0, 0].axvline(x=condition_line, color='k', linestyle='--')

# plot the series resistance values
axs[0, 1].plot(rs_control, marker='.', color='k', linestyle='')
axs[0, 1].plot(rs_drug_plot, marker='.', color='r', linestyle='')
axs[0, 1].set_xlabel('Stimulus Number')
axs[0, 1].set_ylabel('Rs (MOhm)')
axs[0, 1].set_ylim(0, 20)
axs[0, 1].axvline(x=condition_line, color='k', linestyle='--')
axs[0, 0].legend()

''' Plot averaged EPSC trace overlaying all the individual traces '''
# calculate the mean and the SEM of the entire time series
filt_control_subtracted = filt_control - filt_control_baseline
filt_control_mean = filt_control_subtracted.mean(axis=1)
filt_control_sem = filt_control_subtracted.sem(axis=1)
filt_control_std = filt_control_subtracted.std(axis=1)

filt_drug_subtracted = filt_drug - filt_drug_baseline
filt_drug_mean = filt_drug_subtracted.mean(axis=1)
filt_drug_sem = filt_drug_subtracted.sem(axis=1)
filt_drug_std = filt_drug_subtracted.std(axis=1)

# calculate auto y min limit for mean + std
mean_std = (filt_control_mean - filt_control_std)
y_min_mean_std = mean_std[5000:].min()
y_min_mean_lim = y_min_mean_std * 1.1 * 1000

# set up time value for length of traces and window of what to plot
sweep_length = len(control)                  # allow for different sweep length
sweep_time = np.arange(0, sweep_length/10, 0.1)     # time of sweeps in ms

# set up length of line for light stimulation
blue_light = np.arange(500, 550, 0.1)

# plot mean control trace with all traces in gray behind
axs[1, 0].plot(sweep_time, filt_control_subtracted*1000, color='darkgray', linewidth=0.5)
axs[1, 0].plot(sweep_time, filt_control_mean*1000, color='k')
axs[1, 0].hlines(75, 500, 550, color='deepskyblue')
axs[1, 0].set_xlabel('Time (ms)')
axs[1, 0].set_ylabel('Current (pA)')
axs[1, 0].set_xlim(450, 800)
axs[1, 0].set_ylim(y_min_lim, 100)
axs[1, 0].set_title('Control')

# plot mean drug trace with all traces gray behind
axs[1, 1].plot(sweep_time, filt_drug_subtracted*1000, color='darkgray', linewidth=0.5)
axs[1, 1].plot(sweep_time, filt_drug_mean*1000, color='red')
axs[1, 1].hlines(75, 500, 550, color='deepskyblue')
axs[1, 1].set_xlabel('Time (ms)')
axs[1, 1].set_ylabel('Current (pA)')
axs[1, 1].set_xlim(450, 800)
axs[1, 1].set_ylim(y_min_lim, 100)
axs[1, 1].set_title('TTX + 4-AP')

# plot mean control trace with shaded SEM gray behind
axs[2, 0].plot(sweep_time, filt_control_mean*1000, color='k', label='mean')
axs[2, 0].fill_between(sweep_time,
                       (filt_control_mean - filt_control_std) * 1000,
                       (filt_control_mean + filt_control_std) * 1000,
                       color='darkgray',
                       label='st. dev.')
axs[2, 0].hlines(75, 500, 550, color='deepskyblue')
axs[2, 0].set_xlabel('Time (ms)')
axs[2, 0].set_ylabel('Current (pA)')
axs[2, 0].set_xlim(450, 800)
axs[2, 0].set_ylim(y_min_mean_lim, 100)
axs[2, 0].set_title('Control')
axs[2, 0].legend(loc=4)

# plot mean drug trace with shaded SEM gray behind
axs[2, 1].plot(sweep_time, filt_drug_mean*1000, color='red', label='mean')
axs[2, 1].fill_between(sweep_time,
                       (filt_drug_mean - filt_drug_std) * 1000,
                       (filt_drug_mean + filt_drug_std) * 1000,
                       color='darkgray',
                       label='st. dev.')
axs[2, 1].hlines(75, 500, 550, color='deepskyblue')
axs[2, 1].set_xlabel('Time (ms)')
axs[2, 1].set_ylabel('Current (pA)')
axs[2, 1].set_xlim(450, 800)
axs[2, 1].set_ylim(y_min_mean_lim, 100)
axs[2, 1].set_title('TTX + 4-AP')
axs[2, 1].legend(loc=4)
fig
# save figure to file
fig.savefig('{0}{1}_summary.png'.format(figure_dir, cell_id), dpi=300, format='png')

''' Save data to csv files '''
''' Make loopable? '''
''' How to parse directory structure logic? '''
''' Make into autogenerate a jupyter notebook? '''
