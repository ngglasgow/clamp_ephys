# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

"""
# from __main__ import *

import pandas as pd
import numpy as np
import elephant
from neo.io import IgorIO
import os
import platform


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


def std_baseline(data, stim_time, pre_stim=100, sf=10):
    '''
    Find the baseline st. dev. in a given time series
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
    std: float or pandas.Series
        The baseline st. dev. over the defined window
    '''
    start = (stim_time - pre_stim) * sf
    stop = (stim_time - 1) * sf
    window = data.iloc[start:stop]
    std = window.std()

    return std


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
        epsc_peaks = peak_window.max()
    else:
        raise ValueError(
            "polarity must either be + or -"
        )

    return epsc_peaks


''' *********************************************************************** '''

''' ################## Define file structure on server #################### '''
# home_dir will depend on the OS, but the rest will not
# query machine identity and set home_dir from there
machine = platform.uname()[0]

if machine == 'Darwin':
    home_dir = '/Volumes/Urban'

elif machine == 'Linux':
    home_dir = '/run/user/1000/gvfs/smb-share:server=130.49.237.41,share=urban'

elif machine == 'Windows':
    home_dir = os.path.join('Z:', os.sep, 'urban')

else:
    print("OS not recognized. \nPlease see Nate for correction.")

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
data_notes = data_notes.drop(noigor_list)

# update file name list to have only files you want to analyze after logic
file_name_list = data_notes['Cell name'].tolist()
file_path_list = data_notes['File Path'].tolist()

''' ##########################################################################
This is all the analysis, figures, saving
Read in file metadata, open file from igor, convert to pandas
##############################################################################
'''
# make empty lists for response answers for each file
responses_df = pd.DataFrame()

# loop through all the files in file_name_list for plots and saving
for file_name, file_path in zip(file_name_list, file_path_list):
    # set file name from list
    file = file_name

    # open igor file and convert to pandas
    data = igor_to_pandas(file_path, data_dir)

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
    filt_baseline_std = std_baseline(filt_data, 500)
    filt_peaks = epsc_peak(filt_data, filt_baseline, 500, post_stim=50)

    # make binary choice of whether the response is more than 3x std
    mean_std = filt_baseline_std.mean()
    mean_peaks = filt_peaks.mean()

    response_2x = abs(mean_peaks) > mean_std * 2
    response_3x = abs(mean_peaks) > mean_std * 3

    responses = pd.DataFrame({'Cell name': file_name,
                              'Mean Peaks (nA)': mean_peaks,
                              'Mean STD (nA)': mean_std,
                              'Response 2x STD': response_2x,
                              'Response 3x STD': response_3x,
                              },
                              index=range(1))

    responses_df = pd.concat([responses_df, responses], ignore_index=True)

# save responses as csv
responses_df.to_csv(os.path.join(table_dir, 'responses.csv'))


# merge with data_notes and save
responses_data_notes = pd.merge(data_notes, responses_df)
responses_data_notes.to_csv(os.path.join(table_dir, 'responses_data_notes.csv'),
                            float_format='%8.4f', index=False)

''' part to count number of responders for each genotype '''

gg8 = responses_data_notes[responses_data_notes['Genotype'] == 'Gg8']
gg8_yes = gg8[gg8['Response 3x STD'] == True]['Cell name'].tolist()

omp = responses_data_notes[responses_data_notes['Genotype'] == 'OMP']
omp_yes = omp[omp['Response 3x STD'] == True]['Cell name'].tolist()

gg8_yes_ids = []

for name in gg8_yes:
    name_split = name.split('_')
    cell_id = name_split[0] + '_' + name_split[1]
    gg8_yes_ids.append(cell_id)

omp_yes_ids = []

for name in omp_yes:
    name_split = name.split('_')
    cell_id = name_split[0] + '_' + name_split[1]
    omp_yes_ids.append(cell_id)

n_gg8_yes = len(set(gg8_yes_ids))
n_gg8 = len(gg8)
print("Responses in {} out of {} cells in Gg8+ slices".format(n_gg8_yes, n_gg8))

n_omp_yes = len(set(omp_yes_ids))
n_omp = len(omp)
print("Responses in {} out of {} cells in OMP+ slices".format(n_omp_yes, n_omp))

test = responses_data_notes.Genotype.unique()
for genotype in test:
    print(genotype)

responses_data_notes['Cell type'].unique().tolist()
responses_data_notes
