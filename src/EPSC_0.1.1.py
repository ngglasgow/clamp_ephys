# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

"""
#%%
#from __main__ import *
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import elephant
from neo.io import IgorIO
import scipy
import os
% matplotlib

# %%

# define directories
home_dir = os.path.expanduser('~')
project_dir = home_dir + '/ngglasgow@gmail.com/Data_Urban/jane/OSN Stim Analysis/'
data_dir = project_dir + 'ETC/Data/'
figure_dir = project_dir + 'analysis/figures/'
table_dir = project_dir + 'analysis/tables/'

# %% User Inputs
file_name = 'JH190131_c1_drug_1.ibw'
depolarized_sweeps = []

# metadata
file_name_split = file_name.split('_')
cell_id = file_name_split[0]+'_'+file_name_split[1]
cell_num = cell_id[-1:]
cell_tag = file_name.split('.')[0]

# pull out date code and cell numbers
date = '20'+cell_id[2:4]+'-'+cell_id[4:6]+'-'+cell_id[6:8]

# define condition e.g. control, drug1, drug2
if 'drug' in file_name:
    condition = 'TTX+4-AP'
else:
    condition = 'control'


# %%
'''
This cell sorts through the cell_id directory and sorts files only noise vm
files into two lits; 1 for drugs and 1 for contorl. Then it opens the correct
file type for the given file_name, and concatenates vm traces if necessary.
Open igor binary file into neo, read neo.block, convert to np.array, convert
to pd.dataframe, concat with empty data frame noise_vmdf and repeat if concat.
'''

currdf = pd.DataFrame()
temp_curr_raw = IgorIO(filename=data_dir + file_name)
temp_curr_neo = temp_curr_raw.read_block()
temp_currarray = temp_curr_neo.segments[0].analogsignals[0]
temp_currdf = pd.DataFrame(temp_currarray.as_array())
currdf = pd.concat([currdf, temp_currdf], axis=1)

plt.figure()
plt.plot(currdf)

plt.figure()
plt.plot(currdf.iloc[:, 0])

for item in depolarized_sweeps:
    noise_vmdf[item] = np.nan

# %% Exclude sweeps
exclude = depolarized_sweeps


''' Pull out EPSC peak from unfiltered signals '''

# %% Extract window for baseline: 400-500 ms (100 ms preceding blue light)
epsc_baseline_window = currdf.iloc[4000:5000]
epsc_baselines = epsc_baseline_window.mean()

# subtract mean baseline from each trace
epsc_subtracted = currdf - epsc_baselines

# Extract window for epsc peak: 500-750 ms (250 ms after light onset)
epsc_window = epsc_subtracted.iloc[5000:7500]
epsc_peak = epsc_window.min()


''' Pull out EPSC peaks from filtered signals '''

# %% filter signal with butterworth filter at 1 kHz
filt_currdf = elephant.signal_processing.butter(currdf, lowpass_freq=1000.0, fs=10000.0)
filt_currdf = pd.DataFrame(filt_currdf)

# Extract window for baseline: 400-500 ms (100 ms preceding blue light)
filt_epsc_baseline_window = filt_currdf.iloc[4000:5000]
filt_epsc_baselines = filt_epsc_baseline_window.mean()

# subtract mean baseline from each trace
filt_epsc_subtracted = filt_currdf - filt_epsc_baselines

# Extract window for epsc peak: 500-750 ms (250 ms after light onset)
filt_epsc_window = filt_epsc_subtracted.iloc[5000:7500]
filt_epsc_peak = filt_epsc_window.min()


''' Optional peak vs. filter checking '''
peak_delta = epsc_peak - filt_epsc_peak
peak_delta.mean()
peak_delta
plt.figure()
plt.plot(epsc_subtracted.iloc[:, 2])
plt.plot(filt_epsc_subtracted.iloc[:, 2])


''' Calculating Series Resistance (rs) from test pulse (tp) '''
# %% Calculating series resistance.
rs_baseline = currdf.iloc[390:490].mean()   # baseline 10 ms before test pulse
rs_window = currdf.iloc[500:520]            # window to find min peak
rs_peak_current = rs_window.min()           # abs minimum of cap. transient

# set change in I and change in V
delta_i = (rs_peak_current - rs_baseline) * 10**-9  # change in I after baseline subtraction in A
delta_v = -0.005                            # voltage step in V

# calculate Rs with V=IR -> R = V/I
rs = delta_v/delta_i                        # returns Rs in Ohms
rs_Mohms = rs * 10**-6


''' Plot EPSC peaks and Rs over time of experiemnt '''

# set up a dynamic time value that scales with experiment length
n_stimuli = len(filt_epsc_peak)             # number of stimuli in experiment
stim_freq = 10                              # the interstimulus interval in s
duration = n_stimuli*stim_freq              # duration of experiment
time = np.arange(0, duration, stim_freq)    # create time axis for duration

# make a figure with 2 plots
fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
fig.suptitle('Summary for {0} in {1}'.format(cell_id, condition))

# plot the filterd peak currents NOTE: convert peak values to pA
axs[0].plot(time, filt_epsc_peak*1000, marker='.', linestyle='', label='filtered')
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('EPSC Peak (pA)')
axs[0].set_ylim(0, -500)

# plot the series resistance values
axs[1].plot(time, rs_Mohms, marker='.', linestyle='')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Rs (MOhm)')
axs[1].set_ylim(0, 30)

# save figure to file
fig.savefig('{0}epsc_rs{1}.png'.format(figure_dir, cell_tag), dpi=300, format='png')

# optional for plotting unfiltered on same graph for comparison
axs[0].plot(time, epsc_peak*1000, marker='.', linestyle='', label='raw')
axs[0].legend()

# if want filtered vs. raw saved always then just overwrite other file
fig.savefig('{0}epsc_rs{1}.png'.format(figure_dir, cell_tag), dpi=300, format='png')


''' Plot averaged EPSC trace overlaying all the individual traces '''

''' Save data to csv files '''
