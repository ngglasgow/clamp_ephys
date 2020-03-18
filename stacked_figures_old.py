import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from neo.io import IgorIO
import os
from os import listdir
from os.path import isfile, join


def igor_to_pandas(file, data_dir):
    '''This function opens an igor binary file (.ibw), extracts the time
    series data, and returns a pandas DataFrame'''

    file_path = os.path.join(data_dir, file)
    data_raw = IgorIO(filename=file_path)
    data_neo = data_raw.read_block()
    data_neo_array = data_neo.segments[0].analogsignals[0]
    data_df = pd.DataFrame(data_neo_array.as_array())

    return data_df


def mean_baseline(data, stim_time, pre_stim=10, fs=10):
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
    fs: int or float
        The sampling frequency in kHz.

    Returns
    -------
    baseline: float or pandas.Series
        The mean baseline over the defined window
    '''
    start = (stim_time - pre_stim) * fs
    stop = (stim_time - 1) * fs
    window = data.iloc[start:stop]
    baseline = window.mean()

    return baseline


# need to set up the appropriate paths here, and may not need all of these
home_dir = os.path.join('Z:', os.sep)

project_dir = os.path.join(home_dir, 'Huang', 'OSN_OMPvGg8_MTC')
figure_dir = r'C:\Users\Jane Huang\Documents\Grad school\Urban Lab\Data Analysis\Injected_GC_data\RMS_4wpi_4AP\data\CsCl\light100\figure'
table_dir = os.path.join(project_dir, 'tables')
#data_dir = os.path.join(project_dir, 'data')

''' I had to add a r wtf'''
data_dir = r'C:\Users\Jane Huang\Documents\Grad school\Urban Lab\Data Analysis\Injected_GC_data\RMS_4wpi_4AP\data\CsCl\light100'


# your file parsing loop that pulls out waves_mean and makes a waves_mean_df
# goes here

file_name_list = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]
waves_mean_df = pd.DataFrame()


for file in file_name_list:
    cell_waves_df = igor_to_pandas(file, data_dir)
    cell_waves_mean = cell_waves_df.mean(axis=1)
    waves_mean_df = pd.concat([waves_mean_df, cell_waves_mean], axis=1)

test_df = igor_to_pandas(file_name_list[10], data_dir)

for i in range(len(test_df.columns)):
    plt.figure()
    plt.plot(test_df.iloc[:,i])
    
''' I had to rename the columns for indexing'''    
waves_mean_df.columns = range(len(waves_mean_df.columns))

#mean_of_means = all_means.mean(axis=1)

''' ########### stacking figures with no frame or scales ################## '''
# first two lines are for if you want to make plots have unequal dimensions
# gs_kw = dict(height_ratios=[0.5, 0.7, 0.7, 0.7, 4, 1])
# fig, axs = plt.subplots(6, 8, figsize=(20, 5.5), gridspec_kw=gs_kw, constrained_layout=True)

''' some parameters for automating figure sizing and rows '''
n_cells = len(waves_mean_df.columns)
width = 2                   # plot width in inches
height = n_cells * 0.5      # height in inches of figure, scaler is inch per plot
stim_time = 500             # time of stimulus onset in ms
stim_length = 100           # length of stimulus in ms
fs = 10                     # sampling frequency in kHz
light_amplitude = 0.045     # plot amp in nA where you want blue bar
plotx_start = 4500          # index of time axis you want to start plotting
plotx_end = 10000           # index of time axis you want to stop plotting
yscaler = 1.1               # scaler for min/max if exceeding
ploty_max = 0.050           # default scaling values for y
ploty_min = -0.100          # default scaling values for y

# if you want whatever auto ticks are, find x or y tick_list and comment out
x_tick_list = [5000, 10000, 15000]  # set x ticks for making scale bar
y_tick_list = [-0.1, -0.05, 0]        # set y ticks for making scale bar

# create list for storing max and min values for scaling
wave_max_list = []
wave_min_list = []

# generate the figure with all the subplots
fig, axs = plt.subplots(n_cells, 1, figsize=(width, height), constrained_layout=True)

# pull out one wave, subtract baseline, take slice, and plot
for i in range(n_cells):
    # baseline subtract for each wave
    wave = waves_mean_df[i]
    baseline = mean_baseline(wave, stim_time)
    baseline_subtracted_wave = wave - baseline

    # take slice of baseline subtracted wave and calculate max/min for scaling
    wave_slice = baseline_subtracted_wave[plotx_start:plotx_end]

    wave_max = wave_slice.max()
    wave_max_list.append(wave_max)

    wave_min = wave_slice.min()
    wave_min_list.append(wave_min)

    # actually plot the baseline subtracted wave
    axs[i].plot(wave_slice, color='darkgray')

    # remove frame y_ticks and x_ticks for all to make only traces
    axs[i].set_frame_on(False)
    axs[i].set_xticks([])
    axs[i].set_yticks([])

# plot a blue line
axs[0].hlines(light_amplitude, stim_time * fs, (stim_time + stim_length) * fs, color='deepskyblue', linewidth=2)

# calculations and logic for max/min scaling as same for all axes
actual_max = max(wave_max_list)         # find max of all waves
actual_min = min(wave_min_list)         # find min of all waves

if actual_max < ploty_max:
    if actual_min > ploty_min:
        for ax in axs.flat:
            ax.set_ylim(ploty_min, ploty_max)
    else:
        for ax in axs.flat:
            ax.set_ylim(actual_min*yscaler, ploty_max)
else:
    if actual_min > ploty_min:
        for ax in axs.flat:
            ax.set_ylim(ploty_min, actual_max*yscaler)
    else:
        for ax in axs.flat:
            ax.set_ylim(actual_min*yscaler, actual_max*yscaler)
# fig
# save figure to file
fig_save_path = os.path.join(figure_dir, 'CsCl_light100_1000ms_end')
fig.savefig(fig_save_path + '_waves_final.png', dpi=300, format='png')


''' ################### figure with x scale for making axis scale bar ######'''
# generate the figure with all the subplots
fig_xscale, axs = plt.subplots(n_cells, 1, figsize=(width, height), constrained_layout=True)

for i in range(n_cells):

    # baseline subtract for each wave
    wave = waves_mean_df[i]
    baseline = mean_baseline(wave, stim_time)
    baseline_subtracted_wave = wave - baseline

    # take slice of baseline subtracted wave and calculate max/min for scaling
    wave_slice = baseline_subtracted_wave[plotx_start:plotx_end]

    wave_max = wave_slice.max()
    wave_max_list.append(wave_max)

    wave_min = wave_slice.min()
    wave_min_list.append(wave_min)

    # actually plot the baseline subtracted wave
    axs[i].plot(wave_slice, color='darkgray')

    # remove frame y_ticks and x_ticks for all to make only traces
    axs[i].set_frame_on(False)
    axs[i].set_xticks(x_tick_list)
    axs[i].set_yticks([])

# calculations and logic for max/min scaling as same for all axes
actual_max = max(wave_max_list)         # find max of all waves
actual_min = min(wave_min_list)         # find min of all waves

if actual_max < ploty_max:
    if actual_min > ploty_min:
        for ax in axs.flat:
            ax.set_ylim(ploty_min, ploty_max)
    else:
        for ax in axs.flat:
            ax.set_ylim(actual_min*yscaler, ploty_max)
else:
    if actual_min > ploty_min:
        for ax in axs.flat:
            ax.set_ylim(ploty_min, actual_max*yscaler)
    else:
        for ax in axs.flat:
            ax.set_ylim(actual_min*yscaler, actual_max*yscaler)

# fig_xscale

# save figure to file
fig_save_path = os.path.join(figure_dir, 'CsCl_light100_1000ms_end')
fig_xscale.savefig(fig_save_path + '_waves_xscale.png', dpi=300, format='png')


''' ############ figure with y scale for making axis scale bar ############ '''
# generate the figure with all the subplots
fig_yscale, axs = plt.subplots(n_cells, 1, figsize=(width, height), constrained_layout=True)

for i in range(n_cells):

    # baseline subtract for each wave
    wave = waves_mean_df[i]
    baseline = mean_baseline(wave, stim_time)
    baseline_subtracted_wave = wave - baseline

    # take slice of baseline subtracted wave and calculate max/min for scaling
    wave_slice = baseline_subtracted_wave[plotx_start:plotx_end]

    wave_max = wave_slice.max()
    wave_max_list.append(wave_max)

    wave_min = wave_slice.min()
    wave_min_list.append(wave_min)

    # actually plot the baseline subtracted wave
    axs[i].plot(wave_slice, color='darkgray')

    # remove frame y_ticks and x_ticks for all to make only traces
    axs[i].set_frame_on(False)
    axs[i].set_xticks([])
    axs[i].set_yticks(y_tick_list)

# calculations and logic for max/min scaling as same for all axes
actual_max = max(wave_max_list)         # find max of all waves
actual_min = min(wave_min_list)         # find min of all waves

if actual_max < ploty_max:
    if actual_min > ploty_min:
        for ax in axs.flat:
            ax.set_ylim(ploty_min, ploty_max)
    else:
        for ax in axs.flat:
            ax.set_ylim(actual_min*yscaler, ploty_max)
else:
    if actual_min > ploty_min:
        for ax in axs.flat:
            ax.set_ylim(ploty_min, actual_max*yscaler)
    else:
        for ax in axs.flat:
            ax.set_ylim(actual_min*yscaler, actual_max*yscaler)

# fig_yscale

# save figure to file
fig_save_path = os.path.join(figure_dir, 'CsCl_light100_1000ms_end')
fig_yscale.savefig(fig_save_path + '_waves_yscale.png', dpi=300, format='png')
