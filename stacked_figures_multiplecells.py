import clamp_ephys
import pandas as pd
import os

%matplotlib
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
figure_dir = r'C:\Users\Jane Huang\Documents\Grad school\Urban Lab\Data Analysis\Injected_GC_data\GC_attch\figures'
table_dir = os.path.join(project_dir, 'tables')
#data_dir = os.path.join(project_dir, 'data')

''' I had to add a r wtf'''
data_dir = r'C:\Users\Jane Huang\Documents\Grad school\Urban Lab\Data Analysis\Injected_GC_data\New_VC_pairs\JH190829\slice_2\MC\data'


# your file parsing loop that pulls out waves_mean and makes a waves_mean_df
# goes here

file_name_list = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]
waves_mean_df = pd.DataFrame()


for file in file_name_list:
    cell_waves_df = igor_to_pandas(file, data_dir)
    cell_waves_mean = cell_waves_df.mean(axis=1)
    waves_mean_df = pd.concat([waves_mean_df, cell_waves_mean], axis=1, ignore_index=True)
    


'''
test_df = igor_to_pandas(file_name_list[10], data_dir)

for i in range(len(test_df.columns)):
    plt.figure()
    plt.plot(test_df.iloc[:, i])

'''
''' ########### stacking figures with no frame or scales ################## '''
# first two lines are for if you want to make plots have unequal dimensions
# gs_kw = dict(height_ratios=[0.5, 0.7, 0.7, 0.7, 4, 1])
# fig, axs = plt.subplots(6, 8, figsize=(20, 5.5), gridspec_kw=gs_kw, constrained_layout=True)


''' some parameters for automating figure sizing and rows '''
n_sweeps = len(cell_waves_df.columns)
n_to_plot = 5
width = 3                   # plot width in inches
height = n_sweeps * 0.5      # height in inches of figure, scaler is inch per plot
stim_time = 500             # time of stimulus onset in ms
stim_length = 100           # length of stimulus in ms
fs = 10                     # sampling frequency in kHz
light_amplitude = 0.045     # plot amp in nA where you want blue bar
plotx_start = 4500          # index of time axis you want to start plotting
plotx_end = 10000           # index of time axis you want to stop plotting
yscaler = 1.1               # scaler for min/max if exceeding
ploty_max = 0.050           # default scaling values for y
ploty_min = -0.100          # default scaling values for y
x_scalebar = 100            # x scalebar length in ms
y_scalebar = 0.05           # y scalebar height in nA

# if you want whatever auto ticks are, find x or y tick_list and comment out
x_tick_list = [5000, 10000, 15000]  # set x ticks for making scale bar
y_tick_list = [-0.1, -0.05, 0]        # set y ticks for making scale bar

# create list for storing max and min values for scaling
wave_max_list = []
wave_min_list = []

# generate the figure with all the subplots
fig, axs = plt.subplots(n_to_plot, 1, figsize=(width, height), constrained_layout=True)

# pull out one wave, subtract baseline, take slice, and plot
for i in range(n_to_plot):
    # baseline subtract for each wave
    wave = cell_waves_df[i]
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

    # choose whether you want blue line or blue box, this is blue box
    axs[i].axvspan(stim_time * fs, (stim_time + stim_length) * fs, facecolor='cornflowerblue', alpha=0.5)

# plot a blue line, comment out if just want box
# axs[0].hlines(light_amplitude, stim_time * fs, (stim_time + stim_length) * fs, color='deepskyblue', linewidth=2)

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
fig_save_path = os.path.join(figure_dir, 'JH190626_c1_attch_light_1')
fig.savefig(fig_save_path + '_waves_final.png', dpi=300, format='png')

''' ############### same figure as above but with scales ################## '''
# generate the figure with all the subplots
fig_scales, axs = plt.subplots(n_to_plot, 1, figsize=(width, height), constrained_layout=True)

# pull out one wave, subtract baseline, take slice, and plot
for i in range(n_to_plot):
    # baseline subtract for each wave
    wave = cell_waves_df[i]
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

    # add blue shading to all graphs
    axs[i].axvspan(stim_time * fs, (stim_time + stim_length) * fs, facecolor='cornflowerblue', alpha=0.5)

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

# scalebar section makes x and y scales at full scale by drawing lines
x_min = axs[0].get_xlim()[0]
x_max = axs[0].get_xlim()[1]
y_min = axs[0].get_ylim()[0]
y_max = axs[0].get_ylim()[1]

# add y scalebar
axs[0].vlines(x_max - (x_max * 0.01), 0, y_scalebar)

# add x scalebar after some arithmetic
x_scalebar_start = x_min + 3 * (x_scalebar * fs)
x_scalebar_end = x_scalebar_start + (x_scalebar * fs)
axs[0].hlines(y_min, x_scalebar_start, x_scalebar_end)

# add labels
axs[0].text(0.3, 0.9, '{} pA'.format(y_scalebar * 1000), ha='center', va='center', transform=axs[0].transAxes, fontsize=8)
axs[0].text(0.7, 0.9, '{} ms'.format(x_scalebar), ha='center', va='center', transform=axs[0].transAxes, fontsize=8)

# reset axes limits to what they were before lines
axs[0].set_xlim(x_min, x_max)
axs[0].set_ylim(y_min, y_max)

# save figure to file
fig_save_path = os.path.join(figure_dir, 'JH190626_c1_attch_light_1')
fig_scales.savefig(fig_save_path + '_waves_scales.png', dpi=300, format='png')
