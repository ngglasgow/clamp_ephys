import clamp_ephys
import os
import pandas as pd
import matplotlib.pyplot as plt

'''####################### SET THE PROPER PATH YOU WANT ########################### '''
paths = clamp_ephys.workflows.file_structure('local', 'Injected_GC_data/VC_pairs')
p2 = pd.DataFrame()
tables = paths.tables
figures = paths.figures

for root, dirs, files in os.walk(os.path.join(paths.tables, 'p2')):
    for name in files:
        if "mean_subtracted_timeseries" in name:
            path = os.path.join(root, name)
            data = pd.read_csv(path)
            filename = '_'.join(name.split('_')[0:3])
            data.columns = [filename]
            p2 = pd.concat([p2, data], axis=1)

save_path = os.path.join(paths.tables, 'p2_mean_traces.csv')
p2.to_csv(save_path, float_format='%8.4f', index=False)

''' create a DataFrame for light conditions '''
p2_light = p2.filter(regex='light')

''' drop any wanted/irrelevant cells, only show yes respoonses '''
p2_light_MCs = p2_light.drop(columns=['JH200303_c2_light100', 'JH200303_c4_light100', 'JH200303_c6_light100',
                                        'JH200303_c8_light100', 'JH200304_c2_light100', 'JH200304_c3_light100',
                                        'JH200310_c2_light100', 'JH200310_c4_light100', 'JH200311_c2_light100',
                                        'JH200311_c5_light100', 'JH200313_c1_light100', 'JH200313_c5_light100'])

# drop cells without visual responses
p2_light_MCs_responses = p2_light_MCs.drop(columns=['JH200303_c1_light100', 'JH200303_c3_light100', 'JH200304_c1_light100',
                                                    'JH200304_c4_light100', 'JH200310_c1_light100', 'JH200310_c3_light100',
                                                    'JH200313_c4_light100', 'JH200313_c6_light100'])



p2_light_TCs = p2_light.drop(columns=['JH200303_c1_light100', 'JH200303_c3_light100', 'JH200303_c7_light100',
                                        'JH200304_c1_light100', 'JH200304_c4_light100', 'JH200310_c1_light100',
                                        'JH200310_c3_light100', 'JH200311_c1_light100', 'JH200311_c3_light100',
                                        'JH200311_c4_light100', 'JH200311_c6_light100', 'JH200313_c2_light100',
                                        'JH200313_c3_light100', 'JH200313_c4_light100', 'JH200313_c6_light100'])



''' ########### stacking figures with no frame or scales ################## '''
# first two lines are for if you want to make plots have unequal dimensions
# gs_kw = dict(height_ratios=[0.5, 0.7, 0.7, 0.7, 4, 1])
# fig, axs = plt.subplots(6, 8, figsize=(20, 5.5), gridspec_kw=gs_kw, constrained_layout=True)


''' some parameters for automating figure sizing and rows '''
n_cells = len(p2_light_TCs.columns)

# rename columns for indexing
p2_light_TCs.columns = range(n_cells)

width = 2                   # plot width in inches
height = n_cells * 0.5      # height in inches of figure, scaler is inch per plot
stim_time = 500             # time of stimulus onset in ms
stim_length = 100           # length of stimulus in ms
fs = 25                     # sampling frequency in kHz
light_amplitude = 0.045     # plot amp in nA where you want blue bar
plotx_start = 450 * fs      # index of time axis you want to start plotting
plotx_end = 1500 * fs       # index of time axis you want to stop plotting
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
fig, axs = plt.subplots(n_cells, 1, figsize=(width, height), constrained_layout=True)

# pull out one wave, subtract baseline, take slice, and plot
for i in range(n_cells):
    # baseline subtract for each wave
    wave = p2_light_TCs[i]
    
    # take slice of baseline subtracted wave and calculate max/min for scaling
    wave_slice = wave[plotx_start:plotx_end]

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
fig

filename = 'p2_TC_simplifed_traces_noaxes.png'
path = r"C:\Users\jhuang\Documents\phd_projects\Injected_GC_data\VC_pairs\figures\p2"
fig.savefig(path + filename, dpi=300, format='png')
plt.close()


# generate the figure with all the subplots
fig, axs = plt.subplots(n_cells, 1, figsize=(width, height), constrained_layout=True)

# pull out one wave, subtract baseline, take slice, and plot
for i in range(n_cells):
    # baseline subtract for each wave
    wave = p2_light_TCs[i]
    
    # take slice of baseline subtracted wave and calculate max/min for scaling
    wave_slice = wave[plotx_start:plotx_end]

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

fig

filename = 'p2_TC_simplifed_traces_axes.png'
path = r"C:\Users\jhuang\Documents\phd_projects\Injected_GC_data\VC_pairs\figures\p2"
fig.savefig(path + filename, dpi=300, format='png')
plt.close()