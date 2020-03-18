# -*- coding: utf-8 -*-
"""
Created 23 April 2019
mean_traces.py
Version 1
The purpose of this script is to pull all of the mean trace files that were
saved from the initial analysis. These traces are mean subtracted and filtered
and comprise the entire 6 s of recording. The idea here is to open the files
individually, extract the data, save it to a dataframe and compile all of the
files of the same genotype into a dataframe. Then take the mean. Then plot the
means vs. all traces for both OMP and Gg8.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

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


''' ##########################################################################
This is all the analysis, figures, saving
Read in file metadata, open file from igor, convert to pandas
##############################################################################
'''
# grab all files in table_dir
file_list = os.listdir(table_dir)
trace_files = []
cell_ids = []

for file in file_list:
    if 'timeseries' in file:
        trace_files.append(file)
        cell_id = file.split('_')[1] + '_' + file.split('_')[2]
        cell_ids.append(cell_id)

    else:
        continue

traces_df = pd.DataFrame({'file name': trace_files, 'cell id': cell_ids})

# grab data_notes to select out by cell type
analyzed_data_notes = pd.read_csv(os.path.join(table_dir ,'analyzed_data_notes.csv'), index_col=0)
mc_df = analyzed_data_notes[analyzed_data_notes['Cell type'] == 'MC']

# pull out gg8 cells
mc_gg8_df = mc_df[mc_df['Genotype'] == 'Gg8']
mc_gg8_list = mc_gg8_df['Cell name'].tolist()
mc_gg8_list = [name.split('_')[0] + '_' + name.split('_')[1] for name in mc_gg8_list]
mc_gg8_df = pd.DataFrame(mc_gg8_list, columns=['cell id'])
type(mc_gg8_df['Cell name'])
# pull out omp cells
mc_omp_df = mc_df[mc_df['Genotype'] == 'OMP']
mc_omp_list = mc_omp_df['Cell name'].set_minor_locatorist()
mc_omp_list = [name.split('_')[0] + '_' + name.split('_')[1] for name in mc_omp_list]
mc_omp_df = pd.DataFrame(mc_omp_list, columns=['cell id'])

# make list of Gg8 MCs
gg8_mcs = pd.merge(traces_df, mc_gg8_df)
gg8_mc_list = gg8_mcs['file name'].tolist()

# make list of OMP MCs
omp_mcs = pd.merge(traces_df, mc_omp_df)
omp_mc_list = omp_mcs['file name'].tolist()

# create empty dataframes for gg8 and omp cells
gg8_cells = pd.DataFrame()
omp_cells = pd.DataFrame()

# loop through all files, extract data and add to appropriate dataframes
for file in gg8_mc_list:
    # open file and extract data into a new dataframe
    mean_trace = pd.read_csv(os.path.join(table_dir, file), header=None)
    gg8_cells = pd.concat([gg8_cells, mean_trace], axis=1, ignore_index=True)

for file in omp_mc_list:
    # open file and extract data into a new dataframe
    mean_trace = pd.read_csv(os.path.join(table_dir, file), header=None)
    omp_cells = pd.concat([omp_cells, mean_trace], axis=1, ignore_index=True)

# calculate means of gg8 and omp traces
gg8_mean = gg8_cells.mean(axis=1)
omp_mean = omp_cells.mean(axis=1)


'''########## make plots of mean and all traces of omp and gg8 ############'''
sweep_length = len(gg8_mean)             # allow for different sweep length
sweep_time = np.arange(0, sweep_length/10, 0.1)     # time of sweeps in ms
y_lim_low = -200
# make a figure with 2 plots
fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
# fig.suptitle('')

# plot mean data trace with all traces in gray behind
axs[0].plot(sweep_time, omp_cells*1000, color='darkgray', linewidth=0.5)
axs[0].plot(sweep_time, omp_mean*1000, color='k')
axs[0].hlines(75, 500, 550, color='deepskyblue')
axs[0].set_title('OMP')
axs[0].set_xlabel('Time (ms)')
axs[0].set_ylabel('Current (pA)')
axs[0].set_xlim(450, 800)
axs[0].set_ylim(y_lim_low, 100)

axs[1].plot(sweep_time, gg8_cells*1000, color='darkgray', linewidth=0.5)
axs[1].plot(sweep_time, gg8_mean*1000, color='k')
axs[1].hlines(75, 500, 550, color='deepskyblue')
axs[1].set_title('Gg8')
axs[1].set_xlabel('Time (ms)')
axs[1].set_ylabel('Current (pA)')
axs[1].set_xlim(450, 800)
axs[1].set_ylim(y_lim_low, 100)

fig

# save figure to file
fig_save_path = os.path.join(figure_dir, 'OMPvGg8_meantraces.png')
fig.savefig(fig_save_path, dpi=300, format='png')
