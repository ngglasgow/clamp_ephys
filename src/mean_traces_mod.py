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

''' ################## Define file structure on server #################### '''
# home_dir will depend on the OS, but the rest will not
# query machine identity and set home_dir from there
machine = os.uname()[0]
if machine == 'Darwin':
    home_dir = '/Volumes/Urban'

elif machine == 'Linux':
    home_dir = '/run/user/1000/gvfs/smb-share:server=130.49.237.41,share=urban'

else:
    home_dir = 'C:\\Users\\Jane Huang'

project_dir = os.path.join(home_dir, 'Documents', 'Grad school', 'Urban Lab', 'Data Analysis', 'OSN_OMPvGg8_MTC')
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

# pull out MC-only Gg8 cells
noigor_list = np.array(data_notes[data_notes['Igor saved?'] == 'No'].index)
data_notes = data_notes.drop(noigor_list)

data_notes = pd.read_csv(os.path.join(table_dir, 'OSN_Gg8vOMP.csv'))
MC_only = data_notes['Cell type'] == 'MC'
Gg8_only = data_notes['Genotype'] == 'Gg8'

MC_list = data_notes[data_notes['Cell type'] == 'MC']
MC_Gg8_list = MC_list[MC_list['Genotype'] == 'Gg8']


# create empty dataframes for gg8 and omp cells
gg8_cells = pd.DataFrame()
omp_cells = pd.DataFrame()

gg8_example_ctrl = pd.DataFrame()
gg8_example_drug = pd.DataFrame()

# loop through all files, extract data and add to appropriate dataframes

for file in file_list:
    if 'timeseries' in file:
        # open file and extract data into a new dataframe
        mean_trace = pd.read_csv(os.path.join(table_dir, file), header=None)

        # separate out gg8 and omp genotype cells into dataframes
        if 'Gg8' in file:
            if 'c1a' in file: # drops c1a from Gg8 list because it's just ctrl
                continue
            else:                
                gg8_cells = pd.concat([gg8_cells, mean_trace], axis=1, ignore_index=True)
            
        else:
            omp_cells = pd.concat([omp_cells, mean_trace], axis=1, ignore_index=True)

    else:
        continue
    
# Make separate time series for Gg8 example MC cell control and drug traces
    
for file in file_list:
    if 'c1a_mean_timeseries' in file:
        mean_trace = pd.read_csv(os.path.join(table_dir, file), header=None)
        gg8_example_ctrl = pd.concat([gg8_example_ctrl, mean_trace], axis=1, ignore_index=True)
        
    elif 'JH181120_c1_mean_timeseries' in file:
        mean_trace = pd.read_csv(os.path.join(table_dir, file), header=None)
        gg8_example_drug = pd.concat([gg8_example_drug, mean_trace], axis=1, ignore_index=True)

    else:
        continue

# calculate means of gg8 and omp traces
gg8_mean = gg8_cells.mean(axis=1)
omp_mean = omp_cells.mean(axis=1)

# calculate means of example Gg8 traces
gg8_example_ctrl = gg8_example_ctrl.mean(axis=1)
gg8_example_drug = gg8_example_drug.mean(axis=1)


'''########## make plots of mean and all traces of omp and gg8 ############'''
sweep_length = len(gg8_mean)             # allow for different sweep length
sweep_time = np.arange(0, sweep_length/10, 0.1)     # time of sweeps in ms

# sweep length and times for Gg8 example cell - 3s sweep
sweep_length_example = len(gg8_example_ctrl)             # allow for different sweep length
sweep_time_example = np.arange(0, sweep_length_example/10, 0.1)     # time of sweeps in ms

# make a figure with 2 plots
fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
# fig.suptitle('')

# plot mean data trace with all traces in gray behind
''' commented out OMP graph
axs[0].plot(sweep_time, omp_cells*1000, color='darkgray', linewidth=0.5)
axs[0].plot(sweep_time, omp_mean*1000, color='k')
axs[0].hlines(75, 500, 550, color='deepskyblue')
axs[0].set_title('OMP')
axs[0].set_xlabel('Time (ms)')
axs[0].set_ylabel('Current (pA)')
axs[0].set_xlim(450, 800)
axs[0].set_ylim(-1050, 100)
'''
#first graph is Gg8 all averages for MCs only
axs[0].plot(sweep_time, gg8_cells*1000, color='darkgray', linewidth=0.5)
axs[0].plot(sweep_time, gg8_mean*1000, color='k')
axs[0].hlines(75, 500, 550, color='deepskyblue')
axs[0].set_title('Gg8')
axs[0].set_xlabel('Time (ms)')
axs[0].set_ylabel('Current (pA)')
axs[0].set_xlim(450, 800)
axs[0].set_ylim(-400, 100)

#second graph is example Gg8 cell
axs[1].plot(sweep_time_example, gg8_example_ctrl*1000, color='k', label='Ctrl')
axs[1].plot(sweep_time_example, gg8_example_drug*1000, color='red', label='4-AP + TTX')
axs[1].hlines(75, 500, 550, color='deepskyblue')
axs[1].set_title('Gg8')
axs[1].set_xlabel('Time (ms)')
axs[1].set_ylabel('Current (pA)')
axs[1].set_xlim(450, 800)
axs[1].set_ylim(-400, 100)

fig

# save figure to file
fig_save_path = os.path.join(figure_dir, 'OMPvGg8_meantraces.png')
fig.savefig(fig_save_path, dpi=300, format='png')
