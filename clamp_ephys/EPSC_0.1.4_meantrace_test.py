# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

"""
# from __main__ import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import OrderedDict
import math
import platform


''' ################## Define file structure on server #################### '''
# home_dir will depend on the OS, but the rest will not
# query machine identity and set home_dir from there
machine = platform.uname()[0]

if machine == 'Darwin':
    home_dir = '/Volumes/Urban'

elif machine == 'Linux':
    home_dir = '/run/user/1000/gvfs/smb-share:server=130.49.237.41,share=urban'

elif machine == 'Windows':
    home_dir = os.path.join('N:', os.sep, 'urban')

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
data_notes.to_csv(os.path.join(table_dir, 'data_notes_test.csv'))
# drop cells that didn't save to igor
noigor_list = np.array(data_notes[data_notes['Igor saved?'] == 'No'].index)
data_notes = data_notes.drop(noigor_list)

# drop cells that don't have any # of drug sweeps
nodrug_list = np.array(data_notes[data_notes['# of drug sweeps'].isnull() == True].index)
data_notes = data_notes.drop(nodrug_list)

# update file name list to have only files you want to analyze after logic
file_name_list = data_notes['Cell name'].tolist()
data_notes.to_csv(os.path.join(table_dir, 'analyzed_data_notes_test.csv'))
