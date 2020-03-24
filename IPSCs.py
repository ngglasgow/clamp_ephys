# -*- coding: utf-8 -*-
"""
Created 14 November 2019
ipsc_peak_x.y.z.py

Modified from EPSC_0.1.6_meantrace.py

"""
# from __main__ import *


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style(rc = {'figure.facecolor':'white'})
import elephant
from neo.io import IgorIO
import os
import platform
from collections import OrderedDict
import math






''' *********************************************************************** '''


''' ################## Define file structure on server #################### '''

paths = file_structure('local', 'Injected_GC_data/VC_pairs')

# analysis for p2 data set
# Parameters (sampling frequency = 25 kHz, data collected as pA, test pulse is +10 mV jump
#                  and starts at t = 5 ms)
 
file_name_list, data_notes = create_data_notes('p2', 'p2_MC_TC_summary.csv', paths.p2_files)
for file in file_name_list:
    indiv_cell_analysis('p2', file, paths.p2, sf=25, amp_factor=1, peak_factor=-12,
                         tp_start=5, vm_jump=10, pre_tp=3)

# analysis for p14 data set
# p14 parameters (sampling frequency = 10 kHz, data collected as nA, test pulse is -5 mV jump
#                  and starts at t = 50 ms)

file_name_list, data_notes = create_data_notes('p14', 'p14_MC_TC_summary.csv', paths.p14_files)
for file in file_name_list:
    indiv_cell_analysis('p14', file, paths.p14, sf=10, amp_factor=1000, peak_factor=-9,
                         tp_start=50, vm_jump=-5, pre_tp=11)
