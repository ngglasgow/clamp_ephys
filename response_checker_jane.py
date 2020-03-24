# -*- coding: utf-8 -*-
"""
Created 5 March 2019
epsc_peak_x.y.z.py

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

paths = file_structure('local', 'Injected_GC_data/VC_pairs')

# analysis for p2 data set
# Parameters (sampling frequency = 25 kHz, data collected as pA, test pulse is +10 mV jump
#                  and starts at t = 5 ms)
 
file_name_list, data_notes = create_data_notes('p2', 'p2_MC_TC_summary.csv', paths.p2_files)
responses_df = cell_response_checker('p2', file_name_list, paths.p2, sf=25)


# merge with data_notes and save
responses_data_notes = pd.merge(data_notes, responses_df)
responses_data_notes.to_csv(os.path.join(paths.tables, 'p2_responses_data_notes.csv'),
                            float_format='%8.4f', index=False)


# analysis for p14 data set
# Parameters (sampling frequency = 10 kHz, data collected as nA, test pulse is -5 mV jump
#                  and starts at t = 50 ms)
 
file_name_list, data_notes = create_data_notes('p14', 'p14_MC_TC_summary.csv', paths.p14_files)
responses_df = cell_response_checker('14', file_name_list, paths.p14, sf=10)


# merge with data_notes and save
responses_data_notes = pd.merge(data_notes, responses_df)
responses_data_notes.to_csv(os.path.join(paths.tables, 'p14_responses_data_notes.csv'),
                            float_format='%8.4f', index=False)
