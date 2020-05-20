from ipyfilechooser import FileChooser
from IPython.display import display
import ipywidgets as widgets
import os
import clamp_ephys
import numpy as np
from bokeh.io import output_notebook, show
from bokeh.plotting import figure
output_notebook()

''' ################### SET/CHECK THESE PARAMETERS BEFORE RUNNING ################## '''
lowpass_freq = 500  # Hz
stim_time = 500     # ms
post_stim = 250     # ms, amount of time after stimulus to look for max value
tp_start = 5        # ms, time of start of test pulse
vm_jump = 10        # mV, test pulse voltage jump
pre_tp = 3          # ms, amount of time before test pulse start to get baseline
unit_scaler = -12   # unitless, scaler to get back to A, from pA
amp_factor = 1      # scaler for making plots in pA
fs = 25             # kHz, the sampling frequency
timepoint = 'p2'
'''#################### THIS LOOP RUNS THE SINGLE CELL ANALYSIS #################### '''


data_choose = FileChooser()
data_choose.default_path = os.path.expanduser('~')
data_choose.use_dir_icons = True
data_choose.title = '<b>Select Data</b>'

data_notes_choose = FileChooser()
data_notes_choose.default_path = os.path.expanduser('~')
data_notes_choose.use_dir_icons = True
data_notes_choose.title = '<b>Select Data Notes</b>'

sweep_slider = widgets.IntSlider(value=0, min=0, max=10, step=1, continuous_update=False)

open_cell = widgets.Button(description='Open Cell', tooltip='open selected cell file with data notes')
textout = widgets.Output()
plotout = widgets.Output()

def open_cell_click(button):
    file = data_choose.selected_filename
    file_path = data_choose.selected
    data_notes = data_notes_choose.selected
    cell = clamp_ephys.workflows.cell(file_path, fs=fs, path_to_data_notes=data_notes, timepoint=timepoint, amp_factor=amp_factor)
    cell.filter_traces(lowpass_freq)
    trace = cell.traces_filtered[0]
    time = np.arange(0, len(trace) / cell.fs, 1 / cell.fs)
    ntraces = len(cell.traces_filtered.columns)
    sweep_slider.max = ntraces

    with textout:
        print(f'{file} opened')

    with plotout:
        p = figure(plot_width=800, plot_height=400)
        p.line(time, trace)
        show(p)
    

open_cell.on_click(open_cell_click)

display(data_choose, data_notes_choose, open_cell, textout, sweep_slider, plotout)

test_button = widgets.Button(icon='angle-right')
display(test_button)
