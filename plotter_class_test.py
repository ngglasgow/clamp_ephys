from ipyfilechooser import FileChooser
from IPython.display import display
import ipywidgets as widgets
import os
import clamp_ephys
import numpy as np
from bokeh.io import output_notebook, show, push_notebook
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
home_dir = os.path.expanduser('~')
data_path = os.path.join(home_dir, 'urban', 'clamp_ephys', 'test_data')

data_choose = FileChooser()
data_choose.default_path = data_path
data_choose.use_dir_icons = True
data_choose.title = '<b>Select Data</b>'

data_notes_choose = FileChooser()
data_notes_choose.default_path = data_path
data_notes_choose.use_dir_icons = True
data_notes_choose.title = '<b>Select Data Notes</b>'

sweep_slider = widgets.IntSlider(value=0, min=0, max=10, step=1, description='Sweep #: ', continuous_update=False)
sweep_left = widgets.Button(icon='angle-left', tooltip='next sweep', layout=widgets.Layout(width='40px'))
sweep_right = widgets.Button(icon='angle-right', tooltip='previous sweep', layout=widgets.Layout(width='40px'))

open_cell = widgets.Button(description='Open Cell', tooltip='open selected cell file with data notes', layout=widgets.Layout(width='80px'))

textout = widgets.Output()
plotout = widgets.Output()
dropout = widgets.Output()

drop_sweep = widgets.Button(description='Drop', tooltip='adds sweep to drop_list', layout=widgets.Layout(width='63px'))
undrop_sweep = widgets.Button(description='Undrop', tooltip='removes sweep from drop_list', layout=widgets.Layout(width='63px'))
save_drops = widgets.Button(icon='save', tooltip='save dropped sweeps', layout=widgets.Layout(width='32px'))

def update_drop_sweep(button):
    sweep_number = sweep_slider.value
    if sweep_number not in drop_list:
        dropout.clear_output()
        drop_list.append(sweep_number)
        drop_list.sort()
        update_string = f'Dropped sweeps: {drop_list}\nsweep {sweep_number} dropped'
        with dropout:
            print(update_string)

def update_undrop_sweep(button):
    sweep_number = sweep_slider.value
    if sweep_number in drop_list:
        dropout.clear_output()
        drop_list.remove(sweep_number)
        drop_list.sort()
        update_string = f'Dropped sweeps: {drop_list}\nsweep {sweep_number} undropped'
        with dropout:
            print(update_string)

drop_sweep.on_click(update_drop_sweep)
undrop_sweep.on_click(update_undrop_sweep)
dropout = widgets.Output()

def open_cell_click(button):
    file = data_choose.selected_filename
    file_path = data_choose.selected
    data_notes = data_notes_choose.selected
    cell = clamp_ephys.workflows.cell(file_path, fs=fs, path_to_data_notes=data_notes, timepoint=timepoint, amp_factor=amp_factor)
    cell.filter_traces(lowpass_freq)
    trace = cell.traces_filtered[0]
    time = np.arange(0, len(trace) / cell.fs, 1 / cell.fs)
    ntraces = len(cell.traces_filtered.columns)
    sweep_slider.max = ntraces - 1
    
    textout.clear_output()
    plotout.clear_output()

    with textout:
        print(f'{file} opened')

    with plotout:
        p = figure(plot_width=800, plot_height=400, title=f'Sweep #: {sweep_slider.value}', toolbar_location='above')
        p.xaxis.axis_label = 'Time (ms)'
        p.yaxis.axis_label = 'Current (pA)'
        sweep = p.line(time, trace)
        show(p, notebook_handle=True)

    def sweep_change(change):
        with plotout:
            sweep.data_source.data['y'] = cell.traces_filtered[change['new']]
            p.title.text = f'Sweep #: {sweep_slider.value}'
            push_notebook()
            
    def save_dropped_sweeps(button):
        '''idea is to save a single dropped_sweeps.csv for a directory that will be read in by the call to cells if drop=true
         -need to figure out how to drop cols from a pandas read csv
         - need to change this function into a class and do testing'''

    sweep_slider.observe(sweep_change, names='value')
    
    
def sweep_left_click(button):
    if sweep_slider.value > 0:
        sweep_slider.value -= 1

def sweep_right_click(button):
    if sweep_slider.value < sweep_slider.max:
        sweep_slider.value += 1
        


sweep_left.on_click(sweep_left_click)
sweep_right.on_click(sweep_right_click)
open_cell.on_click(open_cell_click)

open_box = widgets.HBox([open_cell, textout])
sweep_box = widgets.HBox([sweep_slider, sweep_left, sweep_right, drop_sweep, undrop_sweep, save_drops, dropout])
widgets.VBox([data_notes_choose, data_choose, open_box, sweep_box, plotout])