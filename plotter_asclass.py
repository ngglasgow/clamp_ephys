
# # Using the Explorer
# The explorer allows sweep exploration and saving sweeps to drop all within the browser.
# - For ease, click the double-forward-arrow icon on top toolbar to run all the cells in this notebook and display the explorer. 
# - You must have bokeh installed: `conda install bokeh`
# - You must have ipyfilechooser installed: `pip install ipyfilechooser`
# - To work in jupyter lab please follow guide here: https://docs.bokeh.org/en/latest/docs/user_guide/jupyter.html

from ipyfilechooser import FileChooser
from IPython.display import display
import ipywidgets as widgets
import os
import pandas as pd
import clamp_ephys
import numpy as np
from bokeh.io import output_notebook, show, push_notebook
from bokeh.plotting import figure
output_notebook()

class physiplot:
    '''
    Wrapper class for the physiplot explorer.
    '''
    def __init__(self):
        self.lowpass_freq_ = widgets.IntText(value=500, description='Lowpass Filter Freq (Hz):', layout=widgets.Layout(width='150px'))
        self.stim_time_ = widgets.IntText(value=500, description='Stim onset time (ms):', layout=widgets.Layout(width='150px'))
        self.post_stim_ = widgets.IntText(value=250, description='Stim window length (ms):', layout=widgets.Layout(width='150px'))
        self.tp_start_ = widgets.IntText(value=5, description='TP onset (ms):', layout=widgets.Layout(width='150px'))
        self.vm_jump_ = widgets.IntText(value=10, description='TP Vm jump (mV):', layout=widgets.Layout(width='150px'))
        self.pre_tp_ = widgets.IntText(value=3, description='TP pre window (ms):', layout=widgets.Layout(width='150px'))
        self.fs_ = widgets.IntText(value=25, description='Sampling Freq (kHz):', layout=widgets.Layout(width='150px'))
        self.amp_factor_ = widgets.IntText(value=1, description='Scaler for plots in pA', layout=widgets.Layout(width='150px'))
        self.timepoint_ = widgets.Dropdown(options=['p2', 'p14'], value='p2', description='Timepoint', layout=widgets.Layout(width='150px'))

        self.data_choose = FileChooser()
        self.data_choose.default_path = os.path.expanduser('~')
        self.data_choose.use_dir_icons = True
        self.data_choose.title = '<b>Select Data</b>'

        self.data_notes_choose = FileChooser()
        self.data_notes_choose.default_path = os.path.expanduser('~')
        self.data_notes_choose.use_dir_icons = True
        self.data_notes_choose.title = '<b>Select Data Notes</b>'

        self.sweep_slider = widgets.IntSlider(value=0, min=0, max=10, step=1, description='Sweep #: ', continuous_update=False)
        self.sweep_left = widgets.Button(icon='angle-left', tooltip='next sweep', layout=widgets.Layout(width='40px'))
        self.sweep_right = widgets.Button(icon='angle-right', tooltip='previous sweep', layout=widgets.Layout(width='40px'))

        self.open_cell = widgets.Button(description='Open Cell', tooltip='open selected cell file with data notes', layout=widgets.Layout(width='80px'))

        self.textout = widgets.Output()
        self.plotout = widgets.Output()
        self.dropout = widgets.Output()
        self.dropout = widgets.Output()

        self.drop_sweep = widgets.Button(description='Drop', tooltip='adds sweep to drop_list', layout=widgets.Layout(width='63px'))
        self.undrop_sweep = widgets.Button(description='Undrop', tooltip='removes sweep from drop_list', layout=widgets.Layout(width='63px'))
        self.save_drops = widgets.Button(icon='save', tooltip='save dropped sweeps', layout=widgets.Layout(width='32px'))
                

        metadata = [self.lowpass_freq_, self.stim_time_, self.post_stim_, self.tp_start_, self.vm_jump_, self.pre_tp_, self.fs_, self.amp_factor_, self.timepoint_]
        metadata_box = widgets.GridBox(metadata, layout=widgets.Layout(grid_template_columns='repeat(2, 160px)'))
        file_box = widgets.VBox([self.data_notes_choose, self.data_choose])
        top_box = widgets.HBox([file_box, metadata_box])
        open_box = widgets.HBox([self.open_cell, self.textout])
        sweep_box = widgets.HBox([self.sweep_slider, self.sweep_left, self.sweep_right, self.drop_sweep, self.undrop_sweep, self.save_drops, self.dropout])
        self.explorer = widgets.VBox([top_box, open_box, sweep_box, self.plotout])

    def open_cell_click(self, button):
        # set metadata from widgets
        lowpass_freq = self.lowpass_freq_.value
        stim_time = self.stim_time_.value
        post_stim = self.post_stim_.value
        tp_start = self.tp_start_.value
        vm_jump = self.vm_jump_.value
        pre_tp = self.pre_tp_.value
        fs = self.fs_.value
        amp_factor = self.amp_factor_.value
        timepoint = self.timepoint_.value

        self.drop_list = []
        self.file = self.data_choose.selected_filename
        self.file_path = self.data_choose.selected
        self.path = self.data_choose.selected_path
        self.data_notes = self.data_notes_choose.selected
        self.cell = clamp_ephys.workflows.cell(self.file_path, fs=fs, path_to_data_notes=self.data_notes, timepoint=timepoint, amp_factor=amp_factor)
        self.cell.filter_traces(lowpass_freq)
        self.trace = self.cell.traces_filtered[0]
        time = np.arange(0, len(self.trace) / self.cell.fs, 1 / self.cell.fs)
        ntraces = len(self.cell.traces_filtered.columns)
        self.sweep_slider.max = ntraces - 1
        TOOLS = ['pan, wheel_zoom, crosshair, reset']
        
        self.textout.clear_output()
        self.plotout.clear_output()
        self.dropout.clear_output()
        
        with self.textout:
            print(f'{self.file} opened')

        with self.plotout:
            self.p = figure(plot_width=800, plot_height=400, title=f'Sweep #: {self.trace.name}', toolbar_location='above',
                    tools=TOOLS, active_scroll="wheel_zoom")
            self.p.xaxis.axis_label = 'Time (ms)'
            self.p.yaxis.axis_label = 'Current (pA)'
            self.sweep = self.p.line(time, self.trace)
            show(self.p, notebook_handle=True)

    def sweep_change(self, change):
        with self.plotout:
            self.sweep.data_source.data['y'] = self.cell.traces_filtered[change['new']]
            self.p.title.text = f'Sweep #: {self.sweep_slider.value}'
            push_notebook()

    def update_drop_sweep(self, button):
        sweep_number = self.sweep_slider.value
        if sweep_number not in self.drop_list:
            self.dropout.clear_output()
            self.drop_list.append(sweep_number)
            self.drop_list.sort()
            update_string = f'Dropped sweeps: {self.drop_list}\nsweep {sweep_number} dropped'
            with self.dropout:
                print(update_string)

    def update_undrop_sweep(self, button):
        sweep_number = self.sweep_slider.value
        if sweep_number in self.drop_list:
            self.dropout.clear_output()
            self.drop_list.remove(sweep_number)
            self.drop_list.sort()
            update_string = f'Dropped sweeps: {self.drop_list}\nsweep {sweep_number} undropped'
        else:
            self.dropout.clear_output()
            update_string = f'Dropped sweeps: {self.drop_list}\nno sweeps undropped'
        with self.dropout:
            print(update_string)
                    

    def save_dropped_sweeps(self, button):
        drop_path = os.path.join(self.path, 'dropped_sweeps.csv')
        new_drops = pd.DataFrame({'dropped': [self.drop_list]}, index=[self.file])

        if os.path.exists(drop_path) == True:
            dropped_sweeps = pd.read_csv(drop_path, index_col=[0])
            if self.file in dropped_sweeps.index:
                dropped_sweeps.drop(index=self.file, inplace=True)

            dropped_sweeps = pd.concat([dropped_sweeps, new_drops])

        else:
            dropped_sweeps = new_drops

        dropped_sweeps.to_csv(drop_path)

        self.dropout.clear_output()
        update_string = f'Dropped sweeps: {self.drop_list}\nSaved in {drop_path}'
        
        with self.dropout:
            print(update_string)

    
    def sweep_left_click(self, button):
        if self.sweep_slider.value > 0:
            self.sweep_slider.value -= 1

    def sweep_right_click(self, button):
        if self.sweep_slider.value < self.sweep_slider.max:
            self.sweep_slider.value += 1
    
    def launch(self):
        self.sweep_slider.observe(self.sweep_change, names='value')
        self.drop_sweep.on_click(self.update_drop_sweep)
        self.undrop_sweep.on_click(self.update_undrop_sweep)
        self.save_drops.on_click(self.save_dropped_sweeps)

        self.sweep_left.on_click(self.sweep_left_click)
        self.sweep_right.on_click(self.sweep_right_click)
        self.open_cell.on_click(self.open_cell_click)

        return self.explorer
