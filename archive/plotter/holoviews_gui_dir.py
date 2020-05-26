import numpy as np
import pandas as pd
import holoviews as hv
from holoviews import opts, streams
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
#from holoviews_functions import *
import skimage.io as io
import xarray as xr
import hvplot.xarray
import panel as pn
import panel.widgets as pnw
import platform
import os
import clamp_ephys
pn.extension()
hv.extension('bokeh')



xpixels = 200
ypixels = 100
aspect_ratio = ypixels / xpixels
frames = 5000
images = np.array([np.random.rand(ypixels, xpixels) for i in range(frames)])
width = 500
height = int(width * aspect_ratio)

ds = hv.Dataset((np.arange(images.shape[0]), np.arange(images.shape[1]), np.arange(images.shape[2]),
                 images.transpose([2, 1, 0])),
                ['Time', 'y', 'x'], 'Fluorescence')

 [markdown]
# ## Play with xarray and my image


testarr = xr.DataArray(images.transpose([2, 1, 0]), dims=['x', 'y', 'time'],
                       coords=dict(x=range(images.shape[2]),
                                   y=range(images.shape[1]),
                                   time=range(images.shape[0])))



# define function for setting directories
def set_home():
    machine = platform.uname()[0]
    
    if remote_local.value == 'local':
        home_dir = os.path.expanduser("~")
        return home_dir
    
    elif remote_local.value == 'server':
        if machine == 'Darwin':
            home_dir = '/Volumes/Urban'
            return home_dir
        
        elif machine == 'Linux':
            home_dir = '/run/user/1000/gvfs/smb-share:server=130.49.237.41,share=urban'
            return home_dir
        
        elif machine == 'Windows':
            home_dir = os.path.join('Z:', os.sep, 'urban')
            return home_dir
        
        
def set_paths(event):
    home_dir = set_home()
    project_dir = os.path.join(home_dir, project_entry.value)
    project_contents = os.listdir(project_dir)
    project_contents.sort()
    
    project_path.value = project_dir
    date_selector.options = project_contents



# set up remote or local selector
remote_local = pnw.RadioBoxGroup(name='remote or local', options=['local', 'server'],
                                 inline=True, height=30, width=100)
# set up set directory button
set_paths_button = pnw.Button(name='Set Paths', width=200)
set_paths_button.on_click(set_paths)

# set up text entry string widget
project_entry = pnw.TextInput(name='Project Directory',
                              placeholder='enter path to project dir from home dir')

# set up text box to store directory path values
project_path = pnw.TextInput(name='Project Path', disabled=True)

dir_entry = pn.Column(pn.Row(remote_local, set_paths_button), project_entry, project_path)

# set up selectors
date_selector = pnw.Select(name='Choose Dir by Date', options=[])
data_selector = pnw.Select(name='Choose Data Directory', options=[])
file_selector = pnw.Select(name='Choose File', options=[])


def set_data_path(event):
    path = os.path.join(project_path.value, event.new)
    listdir = os.listdir(path)
    listdir.sort()
    data_selector.options = listdir

def set_file_path(event):
    path = os.path.join(project_path.value, date_selector.value, event.new)
    listdir = os.listdir(path)
    listdir.sort()
    file_selector.options = listdir

date_selector.param.watch(set_data_path, 'value')
data_selector.param.watch(set_file_path, 'value')

selectors = pn.Column(date_selector, data_selector, file_selector)

dir_logic = pn.Row(dir_entry, selectors)

#dir_logic.show()

 [markdown]
# # Look at actual data


images = io.imread('image.tif')
xpixels = images.shape[2]
ypixels = images.shape[1]
frames = images.shape[0]
aspect_ratio = ypixels / xpixels



#ds = hv.Dataset((np.arange(images.shape[0]), np.arange(images.shape[1]), np.arange(images.shape[2]),
#                 images.transpose([2, 1, 0])),
#                ['Time', 'y', 'x'], 'Fluorescence')

testarr = xr.DataArray(images.transpose([2, 1, 0]), dims=['x', 'y', 'time'],
                       coords=dict(x=range(images.shape[2]),
                                   y=range(images.shape[1]),
                                   time=range(images.shape[0])))



def plot_dff(bg_data, sig_data):
    
    title = 'Background subtracted Df/f for ROIs shown above'
    origin_line = hv.HLine(0)
    origin_line.opts(opts.HLine(line_color='grey', line_alpha=0.5,
                                line_width=1, line_dash='dashed'))
 
    dff = np.zeros_like(range(frames))

    if bg_data and sig_data != None:
        
        if len(bg_data['xs']) != 0:
            bg_indices = dynamic_define_roi(bg_data, xpixels, ypixels)
            sig_indices = dynamic_define_roi(sig_data, xpixels, ypixels)


            bg_roi = ds.select(x=bg_indices[1].tolist(), y=bg_indices[0].tolist())
            bg_mean = bg_roi.aggregate('Time', np.mean, spreadfn=np.std)

            sig_roi = ds.select(x=sig_indices[1].tolist(), y=sig_indices[0].tolist())
            sig_mean = sig_roi.aggregate('Time', np.mean, spreadfn=np.std)

            subtracted = sig_mean['Fluorescence'] - bg_mean['Fluorescence'].mean()
            mean_f = subtracted.mean()
            delta_f = subtracted - mean_f
            dff = delta_f / mean_f


    dff_curve = hv.Curve(data={'Time': range(len(dff)), 'Df/f': dff},
                         kdims=['Time'], vdims='Df/f').opts(title=title)

    dff_layout = dff_curve * origin_line
    return dff_layout

def save_rois(stream):
    xs = stream.data['xs'][0]
    ys = stream.data['ys'][0]
    name = stream.name
    df = pd.DataFrame(dict(xs=xs, ys=ys))
    df.to_csv('roi_{}.csv'.format(name), float_format='%8.4f', index=False)
    
    print('Saved {} ROI to csv'.format(name))
    
def save_dff(dmap):
    element = dmap.get(())
    values = element.dimension_values('Df/f')
    df = pd.DataFrame({'Df/f': values})
    df.to_csv('roi_dff.csv', float_format='%8.4f', index=False)
    
    print('Saved Df/f to csv')



# set defaults for plots
opts.defaults(
    # opts.GridSpace(shared_xaxis=True, shared_yaxis=True),
    opts.Image(cmap='Greys_r',  toolbar='above', xaxis=None, yaxis=None, height=height, width=width, ),
    opts.Labels(text_color='white', text_font_size='8pt', text_align='left', text_baseline='bottom'),
    opts.Path(color='white'),
    opts.Spread(width=600),
    opts.Overlay(show_legend=False),
    opts.Curve(width=width, height=height, toolbar='above',
               default_tools=['save, reset, pan, box_zoom, wheel_zoom'],
               framewise=True)
)



# create objects and streams for choosing ROIs
background = hv.Polygons([])
signal = hv.Polygons([])

background_stream = streams.PolyDraw(source=background,
                                     drag=True,
                                     show_vertices=True,
                                     num_objects=1, linked=True,
                                     rename={'data': 'bg_data'},
                                     name='background')
signal_stream = streams.PolyDraw(source=signal,
                                 drag=True,
                                 show_vertices=True,
                                 num_objects=1, linked=True,
                                 rename={'data': 'sig_data'},
                                 name='signal')

background.opts(opts.Polygons(fill_color='red', fill_alpha=0.5));
signal.opts(opts.Polygons(fill_color='green', fill_alpha=0.5));



# create widget pieces
player = pnw.Player(start=0, end=images.shape[0]-1,
                    value=0, aspect_ratio='auto',
                    width=520, height=100
                    #, step=10  # for faster play
                    #, interval=1 . # for fastest update
                    #, sizing_mode='stretch_width'  # to stretch to width of window
                   )

@pn.depends(player.param.value)
def player_callback(value):
    image = testarr[:, :, value]
    return image.hvplot(x='x', y='y').opts(cmap='Greys_r', colorbar=False)

image_dmap = hv.DynamicMap(player_callback)
dff_dmap = hv.DynamicMap(plot_dff, streams=[background_stream, signal_stream])

x, y = np.meshgrid(np.arange(xpixels), np.arange(ypixels)) # make a canvas with coordinates
x, y = x.flatten(), y.flatten()
points = np.vstack((x, y)).T 





image_pane = pn.panel(image_dmap * background * signal)
dff_pane = pn.panel(dff_dmap)
app = pn.Row(pn.Column(image_pane, player), dff_pane)

dir_app = pn.Column(dir_logic, app)
dir_app.show()





