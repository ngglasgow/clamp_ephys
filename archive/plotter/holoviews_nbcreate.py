'''
2019-08-28
Use papermll to loop through a single days worth of data and generate default
ROI background and cell fluorescence traces that I then modify
'''
import papermill as pm
import os

date = '190912'
# set paths
machine = os.uname()[0]
if machine == 'Darwin':
    home_dir = '/Volumes/Urban'

elif machine == 'Linux':
    home_dir = '/run/user/1000/gvfs/smb-share:server=130.49.237.41,share=urban'

else:
    home_dir = os.path.join('Z:', os.sep)

project_dir = os.path.join(home_dir, 'Glasgow', 'vm_imaging')
date_dir = os.path.join(project_dir, date)

# empty lists for cells and data directories
data_dir_list = [d for d in os.listdir(date_dir) if os.path.isdir(os.path.join(date_dir, d))]

input_path = '/home/ngg1/ngglasgow@gmail.com/Data_Urban/jupyter_interactive/holoviews_gui_template.ipynb'

for dir in data_dir_list[1:]:
    data_path = os.path.join(date_dir, dir)
    output = '/home/ngg1/ngglasgow@gmail.com/Data_Urban/jupyter_interactive/notebooks/{}.ipynb'.format(dir)
    pm.execute_notebook(input_path, output, parameters=dict(data_dir=dir))



data_path
output

data_dir_list

