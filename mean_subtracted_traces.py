import clamp_ephys
import os
import pandas as pd

paths = clamp_ephys.workflows.file_structure('server', 'Injected_GC_data/VC_pairs')
p2 = pd.DataFrame()
p14 = pd.DataFrame()

for root, dirs, files in os.walk(os.path.join(paths.tables, 'p2')):
    for name in files:
        if "mean_subtracted_timeseries" in name:
            path = os.path.join(root, name)
            data = pd.read_csv(path)
            filename = name.split('_')[0]
            data.columns = [filename]
            p2 = pd.concat([p2, data], axis=1)

save_path = os.path.join(paths.tables, 'p2_mean_traces.csv')
p2.to_csv(save_path, float_format='%8.4f', index=False)

for root, dirs, files in os.walk(os.path.join(paths.tables, 'p14')):
    for name in files:
        if "mean_subtracted_timeseries" in name:
            path = os.path.join(root, name)
            data = pd.read_csv(path)
            filename = name.split('_')[0]
            data.columns = [filename]
            p14 = pd.concat([p14, data], axis=1)

save_path = os.path.join(paths.tables, 'p14_mean_traces.csv')
p14.to_csv(save_path, float_format='%8.4f', index=False)
