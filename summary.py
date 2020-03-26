import clamp_ephys
import os
import pandas as pd

paths = clamp_ephys.workflows.file_structure('server', 'Injected_GC_data/VC_pairs')
summary = pd.DataFrame()

for root, dirs, files in os.walk(os.path.join(paths.tables, 'p2'):
    for name in files:
        if "sweepavg" in name:
            path = os.path.join(root, name)
            data = pd.read_csv(path)
            summary = pd.concat([summary, data], ignore_index=True)

save_path = os.path.join(paths.tables, 'p2_summary.csv')
summary.to_csv(save_path, float_format='%8.4f', index=False)