import clamp_ephys

path = "/home/nate/urban/clamp_ephys/test_data/control.ibw"
data = clamp_ephys.cell.data(path, 10)

data.sf