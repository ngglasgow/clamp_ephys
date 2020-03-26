import clamp_ephys

path = "/home/nate/urban/clamp_ephys/test_data/JH200303_c1_light100.ibw"
data_path = "/home/nate/urban/clamp_ephys/test_data/p2_data_notes.csv"
tables = '/home/nate/urban/clamp_ephys/test_data/'
figures = '/home/nate/urban/clamp_ephys/test_data/'

lowpass_freq = 500  # Hz
stim_time = 500     # ms
post_stim = 250     # ms, amount of time after stimulus to look for max value
tp_start = 5        # ms, time of start of test pulse
vm_jump = 10        # mV, test pulse voltage jump
pre_tp = 3          # ms, amount of time before test pulse start to get baseline
unit_scaler = -12   # unitless, scaler to get back to A, from pA
amp_factor = 1      # scaler for making plots in pA

data = clamp_ephys.workflows.cell(path, fs=25, path_to_data_notes=data_path, timepoint='p2')

data.get_raw_peaks(stim_time, post_stim)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)
data.get_series_resistance(tp_start, vm_jump, pre_tp, unit_scaler)
data.get_sweep_data()
data.get_responses(threshold=3.789)

fig = data.plot_peaks_rs(amp_factor)

data.save_fig(figures, fig)
data.save_metadata(tables)
data.save_mean_data(tables)
data.save_mean_filtered_trace(tables)
