import clamp_ephys
import pandas as pd
import os
import scipy
import matplotlib
%matplotlib
import matplotlib.pyplot as plt
import numpy as np

lowpass_freq = 500      # Hz
stim_time = 500         # ms
post_stim = 250         # ms, amount of time after stimulus to look for max value
baseline_start = 3000   # ms, time in the sweep to start taking the baseline
baseline_end = 6000     # ms, time in the sweep at which baseline ends
tp_start = 5            # ms, time of start of test pulse
vm_jump = 10            # mV, test pulse voltage jump
pre_tp = 3              # ms, amount of time before test pulse start to get baseline
unit_scaler = -12       # unitless, scaler to get back to A, from pA
amp_factor = 1          # scaler for making plots in pA
fs = 25                 # kHz, the sampling frequency
width = 3               # ms, the required width for classifying an event as a peak
timepoint = 'p2'

# paths to clamp_ephys/test_data/
project_path = os.path.join(os.path.expanduser('~'), 'Documents', 'phd_projects', 'Injected_GC_data')
notes_path = os.path.join(project_path, 'VC_pairs', 'tables', 'p2_data_notes.csv')
cell_path = os.path.join(project_path, 'VC_pairs', 'data', 'p2', 'JH200311_c5_light100.ibw')

data = clamp_ephys.workflows.cell(cell_path, path_to_data_notes=notes_path, timepoint=timepoint)

data.get_raw_peaks(stim_time, post_stim)
data.filter_traces(lowpass_freq)
data.get_filtered_peaks(stim_time, post_stim)

# data.get_peaks_widths(stim_time, width)   

data.sweeps = data.traces_filtered - data.baseline_filtered
# data.sweeps = pd.DataFrame(data.mean_traces_filtered - data.mean_baseline_filtered)

window_start = (stim_time + 20) * fs
all_peaks_kinetics_df = pd.DataFrame()
all_peaks_kinetics_avg_df = pd.DataFrame()
first3_kinetics_avg_df = pd.DataFrame()



#%% this runs on an individual sweep
sweep = 7
peak_number = 0     


if data.mean_peak_filtered > 0:
    invert = 1
else:
    invert = -1

window_start = (stim_time + 20) * data.fs
all_widths_df = pd.DataFrame()
data.all_widths_df = pd.DataFrame()
test_amp_df = pd.DataFrame()
     
trace = data.sweeps.iloc[window_start:, sweep].values
thresh = 2.5 * trace.std()

peaks, properties = scipy.signal.find_peaks(trace*invert, distance=0.5*data.fs, prominence=thresh, width=width*data.fs)

fig = plt.figure()
fig.suptitle('Sweep {}'.format(sweep))
plt.plot(trace)
plt.plot(peaks, trace[peaks], 'x')

if len(peaks) > 0:
    # calculate 10 to 90% and FWHM
    prominence_data = list(properties.values())[0:3]
    ten_widths, ten_height, ten_left, ten_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=0.9, prominence_data=prominence_data)
    ninety_widths, ninety_height, ninety_left, ninety_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=0.1, prominence_data=prominence_data)
    half_widths, hw_height, hw_left, hw_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=0.5, prominence_data=prominence_data)
    full_widths, fw_height, fw_left, fw_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=1, prominence_data=prominence_data)

    hw_height = hw_height * invert
    fw_height = fw_height * invert
    ten_height = ten_height * invert
    ninety_height = ninety_height * invert
    prominences = properties['prominences'] * invert

    peak_numbers = range(len(peaks))

    all_widths_data = pd.DataFrame({'sweep #': sweep, 'peak #': peak_numbers, 'peaks_index': peaks, 'prominence': prominences,
        'ten_widths': ten_widths, 'ten_height': ten_height, 'ten_left': ten_left, 'ten_right': ten_right, 
        'ninety_widths': ninety_widths, 'ninety_height': ninety_height, 'ninety_left': ninety_left, 
        'ninety_right': ninety_right, 'half_widths': half_widths, 'half_height': hw_height, 
        'half_left': hw_left, 'half_right': hw_right, 'full_widths': full_widths, 'full_height': fw_height, 
        'full_left': fw_left, 'full_right': fw_right})
    
    all_widths_df = pd.concat([all_widths_df, all_widths_data], ignore_index=True)

    # peak_values = trace[peaks]
    # peak_height_est = peak_values - prominences

    # test_amp_data = pd.DataFrame({'sweep #': sweep, 'peak #': peak_numbers, 'peaks_index': peaks,
    #                                 'prominences': prominences, 'peak values': peak_values, 
    #                                 'peak_height': peak_height_est})

    # plt.plot(trace)
    # plt.plot(peaks, trace[peaks], "x")
    # plt.vlines(x=peaks, ymin=peak_values - prominences,
    #         ymax = trace[peaks], color = "C1")
    # plt.show()

else:
    print('No peaks in {} sweep {}'.format(data.file_id, sweep))
    
if len(all_widths_df) > 0:
    data.all_widths_df = all_widths_df.set_index(['sweep #', 'peak #'], inplace=False)
else:
    print('No peaks in {}'.format(data.file_id))

try:
    data.all_widths_df
except:
    print('No peaks in {}'.format(data.file_id))
else:
    print('There are peaks in {}'.format(data.file_id))

# this only pulls out sweeps that had peaks   
trace = data.sweeps.iloc[window_start:, sweep].values

ninety_left = data.all_widths_df.loc[(sweep), 'ninety_left']
ten_left = data.all_widths_df.loc[(sweep), 'ten_left']
ten_to_ninety = (ninety_left - ten_left) / fs
peak_time = data.all_widths_df.loc[(sweep), 'peaks_index'] / fs
hw_time = data.all_widths_df.loc[(sweep), 'half_widths'] / fs
peak_amplitude = trace[peaks]

peak_numbers = range(len(data.all_widths_df.loc[sweep]))

tau_list = []
charge_list = []

for peak_number in peak_numbers:        
    tau = data.get_tau(trace, sweep, peak_number)
    tau_list.append(tau)
    charge = data.get_charge_transf(trace, sweep, peak_number)
    charge_list.append(charge) 

all_peaks_kinetics_data = pd.DataFrame({'sweep #': sweep,
                                        'peak #': peak_numbers,
                                        'peak time (ms)': peak_time,
                                        'peak amplitude (pA)': prominences,
                                        '10 to 90% RT (ms)': ten_to_ninety,
                                        'tau': tau_list,
                                        'half-width (ms)': hw_time,
                                        'charge transferred (pA * s)': charge_list})

all_peaks_kinetics_data = all_peaks_kinetics_data[all_peaks_kinetics_data.tau < 500] # drops peaks with tau values over 500


data.plot_first3_events(all_peaks_kinetics_data, sweep, trace)
delay_to_response = all_peaks_kinetics_data.iloc[0, 2] # gets the time of the first peak

# calculates the average kinetic values for all peaks in the given sweep 
all_peaks_kinetics_avg_data = pd.DataFrame(all_peaks_kinetics_data.mean(axis=0)).T
all_peaks_kinetics_avg_data['sweep #'] = all_peaks_kinetics_avg_data['sweep #'].astype(int) # makes sweep # an int
all_peaks_kinetics_avg_data.drop(['peak #', 'peak time (ms)'], axis=1, inplace=True) #drop peak # and peak time columns
all_peaks_kinetics_avg_data.insert(1, 'delay_to_response (ms)', delay_to_response)

# calculates the average kinetics values for the first three peaks in a sweep   
first3_kinetics_avg_data = pd.DataFrame(all_peaks_kinetics_data.iloc[:3].mean(axis=0)).T
first3_kinetics_avg_data['sweep #'] = first3_kinetics_avg_data['sweep #'].astype(int) # makes sweep # an int
first3_kinetics_avg_data.drop(['peak #', 'peak time (ms)'], axis=1, inplace=True) #drop peak # and peak time columns
first3_kinetics_avg_data.insert(1, 'delay_to_response (ms)', delay_to_response)

all_peaks_kinetics_df = pd.concat([all_peaks_kinetics_df, all_peaks_kinetics_data], ignore_index=True)
all_peaks_kinetics_avg_df = pd.concat([all_peaks_kinetics_avg_df, all_peaks_kinetics_avg_data], ignore_index=True)
first3_kinetics_avg_df = pd.concat([first3_kinetics_avg_df, first3_kinetics_avg_data], ignore_index=True)
    
# this df contains all kinetics values for all the peaks in all the sweeps
data.all_peaks_kinetics_df = all_peaks_kinetics_df.set_index(['sweep #', 'peak #'], inplace=False)

# this df contains avg kinetic values for all the peaks in all the sweeps
data.all_peaks_kinetics_avg_df = all_peaks_kinetics_avg_df.set_index(['sweep #'], inplace=False)

# this df contains avg kinetic values for the first three peaks in all the sweeps
data.first3_kinetics_avg_df = first3_kinetics_avg_df.set_index(['sweep #'], inplace=False)

# this df calculates the average of averages for the first three peaks in all the sweeps
data.avg_first3_kinetics_avg_df = pd.DataFrame(first3_kinetics_avg_df.mean(axis=0)).T        
data.avg_first3_kinetics_avg_df = pd.concat([data.metadata, data.avg_first3_kinetics_avg_df], axis=1)
data.avg_first3_kinetics_avg_df.drop(['sweep #'], axis=1, inplace=True)



# %% this runs on all the sweeps
if data.mean_peak_filtered > 0:
    invert = 1
else:
    invert = -1

window_start = (stim_time + 5) * data.fs
all_widths_df = pd.DataFrame()
data.all_widths_df = pd.DataFrame()

for sweep in range(len(data.sweeps.columns)):        
    trace = data.sweeps.iloc[window_start:, sweep].values
    thresh = 2.5 * trace.std()

    peaks, properties = scipy.signal.find_peaks(trace*invert, distance=0.5*data.fs, prominence=thresh, width=width*data.fs)

    if len(peaks) > 0:
        # calculate 10 to 90% and FWHM
        prominence_data = list(properties.values())[0:3]
        ten_widths, ten_height, ten_left, ten_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=0.9, prominence_data=prominence_data)
        ninety_widths, ninety_height, ninety_left, ninety_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=0.1, prominence_data=prominence_data)
        half_widths, hw_height, hw_left, hw_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=0.5, prominence_data=prominence_data)
        full_widths, fw_height, fw_left, fw_right = scipy.signal.peak_widths(trace*invert, peaks, rel_height=1, prominence_data=prominence_data)

        hw_height = hw_height * invert
        fw_height = fw_height * invert
        ten_height = ten_height * invert
        ninety_height = ninety_height * invert
        prominences = properties['prominences']

        peak_numbers = range(len(peaks))

        all_widths_data = pd.DataFrame({'sweep #': sweep, 'peak #': peak_numbers, 'peaks_index': peaks, 'prominence': prominences,
            'ten_widths': ten_widths, 'ten_height': ten_height, 'ten_left': ten_left, 'ten_right': ten_right, 
            'ninety_widths': ninety_widths, 'ninety_height': ninety_height, 'ninety_left': ninety_left, 
            'ninety_right': ninety_right, 'half_widths': half_widths, 'half_height': hw_height, 
            'half_left': hw_left, 'half_right': hw_right, 'full_widths': full_widths, 'full_height': fw_height, 
            'full_left': fw_left, 'full_right': fw_right})
        
        all_widths_df = pd.concat([all_widths_df, all_widths_data], ignore_index=True)
    
    else:
        print('No peaks in {} sweep {}'.format(data.file_id, sweep))
    
if len(all_widths_df) > 0:
    data.all_widths_df = all_widths_df.set_index(['sweep #', 'peak #'], inplace=False)
else:
    print('No peaks in {}'.format(data.file_id))

try:
    data.all_widths_df
except:
    print('No peaks in {}'.format(data.file_id))
else:
    print('There are peaks in {}'.format(data.file_id))

for sweep in data.all_widths_df.index.levels[0]:     # this only pulls out sweeps that had peaks   
    trace = data.sweeps.iloc[window_start:, sweep].values

    ninety_left = data.all_widths_df.loc[(sweep), 'ninety_left']
    ten_left = data.all_widths_df.loc[(sweep), 'ten_left']
    ten_to_ninety = (ninety_left - ten_left) / fs
    peak_time = data.all_widths_df.loc[(sweep), 'peaks_index'] / fs
    hw_time = data.all_widths_df.loc[(sweep), 'half_widths'] / fs
    
    peak_numbers = range(len(data.all_widths_df.loc[sweep]))

    tau_list = []
    charge_list = []

    for peak_number in peak_numbers:        
        tau = data.get_tau(trace, sweep, peak_number)
        tau_list.append(tau)
        charge = data.get_charge_transf(trace, sweep, peak_number)
        charge_list.append(charge) 
    
    all_peaks_kinetics_data = pd.DataFrame({'sweep #': sweep,
                                            'peak #': peak_numbers,
                                            'peak time (ms)': peak_time,
                                            '10 to 90% RT (ms)': ten_to_ninety,
                                            'tau': tau_list,
                                            'half-width (ms)': hw_time,
                                            'charge transferred (pA * s)': charge_list})

    all_peaks_kinetics_data = all_peaks_kinetics_data[all_peaks_kinetics_data.tau < 500] # drops peaks with tau values over 500

    if len(all_peaks_kinetics_data) == 0: # moves onto the next sweep if no peaks exist after dropping tau values
        continue
    
    data.plot_first3_events(all_peaks_kinetics_data, sweep, trace)
    delay_to_response = all_peaks_kinetics_data.iloc[0, 2] # gets the time of the first peak

    # calculates the average kinetic values for all peaks in the given sweep 
    all_peaks_kinetics_avg_data = pd.DataFrame(all_peaks_kinetics_data.mean(axis=0)).T
    all_peaks_kinetics_avg_data['sweep #'] = all_peaks_kinetics_avg_data['sweep #'].astype(int) # makes sweep # an int
    all_peaks_kinetics_avg_data.drop(['peak #', 'peak time (ms)'], axis=1, inplace=True) #drop peak # and peak time columns
    all_peaks_kinetics_avg_data.insert(1, 'delay_to_response (ms)', delay_to_response)

    # calculates the average kinetics values for the first three peaks in a sweep   
    first3_kinetics_avg_data = pd.DataFrame(all_peaks_kinetics_data.iloc[:3].mean(axis=0)).T
    first3_kinetics_avg_data['sweep #'] = first3_kinetics_avg_data['sweep #'].astype(int) # makes sweep # an int
    first3_kinetics_avg_data.drop(['peak #', 'peak time (ms)'], axis=1, inplace=True) #drop peak # and peak time columns
    first3_kinetics_avg_data.insert(1, 'delay_to_response (ms)', delay_to_response)

    all_peaks_kinetics_df = pd.concat([all_peaks_kinetics_df, all_peaks_kinetics_data], ignore_index=True)
    all_peaks_kinetics_avg_df = pd.concat([all_peaks_kinetics_avg_df, all_peaks_kinetics_avg_data], ignore_index=True)
    first3_kinetics_avg_df = pd.concat([first3_kinetics_avg_df, first3_kinetics_avg_data], ignore_index=True)
    
# this df contains all kinetics values for all the peaks in all the sweeps
data.all_peaks_kinetics_df = all_peaks_kinetics_df.set_index(['sweep #', 'peak #'], inplace=False)

# this df contains avg kinetic values for all the peaks in all the sweeps
data.all_peaks_kinetics_avg_df = all_peaks_kinetics_avg_df.set_index(['sweep #'], inplace=False)

# this df contains avg kinetic values for the first three peaks in all the sweeps
data.first3_kinetics_avg_df = first3_kinetics_avg_df.set_index(['sweep #'], inplace=False)

# this df calculates the average of averages for the first three peaks in all the sweeps
data.avg_first3_kinetics_avg_df = pd.DataFrame(first3_kinetics_avg_df.mean(axis=0)).T        
data.avg_first3_kinetics_avg_df = pd.concat([data.metadata, data.avg_first3_kinetics_avg_df], axis=1)
data.avg_first3_kinetics_avg_df.drop(['sweep #'], axis=1, inplace=True)
