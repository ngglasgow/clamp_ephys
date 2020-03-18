# -*- coding: utf-8 -*-
"""
Created 5 March 2019
IPSC_peak_x.y.z.py

"""
#%%
#from __main__ import *
# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import neo
import elephant
from neo.io import IgorIO
from neo.core import SpikeTrain, AnalogSignal
import quantities as pq
from elephant.sta import spike_triggered_average
import scipy
from neuron import h
import os

#%%
from neuron import h
# %% User Inputs
file_name = 'JH190131_c1_drug_1.ibw'
depolarized_sweeps = []

# metadata
file_name_split = file_name.split('_')
cell_id = file_name_split[0]+'_'+file_name_split[1]
cell_num = cell_id[-1:]

# pull out date code and cell numbers
date = '20'+cell_id[2:4]+'-'+cell_id[4:6]+'-'+cell_id[6:8]

# define condition e.g. control, drug1, drug2
if 'drug' in file_name:
    condition = 'TTX+4-AP'
else:
    condition = 'control'


# %%
'''
This cell sorts through the cell_id directory and sorts files only noise vm
files into two lits; 1 for drugs and 1 for contorl. Then it opens the correct
file type for the given file_name, and concatenates vm traces if necessary.
Open igor binary file into neo, read neo.block, convert to np.array, convert
to pd.dataframe, concat with empty data frame noise_vmdf and repeat if concat.
'''
os.chdir('../ETC/Data/')

currdf = pd.DataFrame()
temp_curr_raw = IgorIO(filename=file_name)
temp_curr_neo = temp_curr_raw.read_block()
temp_currarray = temp_curr_neo.segments[0].analogsignals[0]
temp_currdf = pd.DataFrame(temp_currarray.as_array())
currdf = pd.concat([currdf, temp_currdf], axis=1)
temp_noise_vm_raw = None
temp_noise_vm_block = None
temp_noise_vmarray = None
os.chdir('../..')



for item in depolarized_sweeps:
    noise_vmdf[item] = np.nan
# %%
'''
#Plot all the traces on one plot
plt.figure()
plt.plot(currdf)

# %%
# Plot individual traces
for i in range(len(currdf.columns)):
    plt.figure()
    plt.plot(currdf.iloc[:,i])
'''

#%%
# Add a time column to dataframe for current
t = pd.DataFrame(np.arange(0,6000, 0.1))
currdf = pd.concat([t, currdf], axis=1)
currdf.columns = np.arange(0,51, 1)

#%% Exclude sweeps
exclude = depolarized_sweeps

#%% Extract window for baseline - 260 ms to 490 ms (after test pulse)
baseline_window = currdf.iloc[2600:4900]
mean_baselines = baseline_window.mean()
mean_baselines = pd.DataFrame(mean_baselines)

#%% Extract window for IPSC peak - 500 ms to 2000 ms (light onset to 2000 ms)
IPSC_window = currdf.iloc[5000:10000]
peak_IPSC_list = []
for i in range(len(IPSC_window.columns)):
#    ninetieth = IPSC_window.iloc[:,i].quantile(0.9, interpolation='nearest')
    min_IPSC = IPSC_window.iloc[:,i].min()
    peak_IPSC_list.append(min_IPSC)
peak_IPSC_df = pd.DataFrame(peak_IPSC_list)


    
# subtract baseline
baselineadj_peak_IPSC = peak_IPSC_df - mean_baselines
avg_peak_IPSC = baselineadj_peak_IPSC.mean()
avg_peak_IPSC = avg_peak_IPSC[0]


#%% Calculating series resistance
sealtest_window = currdf.iloc[600:2400]
'''
sealtest_window = currdf.iloc[400:600]

peak_captrans_list = []
for i in range(len(sealtest_window.columns)):
    peak_captrans = sealtest_window.iloc[:,i].max()
    peak_captrans_list.append(peak_captrans)
peak_captrans_df = pd.DataFrame(peak_captrans_list)
delta_current_df = mean_baselines - peak_captrans_df
'''

# using avg current in step of seal test, not transient peak
sealtest_avg_list = []
for i in range(len(sealtest_window.columns)):
    sealtest_avg = sealtest_window.iloc[:,i].mean()
    sealtest_avg_list.append(sealtest_avg)
sealtest_avg_df = pd.DataFrame(sealtest_avg_list)

delta_current_df = mean_baselines - sealtest_avg_df
Rs_df = delta_current_df.rdiv(5)


# %%

sig_st = pd.DataFrame()
for i in range(len(noise_vmdf.columns)):
    tmpvm = h.Vector(noise_vmdf.iloc[:, i].values)
    tmpst = h.Vector()
    tmpst.spikebin(tmpvm, -10)
    tmpstdf = pd.DataFrame(tmpst.to_python())
    sig_st = pd.concat([sig_st, tmpstdf], axis=1)

# Add a time column to all data frames for Vm and for binary spike times
t = pd.DataFrame(np.arange(0, 4000, 0.1))
tst = pd.concat([t, sig_st], axis=1)

neospiketrain_list = []
spiketime_list = []
half_neospiketrain_list = []
half_spiketime_list = []


for i in range(1,len(noise_vmdf.columns)+1):
    fulltst = tst.iloc[:, [0, i]]  # single df with time and one binary spike
    halftst = tst.iloc[7000:32001, [0, i]]
    fulltst = fulltst.iloc[:, 0][fulltst.iloc[:, 1] > 0]
    spiketime_list.append(fulltst)
    fulltst = SpikeTrain(times=fulltst.values*pq.ms, t_start=100, t_stop=3500)
    neospiketrain_list.append(fulltst)
    halftst = halftst.iloc[:, 0][halftst.iloc[:, 1] > 0]
    half_spiketime_list.append(halftst)
    halftst = SpikeTrain(times=halftst.values*pq.ms, t_start=700, t_stop=3200)
    half_neospiketrain_list.append(halftst)

# %% Exclude sweeps
#Exclude sweeps
exclude = depolarized_sweeps
del_neospiketrain_list = list(neospiketrain_list)
del_neospiketrain_list= [spiketrain for i, spiketrain in
                         enumerate(del_neospiketrain_list) if
                         i not in depolarized_sweeps]

del_half_neospiketrain_list = list(half_neospiketrain_list)
del_half_neospiketrain_list = [spiketrain for i, spiketrain in
                         enumerate(del_half_neospiketrain_list) if
                         i not in depolarized_sweeps]

del_half_spiketime_list = list(half_spiketime_list)
del_half_spiketime_list = [spiketimes for i, spiketimes in
                         enumerate(del_half_spiketime_list) if
                         i not in depolarized_sweeps]

nan_half_neospiketrain_list = list(half_neospiketrain_list)
nan_half_spiketime_list = list(half_spiketime_list)

for item in exclude:
    nan_half_neospiketrain_list[item]=np.nan
    nan_half_spiketime_list[item]=np.nan

number_of_trials = len(del_half_neospiketrain_list)

# %% To make a graph raster plot like in elephant tutorial.
plt.figure()
for i, spiketrain in enumerate(del_neospiketrain_list):
    t = spiketrain.rescale(pq.ms)
    plt.plot(t, (i+1)*np.ones_like(t), 'k.', markersize=2)
plt.axis('tight')
plt.xlim(100, 3300)
plt.xticks(np.arange(0, 3250, 500))
plt.xlabel('Time (ms)', fontsize=12)
plt.ylim(0, number_of_trials+1)
plt.ylabel('Trial Number', fontsize=12)
plt.gca().tick_params(axis='both', which='major', labelsize=12)
plt.show()

os.chdir('Analysis/Noise_spike_rasters')
#plt.savefig('Noise_raster_'+file_name+'.png', bbox_inches='tight')
#plt.close()
os.chdir('../..')
print 'Spike raster for ' + file_name + 'saved'

#%%
#doing FI curve
#finding mean firing rate for entire sweep
from elephant.statistics import mean_firing_rate
firing_rate_list = []
for spiketrain in neospiketrain_list:
    temp_firing_rate = mean_firing_rate(spiketrain.rescale(pq.s),t_start=0.5,t_stop=1.5)
    firing_rate_list.append(np.asarray(temp_firing_rate))

# %% conversion of discrete spike times to binary counts
from elephant.conversion import BinnedSpikeTrain
bst_list = BinnedSpikeTrain(del_half_neospiketrain_list,
                            binsize=1.0*pq.ms,
                            t_start=700.0*pq.ms,
                            t_stop=3200.0*pq.ms)

bst_arr = bst_list.to_array()           # export binned spike times to an array
bst_df = pd.DataFrame(bst_arr).T        # turn into a df and transpose (.T)
bst_sum = bst_df.apply(np.sum, axis=1)   # sum by row across columns

# plt.figure()
#plt.plot(bst_sum)

"""
Making PSTH for the whole sweep, without first 500 ms omitted
"""
bst_list_graph = BinnedSpikeTrain(del_neospiketrain_list,
                            binsize=1.0*pq.ms,
                            t_start=200.0*pq.ms,
                            t_stop=3200.0*pq.ms)

bst_arr_graph = bst_list_graph.to_array()           # export binned spike times to an array
bst_df_graph = pd.DataFrame(bst_arr_graph).T        # turn into a df and transpose (.T)
bst_sum_graph = bst_df_graph.apply(np.sum, axis=1)   # sum by row across columns

# make gaussian kernel with 2 ms SD
gauKern = elephant.statistics.make_kernel('GAU', 2*pq.ms, 1*pq.ms)
plt.close()
# convolve kernel with summed spike times
# need to calculate average FR and plot that instead of summed spike counts
averageFR = (bst_sum_graph/(number_of_trials*1))*1000
psth = scipy.convolve(gauKern[0], averageFR)
plt.figure()
plt.title('PSTH')
plt.xlabel('Time (ms)')
plt.xticks(np.arange(0, 3250, 500))
plt.ylabel('Firing Rate Hz')
plt.plot(psth)
#plt.axis('tight')
plt.show()

os.chdir('Analysis/PSTH')  # /PSTH')
plt.savefig('Noise_PSTH_'+file_name+'.png', bbox_inches='tight')
plt.close()
os.chdir('../..')
print 'PSTH for ' + file_name + 'saved'

#%% CV-ISI

"""
The below calculates CV-ISI from only the latter 2.5 s current step, output = half_mean_CV
Method = drop values from neospiketrain_list that occur before 1000 ms, then calculating ISI and CV from that
"""
from elephant.statistics import isi, cv

half_isi_list = [isi(spiketrain) for spiketrain in half_spiketime_list]
half_cv_list = [cv(item) for item in half_isi_list]

half_isi_list = [isis for i, isis in
                         enumerate(half_isi_list) if i not in depolarized_sweeps]

#excluding depolarized sweeps - can't use del_half_spiketime_list bc
#need to sort by current injection later

for item in exclude:
    half_cv_list[item] = np.nan

half_mean_CV = np.nanmean(half_cv_list)

#%% Plot ISI histogram for second 1.5 s half of sweep

# to change list of ISIs (a list of lists) to just a plain list
#half_isi_list = [val for sublist in half_isi_list for val in sublist]
#half_isi_list = [float(i) for i in half_isi_list]
#the np.concatenate method below does the same thing and is simpler

"""
Plot ISI for latter 1.5 s half of sweep
"""

binsrange = np.arange(0, 500, 10)
half_isi_list = np.concatenate(half_isi_list)
plt.figure()
#trying to sum bar heights to 1
weights = np.ones_like(half_isi_list)/float(len(half_isi_list))
plt.hist(half_isi_list, bins = binsrange, weights=weights)

plt.xlabel('ISI (ms)', fontsize=16)
plt.xlim(0, 1500)
plt.ylabel('Probability', fontsize=16)
plt.ylim(0, 0.8)
#fig.tight_layout()
plt.show()
os.chdir('Analysis/ISI')
plt.savefig('Noise_ISI_histogram_'+file_name+'.png', bbox_inches='tight')
plt.close()
os.chdir('..\..')
# %%
# Spike train correlations
# -boxcar(9) is for width 9 ms. In this case width means points wide, and since
#  using 1 ms bins, points is ms. Width 8 ms only gives a width of 7 ms.
# -binnedSpikes is a binary read of spikes in 1 ms bins, I think technically
#  this should not be binned, but still binary at 10 kHz. In that case, the box
#  would be 90 wide instead.
# -df.corr() gives an NxN Pearson correlation matrix. I take the mean of it for
#  a single value.

# take out correlation = 1 AND first 500 ms
box_kernel = scipy.signal.boxcar(9)  # if using neuron, not 9 but 81
spiketime_boxdf = pd.DataFrame()
for i in range(len(bst_df.columns)):
    tmpspiketime_box = scipy.signal.convolve(box_kernel, bst_df.iloc[:, i])
    tmpspiketime_boxdf = pd.DataFrame(tmpspiketime_box)
    spiketime_boxdf = pd.concat([spiketime_boxdf, tmpspiketime_boxdf], axis=1)
st_corr = spiketime_boxdf.corr()
st_corr = st_corr[st_corr < 1]
mean_st_corr = np.mean(st_corr.mean().values)
boxcorr = mean_st_corr

# save the 1 ms binned spiketime dataframe convolved with an 8 ms box to file
spiketime_box_filename = file_name + '_spiketime_box.csv'
os.chdir('Analysis/spiketime')
spiketime_boxdf.to_csv(spiketime_box_filename, float_format='%3f', index=False)
os.chdir('../..')


# %%
"""
Need to exclude the first 500 ms of noise
"""
#open igor binary file, read block, and convert to pd.df (noise = 30 4000 ms-long sweeps, 3000 ms-long noise )

"""
get all the noise_outwaves from Igor and put it in a directory accessible by all cell
files, change the f= location below
"""
os.chdir('Data/' + cell_id)
for file in os.listdir(os.getcwd()):
    if 'noise' in file:
        if cell_id in file:
            noise_file = file
f = neo.io.IgorIO(filename=noise_file)
blk = f.read_block()
'''
# length for new noise current
inj_vmarr = blk.segments[0].analogsignals[0]
inj_vmdf = pd.DataFrame(inj_vmarr.as_array())
additional_inj0s = pd.DataFrame(np.zeros(7000))  # making noise file match length of voltage file
additional_injtime = pd.DataFrame(np.arange(33000, 40000, 1))
inj_vmdf = pd.concat([inj_vmdf, additional_inj0s], axis=0)
inj_vmdf.index = np.arange(0, 40000)
os.chdir('../..')
'''

#'''
# length for old noise currents
inj_vmarr = blk.segments[0].analogsignals[0]
inj_vmdf = pd.DataFrame(inj_vmarr.as_array())
additional_inj0s = pd.DataFrame(np.zeros(5000))  # making noise file match length of voltage file
additional_injtime = pd.DataFrame(np.arange(35000, 40000, 1))
inj_vmdf = pd.concat([inj_vmdf, additional_inj0s], axis=0)
inj_vmdf.index = np.arange(0, 40000)
os.chdir('../..')

wns = inj_vmdf.iloc[7000:32000, :]
wns2 = AnalogSignal(signal=wns.values*pq.nA, sampling_rate=10*pq.kHz,
                    t_start=700.0*pq.ms)

sta_ind = pd.DataFrame()                # starting time column to build on
for i in range(len(del_half_neospiketrain_list)):                     # adding each columb of sta to df
    sta_i = spike_triggered_average(wns2, del_half_neospiketrain_list[i],
                                    (-50*pq.ms, 10.1*pq.ms))
    sta_i_df = pd.DataFrame(sta_i.as_array().flatten())
    sta_ind = pd.concat([sta_ind, sta_i_df], axis=1)
sta_avg = sta_ind.apply(np.mean, axis=1)  # average each row across columns
t_sta = np.arange(-50, 10.1, 0.1)     # time from -50 to 10 ms inclusive by 0.1
#sta_avg.index = t_sta
plt.figure()
plt.plot(t_sta,sta_avg)                  # easier to call as separate x,y
plt.xlabel('Time Before Spike (ms)')
plt.ylabel('Current (pA)')
plt.title('STA')
plt.show()
os.chdir('Analysis/STA')
plt.savefig('STA_'+file_name+'.png', bbox_inches='tight')
plt.close()
os.chdir('../..')

print 'STA for ' + file_name + 'saved'

# %% Fano Factor Calculation
# -The goal is to take a the total spike count over a sliding window for each
#  spike train and calculating the fano factor for it var/mean.
#  Then I take the mean of all fano factors for all spike trains.
#  Since it's binary spikes, if you just sum all the 1's in the window, that's
#  the spike count.
# -I'm using binnedSpikes at 1 ms bins again so the window width here in
#  tmpser.rolling(100) is 100 points wide, in our case, 100 ms wide. S
tfanoser = pd.Series()
for i in (range(len(del_half_neospiketrain_list))):  # 30
    if np.sum(bst_df.iloc[:, i]) == 0:
        continue
    else:
        tmpser = bst_df.iloc[:, i]
        window_spikes = tmpser.rolling(100).sum()
        tmpff = window_spikes.var()/window_spikes.mean()
        tmpffs = pd.Series(tmpff)
        tfanoser = pd.concat([tfanoser, tmpffs])
tfanomean = np.mean(tfanoser)


# %% Create a csv file for summarizing average FI curve stats and append to csv
# create metadata df first, then concat so things are somewhat ordered

noise_metadata_df = pd.DataFrame({'Cell ID': cell_id,
                                  'Cell Number': cell_num,
                                  'Date': date,
                                  'Condition': condition,
                                  'Noise Amplitude (pA)': noise_amp,
                                  'STA Start': '-50 ms',
                                  'STA End': '+10 ms'},
                                 index=range(1))

noise_stats_df = pd.DataFrame({'Spike Train Correlation': boxcorr,
                               'CV-ISI': half_mean_CV,
                               'Fano Factor': tfanomean}, index=range(1))

df_to_append = pd.concat([noise_metadata_df, noise_stats_df], axis=1)
#df_to_append = pd.concat([df_to_append, sta.T], axis=1)

# save noise stats
noise_stats_file_name = 'noise_stats_summary.csv'
os.chdir('Analysis')
dir_list = os.listdir(os.getcwd())
if noise_stats_file_name in dir_list:
    old_file = pd.read_csv(noise_stats_file_name)
    new_file = pd.concat([old_file, df_to_append])
    new_file.to_csv(noise_stats_file_name, float_format='%8.4f', index=False)
    old_file = None
    new_file = None
    df_to_append = None
else:
    df_to_append.to_csv(noise_stats_file_name, float_format='%8.4f', index=False)
    df_to_append = None

sta_to_append = pd.DataFrame({file_name: sta_avg})
sta_file_name = 'noise_sta_summary.csv'
if sta_file_name in dir_list:
    old_file = pd.read_csv(sta_file_name)
    new_file = pd.concat([old_file, sta_to_append], axis=1)
    new_file.to_csv(sta_file_name, float_format='%8.4f', index=False)

else:
    sta_to_append.to_csv(sta_file_name, float_format='%8.4f', index=False)

os.chdir('..')

print 'summary csvs for ' + file_name + ' saved'
