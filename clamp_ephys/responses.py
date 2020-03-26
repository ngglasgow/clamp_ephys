

def cell_response_checker(data, fs, stim_time, baseline, post_stim):
    baseline = mean_baseline(fs, data, stim_time)
    peaks = epsc_peak(fs=25, data=data, baseline=baseline, stim_time=500, post_stim=150)

    '''
        Pull out EPSC peaks from filtered signals
        Baseline 100 ms preceding blue light
        Peak within 250 ms of blue light
    '''

    # filter signal with butterworth filter at 500 Hz for data
    filt_data = elephant.signal_processing.butter(data.T,
                                                lowpass_freq=500.0,
                                                fs=10000.0)
    filt_data = pd.DataFrame(filt_data).T

    filt_baseline = mean_baseline(fs=fs, data=filt_data, stim_time=500)
    filt_baseline_std = std_baseline(fs=fs, data=filt_data, stim_time=500)
    filt_peaks = epsc_peak(fs=fs, data=filt_data, baseline=filt_baseline, stim_time=500, post_stim=150)


    # make binary choice of whether the response is more than 3x std
    mean_std = filt_baseline_std.mean()
    mean_peaks = filt_peaks.mean()

    response_2x = abs(mean_peaks) > mean_std * 2
    response_3x = abs(mean_peaks) > mean_std * 3

    responses = pd.DataFrame({'Cell name': file,
                            'Mean Peaks (pA)': mean_peaks,
                            'Mean STD (pA)': mean_std,
                            'Response 2x STD': response_2x,
                            'Response 3x STD': response_3x,
                            },
                            index=range(1))

    responses_df = pd.concat([responses_df, responses], ignore_index=True)

    return responses_df
