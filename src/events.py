# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import elephant
from neo.io import IgorIO
import os
from collections import OrderedDict
import math
from scipy.signal import *

x = np.linspace(0, 100, 10000)
y = np.sin(x)

start_threshold = 0.5
end_threshold = 0.25

starts = np.where(y > start_threshold, y, 0)
stops = np.where(y > end_threshold, y, 0)

start_indices = np.nonzero(np.diff(np.sign(starts)) == 1)[0] + 1
stop_indices = np.nonzero(np.diff(np.sign(stops)) == -1)[0]

xs, ys, durations = [], [], []

for event in range(len(stop_indices)):
    event_x = x[start_indices[event]:stop_indices[event] + 1]
    event_y = y[start_indices[event]:stop_indices[event] + 1]
    duration = len(event_x) / 10
    xs.append(event_x)
    ys.append(event_y)
    durations.append(duration)

# events = pd.DataFrame({'xs': xs, 'ys': ys, 'durations': durations})

plt.figure()
for event_x, event_y in zip(xs, ys):
    plt.plot(event_x, event_y, c='k')

window_start = 1
window_stop = 2
duration_threshold = 10

for event_x, event_y, duration in zip(xs, ys, durs):
    in_window = np.any(event_x[(event_x > window_start) & (event_x < window_stop)])
    
    if in_window == True:
        if duration >= duration_threshold:
            event = True
            print(event)
        else:
            event = False
            print(event)
    else:
        event = False
        print(event)

