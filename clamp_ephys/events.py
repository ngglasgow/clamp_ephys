# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Trial:
    def __init__(self, x, y):
        self.x = x
        self.y = y


    def find_events(self, start_threshold, end_threshold):
        self.xs = []
        self.ys = []
        self.durs = []

        starts = np.where(self.y > start_threshold, self.y, 0)
        stops = np.where(self.y > end_threshold, self.y, 0)

        start_indices = np.nonzero(np.diff(np.sign(starts)) == 1)[0] + 1
        stop_indices = np.nonzero(np.diff(np.sign(stops)) == -1)[0]

        self.num_events = len(stop_indices)

        for event in range(self.num_events):
            event_x = self.x[start_indices[event]:stop_indices[event] + 1]
            event_y = self.y[start_indices[event]:stop_indices[event] + 1]
            duration = len(event_x) / 10

            self.xs.append(event_x)
            self.ys.append(event_y)
            self.durs.append(duration)
    
    def events_in_window(self, window_start, window_stop, duration_threshold):
        self.is_event_trial = []
        for event in range(self.num_events):
            event_x = self.xs[event]
            event_y = self.ys[event]
            duration = self.durs[event]

            in_window = np.any(
                event_x[(event_x > window_start) & (event_x < window_stop)]
            )

            if in_window == True:
                if duration >= duration_threshold:
                    is_event = 1
                    self.is_event_trial.append(is_event)
                else:
                    is_event = 0
                    self.is_event_trial.append(is_event)
            else:
                is_event = 0
                self.is_event_trial.append(is_event)

    def any_event_in_window(self):
        if sum(self.is_event_trial) >= 1:
            event_trial = 1

        elif sum(self.is_event_trial) == 0:
            event_trial = 0
    
        return event_trial

    def plot_events(self):
        for event in range(self.num_events):
            plt.plot(self.xs[event], self.ys[event], c='k')
    
    def plot_events_vs_signal(self):
        plt.figure()
        plt.plot(x, y, c='k')
        for event in range(self.num_events):
            plt.plot(self.xs[event], self.ys[event], c='r')

    def return_pandas(self):
        event_xs = pd.DataFrame(self.xs).T
        event_ys = pd.DataFrame(self.ys).T
        events = pd.DataFrame({
            'event index': range(self.num_events),
            'durations': self.durs,
            'in window and duration': self.is_event_trial})

        return events


x = np.linspace(0, 100, 10000)
y = np.sin(x)

start_threshold = 0.5
end_threshold = 0.25

window_start = 1
window_stop = 2
duration_threshold = 10

test_trial = Trial(x, y)
test_trial.find_events(start_threshold, end_threshold)
test_trial.events_in_window(window_start, window_stop, duration_threshold)
event_in_trial = test_trial.any_event_in_window()
test_trial.plot_events()
test_trial.plot_events_vs_signal()
trial_events = test_trial.return_pandas()
