#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 14:52:58 2023

@author: chambeca
@author: demeycdr
"""

"""
A set of functions that help to adjust pick times and pick polarities.

"""
# Import packages
from colours import colours
import os
from numbers import Number
from obspy.core.event import Pick, WaveformStreamID, QuantityError, Event
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.pyplot import Figure, Axes
import numpy as np
from obspy import read, read_events, UTCDateTime, Stream

"""
Define the active keys.
"""
# GLOBALS dictionary for setting what keys correspond to.
KEY_MAPPING = {
    "up": "positive",
    "down": "negative",
    "?": "undecidable",
    "b": "fuckup",
    "n": "next",
    "d": "remove",
    "right": 1,
    "left": -1,
    "pageup": 10,
    "pagedown": -10,
}

# GLOBALS dictionary but key-value pairs are inverted.
INVERSE_KEY_MAPPING = {value: key for key, value in KEY_MAPPING.items()}

# Subdictionary of GLOBALS with only the keys to pick polarities.
POLARITY_KEYS = {key for key, value in KEY_MAPPING.items()
                 if value in ["positive", "negative", "undecidable"]}

# SUbdictionary of GLOBALS with only the keys to adjust pick times.
TIME_KEYS = {key for key, value in KEY_MAPPING.items()
             if isinstance(value, Number)}


# Assign colour dictionary for plotting phase arrivals.
COLORS = {"P": "blue", "S": "orange"}

"""
"""
"""
Overall function that is able to adjust pick times, pick polarities, process and plot waveforms.
The function is build upon a set of subfunctions that can largely be subdived into two sets
based on the way they influence the interactive figure.

The first set only contains the "adjust_pick_polarity" function. 
This function changes the pick information using the assigned keys. 
This function also changes all aspects of the figure which relate to the pick time,
e.g. x- and y-axis limits.

The second set of functions contains the processing functions such as filter, zoom and reset functions. 
These functions do not change the pick information, but change the waveform using the defined buttons and sliders.
Activating an assigned key does not affect these functions.

To ensure a smooth transition between the use of both sets of functions occurs, 
it is necessary for most processing steps to be present in both sets of functions
as they both influence x- and y-axes limits.
This makes the overall function seem somewhat clunky and repetitive.
"""
"""
"""
def pick_polarity(
    pick: Pick,
    stream: Stream,
    GaMMA_picks: list,
    GaMMA_times: list,
    GaMMA_events: list,
    pre_pick: float,
    post_pick: float,
    fig: Figure = None,
    ax: Axes = None,
    axlowcut: Axes = None,
    axhighcut: Axes = None,
    filtax_5: Axes = None,
    filtax_10: Axes = None,
    filtax_15: Axes = None,
    filtax_45: Axes = None,
    reset_filtax: Axes = None,
    axzoom: Axes = None,
    reset_filt_zoomax: Axes = None,
    axamp: Axes = None,
    reset_ampax: Axes = None,
):
    """
    Parameters:
        
        - pick: ObsPy pick object.
        - stream: ObsPy stream object.
            * If the pick's phase hint is P, the stream should contain only a single trace e.g. the HHZ channel.
            * If the pick's phase hint is S, the stream should contain two traces e.g. both the HHE and HHN channels.
        - GaMMA_picks: list of ObsPy formatted picks that are associated by GaMMA.
            * If the pick is this set of GaMMA associated picks this is visually indicated on the plot.
        - GaMMA_times: list of the pick times of the different GaMMA associated picks.
        - GaMMA_events: list of the different events that are associated using the picks in GaMMA_picks.
        - pre_pick: value before the pick at which the trace(s) is/are sliced,
                    minimal value of the x-axis limits.
        - post_pick: value after the pick at which the trace(s) is/are sliced,
                     maximal value of the x-axis limits.
        - fig: figure on which is plotted.
        - ax: axis of the figure on which is plotted.
        - axlowcut: axis of the slider to set the lower frequency when filtering.
        - axhighcut: axis of the slider to set the higher frequency when filtering.
        - filtax_5: axis of the button that sets the values of the lowcut and highcut sliders to 1 and 5 Hz,
                    and filters the stream using these set values.
        - filtax_10: axis of the button that sets the values of the lowcut and highcut sliders to 5 and 10 Hz,
                     and filters the stream using these set values.
        - filtax_15: axis of the button that sets the values of the lowcut and highcut sliders to 10 and 15 Hz,
                     and filters the stream using these set values.
        - filtax_45: axis of the button that sets the values of the lowcut and highcut sliders to 1 and 45 Hz,
                     and filters the stream using these set values. These values are used by Earthquake Tranformer
                     to filter the waveform before detecting earthquakes.
        - reset_filtax: axis of the button that resets the values of the lowcut and highcut sliders to their
                        original values.
        - axzoom: axis of the slider that specifies the length of the shown waveform --> horizontal zoom.
        - reset_filt_zoomax: axis of the button that resets the lowcut, highcut and zoom sliders to their
                             original values.
        - axamp: axis of the slider to specify the amplification --> vertical zoom.
        - reset_ampax: axis of the button that resets the slider to its original value.
    """
    
    # Import packages for plotting interactive figures.
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button

    """
    Define the figure.
    """
    # If no figure was previously initialised, make a figure and axes objects.
    if fig is None:
        
        # Makes a figure consisting of vertically stacked subplots.
        # If the studied pick's phase hint is P, only one axis is defined, 
        # if its phase hint is S two axes are defined.
        fig, ax = plt.subplots(len(stream), 1, figsize=(12, 8))
        
    # If a figure was previously initialised, redefine the axes of the figure and clear all objects in the axes.
    else:
        
        # Redefine the figure and axes objects.
        ax = ax or fig.gca()
        
        # Clear the axes of any plotted objects.
        ax.clear()  

    """
    Assign, slice and plot trace for P picks.
    """
    # If the pick's phase hint is P, the stream object should contain only one channel, e.g. HHZ.
    # All initial plotting steps are followed once.
    if pick.phase_hint == "P": 
        
        """
        Assign and slice trace.
        """
        # Assign the single trace to its own object.
        # This trace will be plotted on the figure.
        trace = stream[0]
        
        # Slice a copy of the full trace between the pre- and post-pick time values.
        # This copy will be used to define the y-axis limits of the figure, when it is first displayed.
        sliced_trace = stream[0].slice(pick.time - pre_pick, pick.time + post_pick).copy()
        
        """
        Print warnings and return empty results if problems occur.
        """
        # If the trace contains no data or the pick time is outside the trace,
        # print a warning and return the following:
            # pick = None
            # fuckup = False
            # _quit = False
        if len(trace.data) == 0:
            
            print(f"No data for {trace.id}")
            
            return None, False, False
        
        if pick.time < trace.stats.starttime or pick.time > trace.stats.endtime:
            
            print("Pick outside data, cannot check")
            
            return None, False, False
        
        """
        Plot trace and set x- and y-axis limits for a rolling window.
        """
        # Define the start time of the trace.
        starttime = trace.stats.starttime
        
        # Define the pick time in relative time, i.e. seconds after the trace's start time.
        pick_time = pick.time - starttime
        
        # Define the time difference between subsequent samples in the trace.
        delta = trace.stats.delta
        
        # Define and plot the waveform as a 2D line element based on its relative time since its start (.times()) and
        # its data. This wavefrom will be displayed first
        seismo_line = ax.add_line(
            Line2D(xdata = trace.times(), ydata = trace.data))
        
        # Set the initial x-axis limits of the figure as the time between the pre- and post-pick time values.
        # Setting up the x-axis limits in this way allows for a rolling window if chaning pick times
        # instead of changing the line that indicates pick time.
        # This is usefull for significant pick time adjustments.
        ax.set_xlim((pick_time - pre_pick, pick_time + post_pick))
        
        # Set the intial y-axis limits of the plot as the minimal and maximal values of the trace in the set x-axis limits.
        # This translates to the minimal and maximal values of the sliced trace.
        ax.set_ylim((sliced_trace.data.min(), sliced_trace.data.max()))
        
    """
    Assign, slice and plot traces for S picks.
    """
    # If the pick's phase hint is S, the stream object should contain two traces, e.g. HHE and HHN.
    # All plotting steps are duplicated, so the plotting is conducted for each trace individually.
    if pick.phase_hint == "S":
        
        """
        Assign and slice traces.
        """
        # Assign both traces to their own objects.
        # These traces will be plotted on the figure
        trace_1 = stream[0]
        trace_2 = stream[1]
        
        # Slice a copy of each full trace between the pre- and post-pick time values.
        # These copies will be used to define the y-axis limits of each subplot, when they are first displayed.
        sliced_trace_1 = stream[0].slice(pick.time - pre_pick, pick.time + post_pick)
        sliced_trace_2 = stream[1].slice(pick.time - pre_pick, pick.time + post_pick)
        
        """
        Print warnings and return empty results if problems occur.
        """
        # If the trace contains no data or the pick time is outside the trace,
        # print a warning and return the following:
            # pick = None
            # fuckup = False
            # _quit = False
        if len(trace_1.data) == 0:
            
            print(f"No data for {trace_1.id}")
            
            return None, False, False
        
        if len(trace_2.data) == 0:
            
            print(f"No data for {trace_2.id}")
            
            return None, False, False
        
        if pick.time < trace_1.stats.starttime or pick.time > trace_1.stats.endtime:
            
            print("Pick outside data, cannot check")
            
            return None, False, False

        if pick.time < trace_2.stats.starttime or pick.time > trace_2.stats.endtime:
            
            print("Pick outside data, cannot check")
            
            return None, False, False
        
        """
        Plot traces and set x- and y-axis limits for a rolling window.
        """
        # Define the start time of one of the traces.
        # The code assumes both traces have the same start and end times.
        starttime = trace_1.stats.starttime
        
        # Define the pick time in relative time, i.e. seconds after the trace's start time.
        pick_time = pick.time - starttime
        
        # Define the time difference between subsequent samples in the trace.
        # The code assumes that both traces have the same delta value.
        delta = trace_1.stats.delta

        # Define and plot the intial waveforms of each trace as a 2D line element based on 
        # its relative time since its start (.times()) and its data.
        seismo_line_1 = ax[0].add_line(
            Line2D(xdata = trace_1.times(), ydata = trace_1.data))
        seismo_line_2 = ax[1].add_line(
            Line2D(xdata = trace_2.times(), ydata = trace_2.data))
        
        # Set the intial x-axis limits of each subplot as the time between the pre- and post-pick time values.
        # Setting up the x-axes limits in this way allows for a rolling window if chaning pick times
        # instead of changing the line that indicates pick time.
        # This is usefull for significant pick time adjustments.
        ax[0].set_xlim((pick_time - pre_pick, pick_time + post_pick))
        ax[1].set_xlim((pick_time - pre_pick, pick_time + post_pick))
        
        # Set the initial y-axis limits of each subplot as the minimal and maximal values of the trace
        # in the set x-axis limits. This translates to the minimal and maximal values of the sliced trace.
        ax[0].set_ylim((sliced_trace_1.data.min(), sliced_trace_1.data.max()))
        ax[1].set_ylim((sliced_trace_2.data.min(), sliced_trace_2.data.max()))

    """
    Warning if stream is empty.
    """
    
    # If the stream contains no traces, print a warning and return the following:
        # pick = None
        # fuckup = False
        # _quit = False
    if len(stream) == 0:
        
        print("No traces in stream")
        
        return None, False, False

    """
    Assign pick colour.
    """
    
    # If the pick time of this pick is in the set of provided GaMMA associated pick times,
    # the pick is visually indicated by a red colour.
    if pick.time in GaMMA_times:
        
        c = "red"
    
    # If the pick time is not in the set of provided GaMMA associated pick times,
    # the pick is plotted using its normal colour from the colour dictionary.
    else:
        c = COLORS[pick.phase_hint]
    
    """
    Plot pick line.
    """
    ### Plot the vertical lines indicating pick time.
    
    # If the pick's phase hint is P, indicate the pick time by plotting a vertical line on the figure
    # using the correct plot colour.
    if pick.phase_hint == "P":
        
        # Define and plot a vertical line on the plot. 
        # X-coordinates: relative pick times in seconds.
        # Y-coordinates: extent of the y-axis limits.
        line = ax.add_line(
            Line2D(xdata = [pick_time, pick_time],
                   ydata = list(ax.get_ylim()), color = c))
    
    # If the pick's phase hint is S, indicate the pick time by plotting a vertical line on each subplot
    # using the correct plot colour.
    if pick.phase_hint == "S":
        
        # Define and plot a vertical line on each subplot. 
        # X-coordinates: relative pick times in seconds.
        # Y-coordinates: extent of the y-axis limits.
        line_1 = ax[0].add_line(
            Line2D(xdata = [pick_time, pick_time],
                   ydata = list(ax[0].get_ylim()), color = c))
        line_2 = ax[1].add_line(
            Line2D(xdata = [pick_time, pick_time],
                   ydata = list(ax[1].get_ylim()), color = c))
    
    """
    Assign trace text for title.
    """
    ### Define part of the figure title text indicating the trace ID('s) plotted.
    
    # If the pick time is in the list of provided GaMMA associated pick times,
    # make a list of all associated picks for each GaMMA associated event and
    # define the text that indicates the pick's trace in the figure title using this information.
    if pick.time in GaMMA_times:
        
        # For each event in the list of GaMMA associated events obtain the picks associated to this event,
        # and  make a list of the pick times of these associated picks.
        for event in GaMMA_events:
            
            # Obtain the picks associated to the event.
            event_picks = event.picks
            
            # For each pick in the list of associated picks,
            # see if the pick time of the associated pick equals the pick time of the single pick provided as input.
            for p in event_picks:
                
                # Check if the inputted pick and associated pick have the same pick time
                if pick.time == p.time:
                    
                    # Get the index of the event in the list of associated GaMMA events.
                    # Base value = 1 --> add a 1 to the index
                    event_index = GaMMA_events.index(event) + 1
                    
                    # If the pick's phase hint is P, define the trace text using the trace ID and event index.
                    if pick.phase_hint == "P":
                        
                        trace_text = trace.id + " GaMMA associated to event " + str(event_index)
                    
                    # If the pick's phase hint is S, define the trace text using the trace ID's of both traces
                    # and the event index.
                    else:
                        
                        trace_text = trace_1.id + " and " + trace_2.id + " GaMMA associated to event " + str(event_index)
                        
    # If the pick time is not in the list of the provided GaMMA associated picks,
    # define the text that indicates the pick's trace in the figure title in its simple form.
    else:
        
        # If the pick's phase hint is P, use the single trace ID.
        if pick.phase_hint == "P":
            
            trace_text = trace.id
        
        # If the pick's phase hint is S, combine both trace ID's.
        else:
            trace_text = trace_1.id + " " + trace_2.id
    
    """
    Set subplot title, x-axis label and polarity text.
    """
    ### Add figure title, axes labels and polarity text.
    
    # If the pick's phase hint is P, set the title and x-axis label of the single subplot.
    if pick.phase_hint == "P":
        
        # Set the x-axis label position to the top of the plot.
        ax.xaxis.set_label_position('top')
        
        # Define the x-axis label as the time since the overall trace start time.
        ax.set_xlabel(f"Seconds from {trace.stats.starttime}")
        
        # Set the figure title.
        ax.set_title(f"{trace_text}: pick {INVERSE_KEY_MAPPING['positive']} "
                     f"or {INVERSE_KEY_MAPPING['negative']}, or move "
                     f"{INVERSE_KEY_MAPPING[1]} or {INVERSE_KEY_MAPPING[-1]}. "
                     f"{INVERSE_KEY_MAPPING['next']} to ignore")

        # Add a text line to the figure indicating the picked polarity.
        pol_text = ax.text(0.1, 0.9, f"Polarity: {pick.polarity}", transform = ax.transAxes)
        
    # If the pick's phase hint is S, set the title and x-axis label on the correct subplot.
    elif pick.phase_hint == 'S':
        
        # Set the x-axis label position to the top of the plot.
        ax[0].xaxis.set_label_position('top')
        
        # Define the x-axis label as the time since the overall trace start.
        # Assumes both traces have the same start time
        ax[0].set_xlabel(f"Seconds from {stream[0].stats.starttime}")
        
        ax[0].set_title(f"{trace_text}: pick {INVERSE_KEY_MAPPING['positive']} "
                     f"or {INVERSE_KEY_MAPPING['negative']}, or move "
                     f"{INVERSE_KEY_MAPPING[1]} or {INVERSE_KEY_MAPPING[-1]}. "
                     f"{INVERSE_KEY_MAPPING['next']} to ignore")
        
        # Add a text line to the figure indicating the picked polarity.
        pol_text = ax[0].text(0.1, 0.9, f"Polarity: {pick.polarity}", transform = ax[0].transAxes)

    """
    Add frequency, zoom and amplification buttons and sliders.
    These will be used as the input for the processing functions.
    """
    ### Add the lowcut frequency slider.
    
    # If no lowcut slider was defined, adjust the bottom of the figure so that it can include the slider.
    if axlowcut is None:
        fig.subplots_adjust(bottom=0.25)

    # If an axis element for the lowcut slider was previously defined, clear it.
    if axlowcut:
        
        axlowcut.clear()
    
    # If no axis element for the lowcut slider was previously defined, define it.
    else:
        
        axlowcut = fig.add_axes([0.25, 0.2, 0.65, 0.03])
    
    # Based on the pick's phase hint define the lowcut slider using the sampling rate of the correct trace.
    # Max value: half the trace sampling rate
    # Initial value: 0 Hz
    if pick.phase_hint == "P":
        
        low_slider = Slider(
            ax = axlowcut,
            label = 'Lowcut [Hz]',
            valmin = 0.0,
            valstep = 0.5,
            valmax = trace.stats.sampling_rate / 2,
            valinit = 0.0,
        )
        
    else:
        
        low_slider = Slider(
            ax = axlowcut,
            label = 'Lowcut [Hz]',
            valmin = 0.0,
            valstep = 0.5,
            valmax = trace_1.stats.sampling_rate / 2,
            valinit = 0.0,
        )
        
    ### Add the highcut frequency slider.
    
    # If an axis element for the lowcut slider was previously defined, clear it.
    if axhighcut:
        
        axhighcut.clear()
        
    # If no axis element for the lowcut slider was previously defined, define it.
    else:
        
        axhighcut = fig.add_axes([0.25, 0.15, 0.65, 0.03])
    
    # Based on the pick's phase hint define the lowcut slider using the sampling rate of the correct trace.
    # Max value: half the trace sampling rate
    # Initial value: max value
    if pick.phase_hint == "P":
        
        high_slider = Slider(
            ax = axhighcut,
            label = 'Highcut [Hz]',
            valmin = 0.0,
            valstep = 0.5,
            valmax = trace.stats.sampling_rate / 2,
            valinit = trace.stats.sampling_rate / 2,
        )
        
    else:
        
        high_slider = Slider(
            ax = axhighcut,
            label = 'Highcut [Hz]',
            valmin = 0.0,
            valstep = 0.5,
            valmax = trace_1.stats.sampling_rate / 2,
            valinit = trace_1.stats.sampling_rate / 2,
        )
        
        
    ### Add the button to filter between 1 - 5 Hz.
    
    # If an axis element for the 1 - 5 Hz button was previously defined, clear it.
    if filtax_5:
        
        filtax_5.clear()
    
    # If no axis element for the 1 - 5 Hz button was previously defined, define it.
    else:
        
        filtax_5 = fig.add_axes([0.05, 0.025, 0.09, 0.04])

    # Define the 1 - 5 Hz filter button.
    filt_5_button = Button(filtax_5, "Filter 1-5 Hz", hovercolor='0.975')
    
    ### Add the button to filter between 5 - 10 Hz.
    
    # If an axis element for the 5 - 10 Hz button was previously defined, clear it.
    if filtax_10:
        
        filtax_10.clear()
    
    # If no axis element for the 5 - 10 Hz button was previously defined, define it.
    else:
        
        filtax_10 = fig.add_axes([0.2, 0.025, 0.09, 0.04])

    # Define the 5 - 10 Hz filter button.
    filt_10_button = Button(filtax_10, "Filter 5-10 Hz", hovercolor='0.975')
    
    ### Add the button to filter between 10 - 15 Hz.
    
    # If an axis element for the 10 - 15 Hz button was previously defined, clear it.
    if filtax_15:
        
        filtax_15.clear()
    
    # If no axis element for the 10 - 15 Hz button was previously defined, define it.
    else:
        
        filtax_15 = fig.add_axes([0.35, 0.025, 0.09, 0.04])

    # Define the 10 - 15 Hz filter button.
    filt_15_button = Button(filtax_15, "Filter 10-15 Hz", hovercolor='0.975')

    ### Add the button to filter between 1 - 45 Hz.
    
    # If an axis element for the 1 - 45 Hz button was previously defined, clear it.
    if filtax_45:
        
        filtax_45.clear()
    
    # If no axis element for the 1 - 45 Hz button was previously defined, define it.
    else:
        
        filtax_45 = fig.add_axes([0.5, 0.025, 0.09, 0.04])

    # Define the 1 - 45 Hz filter button.
    filt_45_button = Button(filtax_45, "Filter 1-45 Hz", hovercolor='0.975')
    
    ### Add the button to reset the filter.
    
    # If an axis element for the filter reset button was previously defined, clear it.
    if reset_filtax:
        
        reset_filtax.clear()
    
    # If no axis element for the filter reset button was previously defined, define it.
    else:
        
        reset_filtax = fig.add_axes([0.8, 0.025, 0.09, 0.04])
        
    # Define the reset filter button.
    reset_filt_button = Button(reset_filtax, "Reset filter", hovercolor='0.975')
    
    ### Add the horizontal zoom slider.
    
    # If an axis element for the zoom slider was previously defined, clear it.
    if axzoom:
        
        axzoom.clear()
        
    # If no axis element for the zoom slider was previously defined, define it.
    else:
    
        axzoom = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        
    # Define the zoom slider.
    # The value of the zoom slider defines the total length of the shown trace.
    # Initial value: twice the pre- or post-pick time, assuming equal pre- and post-pick values.
    # Min value: 0.5 seconds, cannot be 0 as otherwise no waveform is plotted.
    # Zooming occurs as follows:
        # The value of the zoom slider is set as the new shown trace length.
        # A copy of the trace is sliced with half of this new trace length before and after the pick time.
        # Decreasing the zoom value thus shortens the shown trace length and zooms in.
    zoom_slider = Slider(
        ax = axzoom,
        label = 'Zoom (s)',
        valmin = 0.5,
        valstep = 0.5,
        valmax = pre_pick*2,
        valinit = pre_pick*2,
    )
    
    ### Add the reset zoom and filter button.
    
    # If an axis element for the reset zoom and filter button was previously defined, clear it.
    if reset_filt_zoomax:
        
        reset_filt_zoomax.clear()
        
    # If no axis element for the reset zoom and filter button was previously defined, define it.
    else:
    
        reset_filt_zoomax = fig.add_axes([0.65, 0.025, 0.09, 0.04])

    # Define the reset zoom and filter button.
    reset_filt_zoom_button = Button(reset_filt_zoomax, "Reset filter and zoom", hovercolor='0.975')

    ### Add the amplification slider.
    
    # If an axis element for the amplification slider was previously defined, clear it.
    if axamp:
        
        axamp.clear()
    
    # If no axis element for the amplification slider was previously defined, define it.
    else:
        
        axamp = fig.add_axes([0.05, 0.25, 0.03, 0.65])
        
    # Define the amplification slider.
    # Slider that defines the amplification magnification factor.
    # Min value: 1, no change in amplification.
    # Amplification occurs by multiplying the trace data by the value of the amplification slider.
    amplification_slider = Slider(
        ax = axamp,
        label = "amplification",
        valmin = 1.0,
        valstep = 0.5,
        valmax = 10.0,
        valinit = 1.0,
        orientation = "vertical"
        )
    
    ### Add the amplification reset button.
    
    # If an axis element for the amplification reset button was previously defined, clear it.
    if reset_ampax:
        
        reset_ampax.clear()
        
    # If no axis element for the amplification reset button was previously defined, define it.
    else:
        
        reset_ampax = fig.add_axes([0.05, 0.15, 0.09, 0.04])

    # Define the amplification reset button.
    amp_reset_button = Button(reset_ampax, "Reset amplification", hovercolor='0.975')
    
    """
    """
    """
    Set of functions that change the plot in-place when polarities are picked, pick times are adjusted or
    the waveforms are filtered, zoomed in and/or amplified.
    """
    """
    """
    # Define the trackers that are assigned during the pick time adjusting and polarity picking.
    fuckup, _quit, polarity_picked, time_adjusted, delete = False, False, False, 0, False
    

    ### Function that is used to pick polarities and adjust pick times.
    ### This function also processes the waveform. The processed waveform is then used to
    ### obtain the correct x- and y-axes limits for the rolling window when pick times are changed.
    def adjust_pick_polarity(event):
        
        # Import parameters from outside the function.
        nonlocal pick, fuckup, _quit, polarity_picked, time_adjusted, delta, pol_text, delete
        
        # Import sliders
        nonlocal low_slider, high_slider, zoom_slider, amplification_slider
        
        # If the pick's phase hint is P, import the trace and line objects for the single trace.
        if pick.phase_hint == "P":
            
            nonlocal trace, line
        
        # If the pick's phase hint is P, import the trace and line objects for both traces.
        else:
            
            nonlocal trace_1, trace_2, line_1, line_2
        
        """
        Get the values used for waveform processing.
        """
        # Obtian the values of the lowcut and highcut frequency, zoom and amplification sliders.
        lowcut = low_slider.val
        highcut = high_slider.val
        zoom = zoom_slider.val
        amp = amplification_slider.val
            
        # If the lowcut value is higher than the highcut value, give a warning.
        if lowcut > highcut:
            
            print("Low greater than high, ignoring")
            
            return
        
        """
        Pick polarities.
        """
        # If one of the polarity keys is activated, change the pick polarity and plot it on the figure.
        if event.key in POLARITY_KEYS:
            
            # Change the pick polarity.
            pick.polarity = KEY_MAPPING.get(event.key)
            
            # Define if the polarity is picked or not.
            polarity_picked = True
            
            # Change the polarity text to the updated polarity.
            pol_text.set_text(f"Polarity: {pick.polarity}")
            
            # Plot this updated text on the figure.
            pol_text.axes.draw_artist(pol_text)
            pol_text.figure.canvas.draw()

            return
        
        """
        Adjust pick times.
        """
        # If one the pick time keys is activated, change the pick time and plot it on the figure.
        if event.key in TIME_KEYS:
            
            # If the pick's phase hint is P, work on the single trace.
            if pick.phase_hint == "P":

                """
                Get the adjusted pick time.
                """
                # Update the pick time.
                # The new pick time is the old pick time + the sampling time difference * the amount of times the time key is pressed.
                pick.time += (KEY_MAPPING[event.key] * delta)
                
                # Define the time difference with which the pick time is adjusted.
                time_adjusted += (KEY_MAPPING[event.key] * delta)
                
                # Calculate the relative pick time from the trace start in seconds.
                pick_time = pick.time - starttime
                
                """
                Waveform processing.
                This waveform processing is conducted to change the x- and y-axes limits for the rolling window.
                These steps use the updated pick time to set the x-axis limits after, e.g. zooming.
                These updated x-axis limits are used slice the trace after it is filtered and/or amplified.
                This sliced trace is used to set the y-axis limits.
                Changing the pick time continuously thus results in smooth transitions in x- and y-axes limits.
                """
                
                ### Filter the raw trace.
                # If the lowcut and highcut sliders are on their original values, reset the trace to its unflitered state.
                if lowcut == low_slider.valinit and highcut == high_slider.valinit:
                    
                    # Set the filtered trace object to its unfiltered state.
                    filtered = trace
                    
                # If the lowcut slider value is its original value, but if the highcut slider value is not its original value,
                # detrend, taper and filter a copy of the raw trace using a lowpass filter.
                elif lowcut == low_slider.valinit:
        
                    filtered = trace.copy().detrend().taper(0.1).filter("lowpass", freq = highcut)
                    
                # If the highcut slider value is its original value, but if the lowcut slider value is not its original value,
                # detrend, taper and filter a copy of the raw trace using a highpass filter.
                elif highcut == high_slider.valinit:
                    
                    filtered = trace.copy().detrend().taper(0.1).filter("highpass", freq = lowcut)
                
                # If both the lowcut and highcut slider values or not their original values,
                # detrend, taper and filter a copy of the raw trace using a bandpass filter.
                else:
                    
                    filtered = trace.copy().detrend().taper(0.1).filter("bandpass", freqmin = lowcut, freqmax = highcut)
                
                ### Slice the filtered trace depending on the zoom value. This sliced trace is used to update x-axis limits.
                # Set the x-axis limits for the rolling window. 
                # Changing the pick time will also shift the x-axis limits.
                if zoom == zoom_slider.valinit:
                    
                    # Slice the trace again using the updated pick time. 
                    # If no zoom is applied, use the pre- post pick time values.
                    sliced_trace = filtered.slice(pick.time - pre_pick, pick.time + post_pick).copy()
                    
                    # If no zoom is applied set the rolling window to the initial window length defined by the pre- and post-pick values.
                    x_lim = (pick_time - pre_pick, pick_time + post_pick)
                
                # If a zoom is applied set the rolling window length to half the zoom value before and after the pick time.
                else:
                    
                    # Slice the trace again using the updated pick time. 
                    # If zoomed, use the length of the zoomed window.
                    sliced_trace = filtered.slice(pick.time - zoom / 2, pick.time + zoom / 2).copy()
                    
                    # If zoomed set the rolling window to the length defined by the zoom slider value.
                    x_lim = (pick_time - zoom / 2, pick_time + zoom / 2)
                
                ### Amplify the filtered data.
                # Multiply the data of the filtered trace with the value of the amplification slider.
                amped_data = filtered.data * amp
                
                # Define a amplified trace object wich is a copy of the filtered trace object, containing the amplified data
                amped = filtered.copy()
                amped.data = amped_data
                
                # Set the y-axis depending on the amplification.
                if amp == amplification_slider.valinit:
                    
                    # If no amplification is applied, set the y-axis limits as the minimal and maximal values of the sliced data.
                    y_lim = (sliced_trace.data.min(), sliced_trace.data.max())
                    
                else:
                    
                    # When amplification occurs, use the same y-axis limits width, i.e. the difference between the
                    # minimal and maximal y-axis limits is the same as for the unamplified data.
                    # Otherwise the influence of amplification is not observed.
                    # However, the data value at the updated pick time may lay outside of the initial y-axis limits.
                    # Therefore, determine the data value of the amplified waveform at the adjusted pick time.
                    pick_data = np.where(amped.times() == pick_time)
                    
                    # Set the updated y-axis limits with the same width as the initial y-axis limits,
                    # i.e. the second term is the difference of the two min and max limits, around the data value of the amplified trace
                    # on the adjusted pick time.
                    y_lim = (amped.data[pick_data[0]] - (sliced_trace.data.max() - sliced_trace.data.min()) / 2,  
                             amped.data[pick_data[0]] + (sliced_trace.data.max() - sliced_trace.data.min()) / 2)
                
                # Set the x- and y-axes limits.
                ax.set_xlim(x_lim)
                ax.set_ylim(y_lim)
                
                # Update the line element so that is has the new pick time and extends to the updated y-axis limits.
                line.set_data([pick_time, pick_time], list(line.axes.get_ylim()))
                
                # Draw the updated pick time line on the figure.
                line.axes.draw_artist(line)
                line.figure.canvas.draw()
            
            # If the pick's phase hint is S, work on both traces.
            else:
                
                """
                Get the adjusted pick time.
                """
                # Update the pick time.
                # The new pick time is the old pick time + the sampling time difference * the amount of times the time key is pressed.
                pick.time += (KEY_MAPPING[event.key] * delta)
                
                # Define the time difference with which the pick time is adjusted.
                time_adjusted += (KEY_MAPPING[event.key] * delta)
                
                # Calculate the relative pick time from the trace start in seconds.
                pick_time = pick.time - starttime
                

                """
                Waveform processing.
                This waveform processing is conducted to change the x- and y-axes limits for the rolling window.
                These steps use the updated pick time to set the x-axis limits after, e.g. zooming.
                These updated x-axis limits are used slice the trace after it is filtered and/or amplified.
                This sliced trace is used to set the y-axis limits.
                Changing the pick time continuously thus results in smooth transitions in x- and y-axes limits.
                """
                
                ### Filter the raw trace.
                # If the lowcut and highcut sliders are on their original values, reset the trace to its unflitered state.
                if lowcut == low_slider.valinit and highcut == high_slider.valinit:
                    
                    # Set the filtered trace object to its unfiltered state.
                    filtered_1 = trace_1
                    filtered_2 = trace_2
                    
                # If the lowcut slider value is its original value, but if the highcut slider value is not its original value,
                # detrend, taper and filter a copy of the raw trace using a lowpass filter.
                elif lowcut == low_slider.valinit:
        
                    filtered_1 = trace_1.copy().detrend().taper(0.1).filter("lowpass", freq = highcut)
                    filtered_2 = trace_2.copy().detrend().taper(0.1).filter("lowpass", frea = highcut)
                    
                # If the highcut slider value is its original value, but if the lowcut slider value is not its original value,
                # detrend, taper and filter a copy of the raw trace using a highpass filter.
                elif highcut == high_slider.valinit:
                    
                    filtered_1 = trace_1.copy().detrend().taper(0.1).filter("highpass", freq = lowcut)
                    filtered_2 = trace_2.copy().detrend().taper(0.1).filter("highpass", freq = lowcut)
                
                # If both the lowcut and highcut slider values or not their original values,
                # detrend, taper and filter a copy of the raw trace using a bandpass filter.
                else:
                    
                    filtered_1 = trace_1.copy().detrend().taper(0.1).filter("bandpass", freqmin = lowcut, freqmax = highcut)
                    filtered_2 = trace_2.copy().detrend().taper(0.1).filter("bandpass", freqmin = lowcut, freqmax = highcut)
                
                ### Slice the filtered trace depending on the zoom value. This sliced trace is used to update x-axis limits.
                # Set the x-axis limits for the rolling window. 
                # Changing the pick time will also shift the x-axis limits.
                if zoom == zoom_slider.valinit:
                    
                    # Slice the trace again using the updated pick time. 
                    # If no zoom is applied, use the pre- post pick time values.
                    sliced_trace_1 = filtered_1.slice(pick.time - pre_pick, pick.time + post_pick).copy()
                    sliced_trace_2 = filtered_2.slice(pick.time - pre_pick, pick.time + post_pick).copy()
                    
                    # If no zoom is applied set the rolling window to the initial window length defined by the pre- and post-pick values.
                    x_lim = (pick_time - pre_pick, pick_time + post_pick)
                
                # If a zoom is applied set the rolling window length to half the zoom value before and after the pick time.
                else:
                    
                    # Slice the trace again using the updated pick time. 
                    # If zoomed, use the length of the zoomed window.
                    sliced_trace_1 = filtered_1.slice(pick.time - zoom / 2, pick.time + zoom / 2).copy()
                    sliced_trace_2 = filtered_2.slice(pick.time - zoom /2, pick.time + zoom / 2).copy()
                    
                    # If zoomed set the rolling window to the length defined by the zoom slider value.
                    x_lim = (pick_time - zoom / 2, pick_time + zoom / 2)
                
                ### Amplify the filtered data.
                # Multiply the data of the filtered trace with the value of the amplification slider.
                amped_data_1 = filtered_1.data * amp
                amped_data_2 = filtered_2.data * amp
                
                # Define a amplified trace object wich is a copy of the filtered trace object, containing the amplified data
                amped_1 = filtered_1.copy()
                amped_1.data = amped_data_1
                amped_2 = filtered_2.copy()
                amped_2.data = amped_data_2
                
                # Set the y-axis depending on the amplification.
                if amp == amplification_slider.valinit:
                    
                    # If no amplification is applied, set the y-axis limits as the minimal and maximal values of the sliced data.
                    y_lim_1 = (sliced_trace_1.data.min(), sliced_trace_1.data.max())
                    y_lim_2 = (sliced_trace_2.data.min(), sliced_trace_2.data.max())
                    
                else:
                    
                    # When amplification occurs, use the same y-axis limits width, i.e. the difference between the
                    # minimal and maximal y-axis limits is the same as for the unamplified data.
                    # Otherwise the influence of amplification is not observed.
                    # However, the data value at the updated pick time may lay outside of the initial y-axis limits.
                    # Therefore, determine the data value of the amplified waveform at the adjusted pick time.
                    pick_data_1 = np.where(amped_1.times() == pick_time)
                    pick_data_2 = np.where(amped_2.times() == pick_time)
                    
                    # Set the updated y-axis limits with the same width as the initial y-axis limits,
                    # i.e. the second term is the difference of the two min and max limits, around the data value of the amplified trace
                    # on the adjusted pick time.
                    y_lim_1 = (amped_1.data[pick_data_1[0]] - (sliced_trace_1.data.max() - sliced_trace_1.data.min()) / 2,  
                             amped_1.data[pick_data_1[0]] + (sliced_trace_1.data.max() - sliced_trace_1.data.min()) / 2)
                    y_lim_2 = (amped_2.data[pick_data_2[0]] - (sliced_trace_2.data.max() - sliced_trace_2.data.min()) / 2,
                               amped_2.data[pick_data_2[0]] + (sliced_trace_2.data.max() - sliced_trace_2.data.min()) / 2)
                
                # Set the x- and y-axes limits.
                ax[0].set_xlim(x_lim)
                ax[0].set_ylim(y_lim_1)
                ax[1].set_xlim(x_lim)
                ax[1].set_ylim(y_lim_2)
                
                # Update the line elements so that they have the new pick time and extends to the updated y-axis limits.
                line_1.set_data([pick_time, pick_time], list(line_1.axes.get_ylim()))
                line_2.set_data([pick_time, pick_time], list(line_2.axes.get_ylim()))
                
                
                # Draw the updated pick time lines on the figure.
                line_1.axes.draw_artist(line_1)
                line_1.figure.canvas.draw()
                line_2.axes.draw_artist(line_2)
                line_2.figure.canvas.draw()
        
        # If the key is a single number, set the lowcut slider value to the value of the pressed key.
        # Seems redundant with the dedicated filter buttons.
        elif event.key in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
            
            low_slider.set_val(float(event.key))
        
        # If the pressed key is "n", go to the next pick.
        elif event.key == INVERSE_KEY_MAPPING["next"]:
            
            # Go to the next pick.
            fig.canvas.stop_event_loop()
            return
        
        # If the pressed key is "d", remove the pick and move to the next one.
        elif event.key == INVERSE_KEY_MAPPING["remove"]:
            
            # Set the remove statement to True.
            delete = True
            
            # Go to the next pick.
            fig.canvas.stop_event_loop()
            
            return
        
        # If the pressed key is "b", go back.
        elif event.key == INVERSE_KEY_MAPPING["fuckup"]:
            
            # Go back to the previous one...
            print("Fucked-up, going back")
            
            # Set the go back statement is True.
            fuckup = True
            
            # Go to the previous pick
            fig.canvas.stop_event_loop()
            
            return
        
        # If the pressed key is "q", close the program.
        elif event.key == "q":
            
            print("Quitting")
            
            # Set the quit statement to True.
            _quit = True
            
            # Exit the event.
            fig.canvas.stop_event_loop()
            return
        else:
            print(f"{event.key} is not bound")

    ### Function that actively changes the plotted waveform, as it apllies a filter, zoom and amplification.
    def zoom_and_filter_data(event):
        
        # Import parameters.
        nonlocal high_slider, low_slider, zoom_slider, pick, time_adjusted, amplification_slider, ax
        
        # If the pick's phase hint is P, import the waveform, pick line and trace objects for the single trace.
        if pick.phase_hint == "P":
            
            nonlocal seismo_line, line, trace
            
        # If the pick's phase hint is S, import the waveform, pick line and trace objects for both traces.
        else:
            nonlocal seismo_line_1, seismo_line_2, line_1, line_2, trace_1, trace_2
        
        # Obtian the values of the lowcut and highcut frequency, zoom and amplification sliders.
        lowcut = low_slider.val
        highcut = high_slider.val
        zoom = zoom_slider.val
        amp = amplification_slider.val
            
        # If the lowcut value is higher than the highcut value, give a warning.
        if lowcut > highcut:
            
            print("Low greater than high, ignoring")
            
            return
         
        # If the pick's phase hint is P, work on the single trace.
        # Define a filtered trace object that is subsequently, sliced, amplified and zoomed in.
        if pick.phase_hint == "P":
            
            # If the lowcut and highcut sliders are on their original values, reset the trace to its unflitered state.
            if lowcut == low_slider.valinit and highcut == high_slider.valinit:
                
                # Set the filtered trace object to its unfiltered state.
                filtered = trace
                
            # If the lowcut slider value is its original value, but if the highcut slider value is not its original value,
            # detrend, taper and filter a copy of the raw trace using a lowpass filter.
            elif lowcut == low_slider.valinit:
    
                filtered = trace.copy().detrend().taper(0.1).filter("lowpass", freq = highcut)
                
            # If the highcut slider value is its original value, but if the lowcut slider value is not its original value,
            # detrend, taper and filter a copy of the raw trace using a highpass filter.
            elif highcut == high_slider.valinit:
                
                filtered = trace.copy().detrend().taper(0.1).filter("highpass", freq = lowcut)
            
            # If both the lowcut and highcut slider values or not their original values,
            # detrend, taper and filter a copy of the raw trace using a bandpass filter.
            else:
                
                filtered = trace.copy().detrend().taper(0.1).filter("bandpass", freqmin = lowcut, freqmax = highcut)
                
            # Define a sliced trace object wich is a slice copy of the filtered trace object, with a window length equal to the zoom length.
            sliced = filtered.slice(pick.time - zoom / 2, pick.time + zoom / 2).copy()
            
            # Obtain the start time of the trace
            starttime = trace.stats.starttime
            
            # Define the relative pick time from the trace start time.
            pick_time = pick.time - starttime
            
            # Set the x-axis limit.
            x_lim = ((pick_time - zoom / 2, pick_time + zoom / 2))
            
            # Multiply the data of the filtered trace with the value of the amplification slider.
            amped_data = filtered.data * amp
            
            # Define a amplified trace object wich is a copy of the filtered trace object, containing the amplified data
            amped = filtered.copy()
            amped.data = amped_data
            
            # Set the y-axis depending on the amplification.
            if amp == amplification_slider.valinit:
                
                # If no amplification is applied, set the y-axis limits as the minimal and maximal values of the sliced data.
                y_lim = (sliced.data.min(), sliced.data.max())
                
            else:
                
                # When amplification occurs, use the same y-axis limits width, i.e. the difference between the
                # minimal and maximal y-axis limits is the same as for the unamplified data.
                # Otherwise the influence of amplification is not observed.
                # However, the data value at the updated pick time may lay outside of the initial y-axis limits.
                # Therefore, determine the data value of the amplified waveform at the adjusted pick time.
                pick_data = np.where(amped.times() == pick_time)
                
                # Set the updated y-axis limits with the same width as the initial y-axis limits,
                # i.e. the second term is the difference of the two min and max limits, around the data value of the amplified trace
                # on the adjusted pick time.
                y_lim = (amped.data[pick_data[0]] - (sliced.data.max() - sliced.data.min()) / 2,  
                         amped.data[pick_data[0]] + (sliced.data.max() - sliced.data.min()) / 2)
                
                
            # Set the x- and y-axes limits.
            ax.set_xlim(x_lim)
            ax.set_ylim(y_lim)

            # Update the line element so that is has the new pick time and extends to the updated y-axis limits.
            line.set_data([pick_time, pick_time], list(line.axes.get_ylim()))
            
            # Draw the updated pick time line on the figure.
            line.axes.draw_artist(line)
            line.figure.canvas.draw()
            
            # Update the waveform element with the new, filtered and amplified data. This is displayed when processing occurs.
            seismo_line.set_data(amped.times(), amped.data)
            
            # Draw the updated waveform element
            ax = seismo_line.axes
            ax.draw_artist(seismo_line)
            seismo_line.figure.canvas.draw()
    
        # If the pick's phase hint is S, work on the both traces.
        # Define a filtered traces object that is subsequently, sliced, amplified and zoomed in.
        else:
            
            # If the lowcut and highcut sliders are on their original values, reset the trace to its unflitered state.
            if lowcut == low_slider.valinit and highcut == high_slider.valinit:
                
                # Set the filtered trace object to its unfiltered state.
                filtered_1 = trace_1
                filtered_2 = trace_2
                
            # If the lowcut slider value is its original value, but if the highcut slider value is not its original value,
            # detrend, taper and filter a copy of the raw trace using a lowpass filter.
            elif lowcut == low_slider.valinit:
    
                filtered_1 = trace_1.copy().detrend().taper(0.1).filter("lowpass", freq = highcut)
                filtered_2 = trace_2.copy().detrend().taper(0.1).filter("lowpass", freq = highcut)
                
            # If the highcut slider value is its original value, but if the lowcut slider value is not its original value,
            # detrend, taper and filter a copy of the raw trace using a highpass filter.
            elif highcut == high_slider.valinit:
                
                filtered_1 = trace_1.copy().detrend().taper(0.1).filter("highpass", freq = lowcut)
                filtered_2 = trace_2.copy().detrend().taper(0.1).filter("highpass", freq = lowcut)
            
            # If both the lowcut and highcut slider values or not their original values,
            # detrend, taper and filter a copy of the raw trace using a bandpass filter.
            else:
                
                filtered_1 = trace_1.copy().detrend().taper(0.1).filter("bandpass", freqmin = lowcut, freqmax = highcut)
                filtered_2 = trace_2.copy().detrend().taper(0.1).filter("bandpass", freqmin = lowcut, freqmax = highcut)
                
            # Define a sliced trace object wich is a slice copy of the filtered trace object, with a window length equal to the zoom length.
            sliced_1 = filtered_1.slice(pick.time - zoom / 2, pick.time + zoom / 2).copy()
            sliced_2 = filtered_2.slice(pick.time - zoom / 2, pick.time + zoom / 2).copy()
            
            # Obtain the start time of the trace
            starttime = trace_1.stats.starttime
            
            # Define the relative pick time from the trace start time.
            pick_time = pick.time - starttime
            
            # Set the x-axis limit.
            x_lim = ((pick_time - zoom / 2, pick_time + zoom / 2))
            
            # Multiply the data of the filtered trace with the value of the amplification slider.
            amped_data_1 = filtered_1.data * amp
            amped_data_2 = filtered_2.data * amp
            
            # Define a amplified trace object wich is a copy of the filtered trace object, containing the amplified data
            amped_1 = filtered_1.copy()
            amped_1.data = amped_data_1
            amped_2 = filtered_2.copy()
            amped_2.data = amped_data_2
            
            # Set the y-axis depending on the amplification.
            if amp == amplification_slider.valinit:
                
                # If no amplification is applied, set the y-axis limits as the minimal and maximal values of the sliced data.
                y_lim_1 = (sliced_1.data.min(), sliced_1.data.max())
                y_lim_2 = (sliced_2.data.min(), sliced_2.data.max())
                
            else:
                
                # When amplification occurs, use the same y-axis limits width, i.e. the difference between the
                # minimal and maximal y-axis limits is the same as for the unamplified data.
                # Otherwise the influence of amplification is not observed.
                # However, the data value at the updated pick time may lay outside of the initial y-axis limits.
                # Therefore, determine the data value of the amplified waveform at the adjusted pick time.
                pick_data_1 = np.where(amped_1.times() == pick_time)
                pick_data_2 = np.where(amped_2.times() == pick_time)
                
                # Set the updated y-axis limits with the same width as the initial y-axis limits,
                # i.e. the second term is the difference of the two min and max limits, around the data value of the amplified trace
                # on the adjusted pick time.
                y_lim_1 = (amped_1.data[pick_data_1[0]] - (sliced_1.data.max() - sliced_1.data.min()) / 2,  
                           amped_1.data[pick_data_1[0]] + (sliced_1.data.max() - sliced_1.data.min()) / 2)
                y_lim_2 = (amped_2.data[pick_data_2[0]] - (sliced_2.data.max() - sliced_2.data.min()) / 2,  
                           amped_2.data[pick_data_2[0]] + (sliced_2.data.max() - sliced_2.data.min()) / 2)
                
                
            # Set the x- and y-axes limits.
            ax[0].set_xlim(x_lim)
            ax[0].set_ylim(y_lim_1)
            ax[1].set_xlim(x_lim)
            ax[1].set_ylim(y_lim_2)

            # Update the line element so that is has the new pick time and extends to the updated y-axis limits.
            line_1.set_data([pick_time, pick_time], list(line_1.axes.get_ylim()))
            line_2.set_data([pick_time, pick_time], list(line_2.axes.get_ylim()))
            
            # Draw the updated pick time line on the figure.
            line_1.axes.draw_artist(line_1)
            line_1.figure.canvas.draw()
            line_2.axes.draw_artist(line_2)
            line_2.figure.canvas.draw()
            
            # Update the waveform element with the new, filtered and amplified data. This is displayed when processing occurs.
            seismo_line_1.set_data(amped_1.times(), amped_1.data)
            seismo_line_2.set_data(amped_2.times(), amped_2.data)
            
            # Draw the updated waveform element
            ax[0] = seismo_line_1.axes
            ax[0].draw_artist(seismo_line_1)
            seismo_line_1.figure.canvas.draw()
            ax[1] = seismo_line_2.axes
            ax[1].draw_artist(seismo_line_2)
            seismo_line_2.figure.canvas.draw()
    
    ### Fucntion that filters the data between 1 - 5 Hz after the corresponding button input.
    def filter_data_5(event):
        
        # Import nonlocal parameters.
        nonlocal high_slider, low_slider

        # Define the filtering values
        lowcut = 1
        highcut = 5
        
        # Reset the highcut and lowcut sliders if they are not on their original positions.
        high_slider.reset()
        low_slider.reset()
        
        # Set the highcut and lowcut sliders to the filtering values.
        high_slider.set_val(highcut)
        low_slider.set_val(lowcut)
        
        # Filter the data using the corresponding filter.
        zoom_and_filter_data(event)
        
    ### Fucntion that filters the data between 5 - 10 Hz after the corresponding button input.  
    def filter_data_10(event):
        
        # Import nonlocal parameters.
        nonlocal high_slider, low_slider

        # Define the filtering values
        lowcut = 5
        highcut = 10
        
        # Reset the highcut and lowcut sliders if they are not on their original positions.
        high_slider.reset()
        low_slider.reset()
        
        # Set the highcut and lowcut sliders to the filtering values.
        high_slider.set_val(highcut)
        low_slider.set_val(lowcut)
        
        # Filter the data using the corresponding filter.
        zoom_and_filter_data(event)
            
    ### Fucntion that filters the data between 10 - 15 Hz after the corresponding button input.
    def filter_data_15(event):
        
        # Import nonlocal parameters.
        nonlocal high_slider, low_slider

        # Define the filtering values
        lowcut = 10
        highcut = 15
        
        # Reset the highcut and lowcut sliders if they are not on their original positions.
        high_slider.reset()
        low_slider.reset()
        
        # Set the highcut and lowcut sliders to the filtering values.
        high_slider.set_val(highcut)
        low_slider.set_val(lowcut)
        
        # Filter the data using the corresponding filter.
        zoom_and_filter_data(event)
            
    ### Fucntion that filters the data between 1 - 45 Hz after the corresponding button input.        
    def filter_data_45(event):
        
        # Import nonlocal parameters.
        nonlocal high_slider, low_slider

        # Define the filtering values
        lowcut = 1
        highcut = 45
        
        # Reset the highcut and lowcut sliders if they are not on their original positions.
        high_slider.reset()
        low_slider.reset()
        
        # Set the highcut and lowcut sliders to the filtering values.
        high_slider.set_val(highcut)
        low_slider.set_val(lowcut)
        
        # Filter the data using the corresponding filter.
        zoom_and_filter_data(event)
    
    ### Fucntion that resets the zoom and filter sliders to its original position and processes the waveform as such.
    def unzoom(event):

        # Import nonlocal parameters
        nonlocal high_slider, low_slider, zoom_slider
        
        print("Reseting")
        
        # Reset the sliders to their original positions.
        high_slider.reset()
        low_slider.reset()
        zoom_slider.reset()
        
        # Process the waveform with these new slider values.
        zoom_and_filter_data(event)
        
    ### Fucntion that resets the amplification slider and processes the waveform as such.
    def deamp(event):
        
        # Import nonlocal parameters.
        nonlocal amplification_slider
        
        print("Reseting")
        
        # Reset the amplification slider to its original position.
        amplification_slider.reset()
        
        # Process the waveform with the new value of the amplification slider.
        zoom_and_filter_data(event)
        
    ### Function that resets the frequency sliders to their original positions and filters the data as such.
    def unfilter(event):
        
        # Import nonlocal parameters
        nonlocal high_slider, low_slider
        
        print("Reseting")
        
        # Reset the frequency sliders to their original positions.
        high_slider.reset()
        low_slider.reset()
        
        # Filter the data using the new values of the filter sliders.
        zoom_and_filter_data(event)
        
    # Define which functions are activated when the different buttons are pressed or the sliders changed.
    filt_5_button.on_clicked(filter_data_5)
    filt_10_button.on_clicked(filter_data_10)
    filt_15_button.on_clicked(filter_data_15)
    filt_45_button.on_clicked(filter_data_45)
    low_slider.on_changed(zoom_and_filter_data)
    high_slider.on_changed(zoom_and_filter_data)
    zoom_slider.on_changed(zoom_and_filter_data)
    amplification_slider.on_changed(zoom_and_filter_data)
    reset_filt_zoom_button.on_clicked(unzoom)
    reset_filt_button.on_clicked(unfilter)
    amp_reset_button.on_clicked(deamp)
    
    # Attach responder
    cid = fig.canvas.mpl_connect('key_press_event', adjust_pick_polarity)

    manager = plt.get_current_fig_manager()
    manager.full_screen_toggle()
    fig.show()
    fig.canvas.draw()  # Redraw
    fig.canvas.start_event_loop()
    # Picking happens here
    fig.canvas.mpl_disconnect(cid)
    # If we want to delete the pick
    if delete:
        plt.close()
        return None, fuckup, _quit
    if not fuckup and not _quit:
        if pick.phase_hint == "P":
            if polarity_picked:
                print(f"{trace.id} picked as {pick.polarity}")
                print(f"{trace.id} adjusted by {time_adjusted}")
        else:
            if polarity_picked:
                print(f"{trace_1.id} and {trace_2.id} picked as {pick.polarity}")
            if time_adjusted:
                print(f"{trace_1.id} and {trace_2.id} adjusted by {time_adjusted}")
    plt.close()
    return pick, fuckup, _quit

### Function that reads in a pick CSV file and converts them to ObsPy formatted picks.
def csv_2_ObsPy(filename):
    
    # Define an empty list where the pick information will be stored.
    picks = []
    
    # Open the CSV file using Pandas.
    csv_file = pd.read_csv("/Users/home/demeycdr/Downloads/picks_CSV//" + filename)
    
    # For each row in the csv_file, append the information to the pick dataframe in the format: id - timestamp - prob - type.
    for i in range(len(csv_file["Network"])):

        pick = Pick(time = csv_file["Time"][i],
            waveform_id = WaveformStreamID(station_code = csv_file["Station"][i],
                channel_code = None, network_code = csv_file["Network"][i],
                location_code = csv_file["Location"][i]),
            phase_hint = csv_file["Phase"][i],
            evaluation_mode = "automatic",
            creation_info = csv_file["Creation_info"][i],
            time_errors = QuantityError(confidence_level = csv_file["Confidence"][i],
                lower_uncertainty = csv_file["Lower_Uncertainty"][i],
                upper_uncertainty = csv_file["Upper_Uncertainty"][i],
                uncertainty = csv_file["Uncertainty"][i]))

        picks.append(pick)

    # Return the formatted picks.
    return picks

wav_directory = "/Users/home/demeycdr/Downloads/wav//"

# Define the path to the GaMMA associated event QuakeML files.
QuakeML_directory = "/Users/home/demeycdr/Downloads/EQT_QuakeML//" 

# Convert the files in this directory to a list.
QuakeML_list = os.listdir(QuakeML_directory)

wav_list = os.listdir(wav_directory)

# For some reason the above command does not read in the files in alfphabetical order, i.e. increasing time, and this ticks me off.
# Make an empty list that will hold the waveform start times as UTCDateTime objects.
wav_time_list = []

# For each file append the UTCDateTime object to this list.
for filename in wav_list:
    
    wav_time_list.append(UTCDateTime(int(filename[:4]), int(filename[5:7]), int(filename[8:10]), 
                                     int(filename[11:13]), int(filename[14:16]), int(filename[17:])))

# Sort both the UTCDateTime and waveform filename list according to increasing start time.
wav_list = [x for _,x in sorted(zip(wav_time_list, wav_list))] 

for filename in wav_list:
    
        if filename + ".xml" not in os.listdir("/Users/home/demeycdr/Downloads/Adjusted_QuakeML//"):
        
            print(colours.GREEN + "Working on earthquake " + filename + colours.ENDC)
        
            st = read(wav_directory + filename)
            
            EQT_picks = csv_2_ObsPy(filename + "_avg.csv")
            
            stations = []
            for tr in st:
                if tr.stats.station not in stations:
                    stations.append(tr.stats.station)
                    
            earliest_picks = []
            corrected_stations = []
            for station in stations:
                picks = []
                for pick in EQT_picks:
                    if pick.waveform_id.station_code == station:
                        picks.append(pick.time)
                if len(picks) > 0:
                    earliest_picks.append(min(picks))
                    corrected_stations.append(station)
                
            corrected_stations = [x for _,x in sorted(zip(earliest_picks, corrected_stations))]       
            
            print(corrected_stations)
            
            picks_GaMMA = []
            events_GaMMA = []
            times_GaMMA = []
            
            # For each file in this QuakeML list look if the name corresponds to the waveform filename.
            for QuakeML_file in QuakeML_list:
                
                # If the QuakeML filename coincides with the waveform filename.
                if QuakeML_file[:19] == filename and QuakeML_file[-7:] == "avg.xml":
                    
                    print(QuakeML_file)
                    
                    # Read the QuakeML file.
                    event = read_events(QuakeML_directory + QuakeML_file)[0]
                    
                    # Append the event to the list of associated events
                    events_GaMMA.append(event)
                    
                    # Get the picks from the event.
                    picks = event.picks
                    
                    # For each pick change the channel code so that it holds information on which event it was associated to.
                    for pick in picks:
                        
                        # Append the pick to the list of GaMMA associated picks.
                        picks_GaMMA.append(pick)
                        times_GaMMA.append(pick.time)
                    
            adjusted_picks = []
                    
            for station in corrected_stations:
                print(station)
                
                for pick in EQT_picks:
                    
                    if pick.waveform_id.station_code == station:
                        
                        sub_st = st.select(station = station)
                        
                        if pick.phase_hint == "P":
                            
                            pick_st = sub_st.select(channel = "HHZ")
                            
            
                            pick, fuckup, _quit = pick_polarity(pick = pick, stream = pick_st,
                                                                GaMMA_picks = picks_GaMMA,
                                                                GaMMA_times = times_GaMMA,
                                                                GaMMA_events = events_GaMMA,
                                                                pre_pick = 5, post_pick = 5)
        
                            if pick != None:
                                adjusted_picks.append(pick)
                                
                        elif pick.phase_hint == "S":
                            pick_st_1 = sub_st.select(channel = "HHN")
                            pick_st_2 = sub_st.select(channel = "HHE")
                            pick_st_1_alt = sub_st.select(channel = "HH1")
                            pick_st_2_alt = sub_st.select(channel = "HH2")
                            if len(pick_st_1) != 0:
                                pick_st = pick_st_1 + pick_st_2
                            else:
                                pick_st = pick_st_1_alt + pick_st_2_alt
                            
            
                            pick, fuckup, _quit = pick_polarity(pick = pick, stream = pick_st,
                                                                GaMMA_picks = picks_GaMMA,
                                                                GaMMA_times = times_GaMMA,
                                                                GaMMA_events = events_GaMMA,
                                                                pre_pick = 5, post_pick = 5)
                            
                            if pick != None:
                                adjusted_picks.append(pick)
            
            adjusted_event = Event(picks = adjusted_picks)
            
            adjusted_event.write("/Users/home/demeycdr/Downloads//Adjusted_QuakeML//" + filename + ".xml", format = "QUAKEML")
            
            print(colours.CYAN + " Adjusted picks" + colours.ENDC)
            for pick in adjusted_picks:
                print(pick.waveform_id.station_code, pick.phase_hint, pick.time)
            
            print(colours.GREEN + "Next eartquake is " + wav_list[wav_list.index(filename) + 1] + colours.ENDC)
            inp = input("Do you want to continue?")
            
            if inp == "Yes" or inp == "yes" or inp == "y":
                continue
            elif inp == "No" or inp == "no" or inp == "n":
                break
            