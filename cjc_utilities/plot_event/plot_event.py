"""
Functions for plotting events.

Calum Chamberlain
"""

import numpy as np


def plot_event_from_client(event, client, length=60, size=(10.5, 10.5),
                           all_channels=False, filt=None, ignore_rotated=True,
                           return_stream=False, fig=None):
    """
    Plot the waveforms for an event with pick and calculated arrival times.

    :type event: `obspy.core.event.Event`
    :param event: Event to plot
    :param client: An obspy client with `get_waveforms_bulk` method
    :type length: float
    :param length: Length to plot, from origin time.
    :type all_channels: bool
    :param all_channels: Whether to download all channels from that sensor.
    """
    from obspy import Stream
    from obspy.clients.fdsn.client import FDSNNoDataException, FDSNException

    try:
        origin_time = event.preferred_origin().time or event.origins[0].time
    except AttributeError:
        # If there isn't an origin time, use the start of the stream
        origin_time = min([pick.time for pick in event.picks])
    t1 = origin_time - length / 10
    t2 = t1 + length
    bulk = []
    for pick in event.picks:
        if all_channels and pick.waveform_id.channel_code:
            channel = "{0}?".format(pick.waveform_id.channel_code[0:2])
        else:
            channel = pick.waveform_id.channel_code or "*"
        chan_info = (pick.waveform_id.network_code or "*", 
                     pick.waveform_id.station_code or "*",
                     pick.waveform_id.location_code or "*", channel,
                     t1, t2)
        if ignore_rotated and channel[-1] in ["R", "T"]:
            continue
        if chan_info not in bulk:
            bulk.append(chan_info)
    try:
        st = client.get_waveforms_bulk(bulk)
    except (FDSNNoDataException, FDSNException) as e:
        print("No data exception - trying individual channels")
        st = Stream()
        for chan_info in bulk:
            try:
                st += client.get_waveforms(
                    network=chan_info[0], station=chan_info[1],
                    location=chan_info[2], channel=chan_info[3],
                    starttime=chan_info[4], endtime=chan_info[5])
                print("Downloaded for {0}.{1}.{2}.{3}".format(*chan_info[0:4]))
            except FDSNNoDataException:
                print("No data for {0}.{1}.{2}.{3}".format(*chan_info[0:4]))
    if filt:
        st.detrend().filter('bandpass', freqmin=filt[0], freqmax=filt[1])
    if return_stream:
        return plot_event(event, st, length=length, size=size, fig=fig), st
    return plot_event(event, st, length=length, size=size, fig=fig)


def plot_event(event, st, length=60., size=(10.5, 10.5), fig=None):
    """
    Plot the waveforms for an event with pick and calculated arrival times.

    :type event: `obspy.core.event.Event`
    :param event: Event to plot
    :type st: `obspy.core.stream.Stream`
    :param st: Obspy Stream for this event
    :type length: float
    :param length: Length to plot, from origin time.
    """
    import matplotlib.pyplot as plt

    try:
        origin_time = event.preferred_origin().time or event.origins[0].time
    except AttributeError:
        # If there isn't an origin time, use the start of the stream
        origin_time = st[0].stats.starttime
    if length:
        st = st.slice(origin_time, origin_time + length)
    # Trim the event around the origin time
    fig = fig or plt.figure()
    fig.set_size_inches(size)
    axes = fig.subplots(len(st), 1, sharex=True)
    if len(st) == 1:
        axes = [axes]
    lines, labels = ([], [])
    min_x = []
    max_x = []
    for ax, tr in zip(axes, st):
        ax.cla()
        picks, arrivals = ([], [])
        for pick in event.picks:
            if pick.waveform_id.station_code == tr.stats.station:
                picks.append(pick)
        try:
            origin = event.preferred_origin() or event.origins[0]
            for arrival in origin.arrivals:
                if arrival.pick_id.get_referred_object(
                        ).waveform_id.station_code == tr.stats.station:
                    arrivals.append(arrival)
        except IndexError:
            pass
        lines, labels, chan_min_x, chan_max_x = _plot_channel(
            ax=ax, tr=tr, picks=picks, arrivals=arrivals, lines=lines, 
            labels=labels)
        min_x.append(chan_min_x)
        max_x.append(chan_max_x)
    axes[-1].set_xlim([np.min(min_x), np.max(max_x)])
    axes[-1].set_xlabel("Time")
    fig.subplots_adjust(hspace=0, wspace=0)
    fig.legend(lines, labels)
    return fig


def _plot_channel(ax, tr, picks=[], arrivals=[], lines=[], labels=[]):
    """ Plot a single channel into an axis object. """
    x = np.arange(0, tr.stats.endtime - tr.stats.starttime + tr.stats.delta,
                  tr.stats.delta)
    y = tr.data
    if len(x) > len(y):
        x = x[0:len(y)]
    elif len(x) < len(y):
        last_x = x[-1]
        for i in range(len(y) - len(x)):
            x.append(last_x + (tr.stats.delta * i))
    x = np.array([(tr.stats.starttime + _x).datetime for _x in x])
    min_x, max_x = (x[0], x[-1])
    ax.plot(x, y, 'k', linewidth=1.2)
    for pick in picks:
        if not pick.phase_hint:
            pcolor = 'k'
            label = 'Unknown pick'
        elif 'P' in pick.phase_hint.upper():
            pcolor = 'red'
            label = 'P-pick'
        elif 'S' in pick.phase_hint.upper():
            pcolor = 'blue'
            label = 'S-pick'
        else:
            pcolor = 'k'
            label = 'Unknown pick'
        line = ax.axvline(x=pick.time.datetime, color=pcolor, linewidth=2,
                          linestyle='--', label=label)
        if label not in labels:
            lines.append(line)
            labels.append(label)
        if pick.time.datetime > max_x:
            max_x = pick.time.datetime
        elif pick.time.datetime < min_x:
            min_x = pick.time.datetime
    for arrival in arrivals:
        if not arrival.phase:
            pcolor = 'k'
            label = 'Unknown arrival'
        elif 'P' in arrival.phase.upper():
            pcolor = 'red'
            label = 'P-arrival'
        elif 'S' in arrival.phase.upper():
            pcolor = 'blue'
            label = 'S-arrival'
        else:
            pcolor = 'k'
            label = 'Unknown arrival'
        arrival_time = (
            arrival.pick_id.get_referred_object().time + arrival.time_residual)
        line = ax.axvline(x=arrival_time.datetime, color=pcolor, linewidth=2,
                          linestyle='-', label=label)
        if label not in labels:
            lines.append(line)
            labels.append(label)
        if arrival_time.datetime > max_x:
            max_x = arrival_time.datetime
        elif arrival_time.datetime < min_x:
            min_x = arrival_time.datetime
    ax.set_ylabel(tr.id, rotation=0, horizontalalignment="right")
    ax.yaxis.set_ticks([])
    return lines, labels, min_x, max_x


if __name__ == "__main__":
    from obspy.clients.fdsn import Client
    import matplotlib.pyplot as plt
    import argparse

    parser = argparse.ArgumentParser(
        description="Download and plot waveforms and picks for an FDSN event")
    parser.add_argument(
        "-c", "--client", help="Client name parsable by obspy", required=False,
        default="GEONET", type=str)
    parser.add_argument(
        "-i", "--eventid", help="FDSN event id", required=True, type=str)
    parser.add_argument(
        "-l", "--length", help="Length of data stream in seconds", 
        required=False, type=float, default=60.)
    parser.add_argument(
        "-n", "--n-stations", help="Limit the number of stations to download",
        required=False, type=int, default=None)
    parser.add_argument(
        "-a", "--all-channels", 
        help="Flag to download all (not just the picked) channels", 
        required=False, action="store_true")
    
    args = vars(parser.parse_args())

    client = Client(args["client"])
    event = client.get_events(eventid=args["eventid"])[0]

    if args["n_stations"] is not None:
        sorted_picks = sorted(event.picks, key=lambda p: p.time)
        stations_to_use = []
        n_stations = 0
        for pick in sorted_picks:
            if pick.waveform_id.station_code not in stations_to_use:
                stations_to_use.append(pick.waveform_id.station_code)
                n_stations += 1
            if n_stations == args["n_stations"]:
                break
        else:
            print("Using all stations")
        event.picks = [p for p in event.picks 
                       if p.waveform_id.station_code in stations_to_use]
    
    fig, st = plot_event_from_client(
        event=event, client=client, length=args["length"], return_stream=True,
        all_channels=args["all_channels"])
    st.write("{0}.ms".format(args["eventid"]), format="MSEED")
    event.write("{0}.xml".format(args["eventid"]), format="QUAKEML")
    plt.show()