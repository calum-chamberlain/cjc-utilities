"""
Functions to halp assess and adjust picks.

"""

from typing import Tuple
import progressbar

from numbers import Number

from obspy import read_events, Trace, Stream
from obspy.core.event import Pick, Event, Catalog
from obspy.clients.fdsn import Client

from matplotlib.lines import Line2D
from matplotlib.pyplot import Figure


# GLOBALS for setting what keys correspond to.
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
INVERSE_KEY_MAPPING = {value: key for key, value in KEY_MAPPING.items()}
POLARITY_KEYS = {key for key, value in KEY_MAPPING.items()
                 if value in ["positive", "negative", "undecidable"]}
TIME_KEYS = {key for key, value in KEY_MAPPING.items()
             if isinstance(value, Number)}

COLORS = {"P": "red", "S": "yellow"}


def pick_polarity(
    trace: Trace,
    pick: Pick,
    pre_pick: float,
    post_pick: float,
    fig: Figure = None,
):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import timezone

    utc = timezone.utc

    if fig is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    else:
        ax = fig.gca()
        ax.clear()  # Remove the old plot

    trace = trace.slice(pick.time - pre_pick, pick.time + post_pick)
    starttime = trace.stats.starttime
    delta = trace.stats.delta
    phase_type = pick.phase_hint[0]

    ax.plot(trace.times(), trace.data)

    pick_time = pick.time - starttime
    line = ax.add_line(
        Line2D(xdata=[pick_time, pick_time],
               ydata=list(ax.get_ylim()), color=COLORS[phase_type]))
    ax.set_title(f"{trace.id}: pick {INVERSE_KEY_MAPPING['positive']} "
                 f"or {INVERSE_KEY_MAPPING['negative']}, or move "
                 f"{INVERSE_KEY_MAPPING[1]} or {INVERSE_KEY_MAPPING[-1]}. "
                 f"{INVERSE_KEY_MAPPING['next']} to ignore")
    ax.set_xlabel(f"Seconds from {trace.stats.starttime}")
    pol_text = ax.text(0.1, 0.9, f"Polarity: {pick.polarity}", transform=ax.transAxes)

    # Trackers
    fuckup, _quit, polarity_picked, time_adjusted, delete = False, False, False, 0, False

    def updown(event):
        nonlocal fuckup, _quit, polarity_picked, time_adjusted, line, starttime, delta, pol_text, delete
        # Set polarity
        if event.key in POLARITY_KEYS:
            pick.polarity = KEY_MAPPING.get(event.key)
            polarity_picked = True
            pol_text.set_text(f"Polarity: {pick.polarity}")
            pol_text.axes.draw_artist(pol_text)
            pol_text.figure.canvas.draw()
            # plt.close()
            #fig.canvas.stop_event_loop()
            return
        # Adjust the time.
        elif event.key in TIME_KEYS:
            pick.time += (KEY_MAPPING[event.key] * delta)
            time_adjusted += (KEY_MAPPING[event.key] * delta)
            pick_time = pick.time - starttime
            line.set_data([pick_time, pick_time], list(line.axes.get_ylim()))
            line.axes.draw_artist(line)
            line.figure.canvas.draw()
        # Go to next event
        elif event.key == INVERSE_KEY_MAPPING["next"]:
            fig.canvas.stop_event_loop()
            return
        elif event.key == INVERSE_KEY_MAPPING["remove"]:
            delete = True
            fig.canvas.stop_event_loop()
            return
        # Go back one event
        elif event.key == INVERSE_KEY_MAPPING["fuckup"]:
            # Go back to the previous one...
            print("Fucked-up, going back")
            fuckup = True
            # plt.close()
            fig.canvas.stop_event_loop()
            return
        # Exit the whole sehbang
        elif event.key == "q":
            print("Quitting")
            _quit = True
            fig.canvas.stop_event_loop()
            return
        else:
            print(f"{event.key} is not bound")

    # Attach responder
    cid = fig.canvas.mpl_connect('key_press_event', updown)

    fig.show()
    fig.canvas.draw()  # Redraw
    fig.canvas.start_event_loop()
    # Picking happens here
    fig.canvas.mpl_disconnect(cid)
    # If we want to delete the pick
    if delete:
        return None, fuckup, _quit
    if not fuckup and not _quit:
        if polarity_picked:
            print(f"{trace.id} picked as {pick.polarity}")
        if time_adjusted:
            print(f"{trace.id} adjusted by {time_adjusted}")
    return pick, fuckup, _quit


def check_event(
    st: Stream,
    event: Event,
    pre_pick: float=1.0,
    post_pick: float=1.0,
):
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")  # Easier on my eyes

    # Keep all the non P-picks
    picks_to_check = [p for p in event.picks if p.phase_hint.startswith(("P", "S"))]
    picks_to_keep = [p for p in event.picks if not p.phase_hint.startswith(("P", "S"))]

    # Get available stations
    stations = {tr.stats.station for tr in st}
    checked_picks = []
    i = 0  # Using a while loop to allow us to go back if we fuckup.
    fig, ax = plt.subplots(figsize=(12, 8))
    while i < len(picks_to_check):
        pick = picks_to_check[i]
        tr = st.select(station=pick.waveform_id.station_code)
        if len(tr) == 0:
            print(f"No data for {pick.waveform_id.station_code}")
            continue
            # pick.polarity = "undecidable"
        for _tr in tr.merge():
            checked_pick, mistake, _quit = pick_polarity(
                _tr, pick.copy(), fig=fig, pre_pick=pre_pick, post_pick=post_pick)
            if checked_pick is None:
                print("Removing pick")
                continue
            if _quit:
                print("Quitting")
                return None
            if mistake is True:
                print("You fucked up - removing last pick")
                if i == 0:
                    print("Already at zeroth pick, try again")
                else:
                    i -= 1  # Go back.
                    checked_picks = checked_picks[0:-1]  # Remove the last pick.
                break
            checked_pick.waveform_id.channel_code = _tr.stats.channel
            checked_pick.waveform_id.network_code = _tr.stats.network
            checked_pick.waveform_id.location_code = _tr.stats.location
            checked_picks.append(checked_pick)
        else:
            # If we don't break, continue
            i += 1
    event_picked = event.copy()
    event_picked.picks = picks_to_keep + checked_picks
    plt.close(fig)
    return event_picked


def pick_geonet_event(eventid: str) -> Event:
    from cjc_utilities.get_data.get_data import get_event_data

    client = Client("GEONET")

    event, st = get_event_data(client=client, eventid=eventid, length=20,
                               all_channels=False, ignore_rotated=True)
    checked_event = check_event(st=st, event=event)
    return checked_event


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eventid", type=str, required=True)
    parser.add_argument("-o", "--outfile", type=str, default="Checked_event.xml")

    args = parser.parse_args()
    checked_event = pick_geonet_event(args.eventid)
    check_event.write(args.outfile, format="QUAKEML")
