"""
Script for QC-ing earthquake detections.

Calum Chamberlain: 16/07/2018
Updated: 13/05/2024
"""
import glob
import os
import numpy as np
import matplotlib.pyplot as plt

from typing import Union

from obspy import read_events, Catalog, Stream, read
from obspy.core.event import Event
from obspy.clients.fdsn import Client

from obsplus.bank import WaveBank


# Global figure naming to get consistent figure lookup
FIG_NAME = "{plot_dir}/{ori_time}_{rid}.png"
DPI = 1200

def manual_check(
    catalog: Catalog, 
    client: Union[Client, WaveBank],
    plot_dir: str = None,
    checked_dict: dict = None,
    save_progress: bool = True,
    fig: plt.Figure = None,
    check: bool = True,
    length: float = 120.0,
) -> dict:
    """ Perform manual checks of the detections. """
    import json
    fig = None
    checked_dict = checked_dict or dict()

    if catalog is None:
        print(f"No catalog given, just checking events in {plot_dir}")
        plot_files = glob.glob(f"{plot_dir}/*.png")
        total_events = len(plot_files)
        for i, plot_file in enumerate(plot_files):
            try:
                status, fig = check_event(
                        fig=fig, check=check, event_no=i, total_events=total_events,
                        plot_file=plot_file)
            except Exception as e:
                print(f"Could not check event due to {e}")
        return checked_dict

    catalog.events.sort(key=lambda e: e.picks[0].time)  # Sort by pick time
    total_events = len(catalog)
    for i, event in enumerate(catalog):
        if event.resource_id.id.split('/')[-1] in checked_dict.keys():
            continue
        status, fig = check_event(
            event=event, client=client, fig=fig, check=check,
            event_no=i, total_events=total_events, plot_dir=plot_dir)
        if check:
            checked_dict.update({event.resource_id.id.split('/')[-1]: status})
            if save_progress:
                with open("manual_check_progress.json", "w") as f:
                    json.dump(checked_dict, f)
    return checked_dict, fig


def check_event(
    event: Event = None, 
    client: Union[Client, WaveBank] = None,
    plot_dir: str = None,
    plot_file: str = None,
    fig: plt.Figure = None, 
    min_stations: int = 4,
    check: bool = True,
    event_no: int = 1,
    total_events: int = 1,
    min_p_picks: int = 0,
    length: float = 120.0,
) -> str:
    """ Check a single event. """
    import matplotlib.image as img
    from cjc_utilities.plot_event.plot_event import plot_event_from_client

    if event:
        try:
            event_time = (event.preferred_origin() or event.origins[0]).time
        except IndexError:
            event_time = sorted(event.picks, key=lambda p: p.time)[0].time
    else:
        event_time = None
    status_mapper = {"G": "good", "B": "bad", "U": "Undecided"}
    status = None
    fig = fig or plt.figure()
    if fig is not None:
        fig.clf()
    fig_name = None
    if plot_dir and not plot_file:
        fig_name = FIG_NAME.format(
            plot_dir=plot_dir, 
            ori_time=(event.preferred_origin() or event.origins[-1]).time,
            rid=event.resource_id.id.split('/')[-1])
    elif plot_file:
        fig_name = plot_file
    if fig_name:
        print(f"Looking for figure: {fig_name}")
    if fig_name and os.path.isfile(fig_name):
        # Reuse old figure
        print("Reusing old figure")
        old_plot = img.imread(fig_name)
        ax = fig.gca()
        ax.imshow(old_plot)
        ax.axis("off")
    else:        
        print("Downloading new data")
        picked_stations = {pick.waveform_id.station_code for pick in event.picks}
        if len(picked_stations) < min_stations:
            return "bad", fig
        p_picks = [p for p in event.picks if p.phase_hint[0] == "P"]
        if len(p_picks) < min_p_picks:
            print("Fewer than {0} P picks for event {1}, "
                  "considered bad".format(min_p_picks, event_no))
            return "bad", fig
        fig = plot_event_from_client(
                client=client, event=event, fig=fig, size=(8.5, 8.5),
                length=length)
    fig.canvas.draw()
    fig.show()
    status = None
    while not status:
        status_short = input(
            "Event {0} of {1} at {2}\tType your verdict: "
            "(G=good, B=bad, U=undecided)".format(
                event_no, total_events, event_time))
        if status_short.upper() in status_mapper.keys():
            status = status_mapper[status_short.upper()]
        else:
            print("Unknown status {0}, try again".format(status_short))
            continue
    return status, fig


def plot_for_all_events(
    catalog: Catalog, 
    client: Union[Client, WaveBank], 
    plot_dir: str,
    length: float = 120.0,
    overwrite: bool = False,
):
    from cjc_utilities.plot_event.plot_event import plot_event_from_client

    fig = plt.figure()
    for i, event in enumerate(catalog):
        fig_name = FIG_NAME.format(
            plot_dir=plot_dir, 
            ori_time=(event.preferred_origin() or event.origins[-1]).time,
            rid=event.resource_id.id.split('/')[-1])
        if os.path.isfile(fig_name) and not overwrite:
            continue
        fig.clf()
        print("Working on event {0} of {1}".format(i, len(catalog)))
        fig = plot_event_from_client(
                client=client, event=event, fig=fig, size=(8.5, 8.5),
                length=length)
        fig.savefig(fig_name, dpi=DPI)

    return


def main():
    import json
    import glob

    from argparse import ArgumentParser

    parser = ArgumentParser(description="QC detections")

    parser.add_argument("-p", "--plot", action="store_true", 
                        help="Plot all events in advance, requires --plot-dir argument")
    parser.add_argument("--plot-dir", type=str, help="Directory to save plots to, or read from",
                        default=None)
    parser.add_argument("-c", "--check", action="store_true", help="Check detections")
    parser.add_argument("--catalog", type=str, required=False, help="Catalog to read events from")
    parser.add_argument("--client", type=str, help="FDSN client URL or ID to get waveforms from", default=None)
    parser.add_argument("--wavebank", type=str, help="WaveBank path to get data from", default=None)
    parser.add_argument("--length", type=float, help="Length of waveform to plot", default=120.0)
    parser.add_argument("--overwrite", action="store_true", help="Overwrite old figures")

    args = parser.parse_args()

    if args.plot:
        if args.plot_dir:
            assert os.path.isdir(args.plot_dir), f"{args.plot_dir} does not exist"
        else:
            raise IOError("Requires --plot-dir")

    assert args.wavebank or args.client, "Requires either --client or --wavebank"
    print(args)
    if args.client:
        client = Client(args.client)
    elif args.wavebank:
        client = WaveBank(args.wavebank)


    if args.plot:
        cat = read_events(args.catalog)
        cat.events.sort(key=lambda ev: ev.origins[-1].time)
        print("There are {0} events in this file".format(len(cat)))
        fig = None
        plot_for_all_events(
                catalog=cat, client=client, plot_dir=args.plot_dir, 
                length=args.length, overwrite=args.overwrite)

    if args.check:
        print("Running manual check")
        if os.path.isfile("manual_check_progress.json"):
            with open("manual_check_progress.json", "rb") as f:
                check_dict = json.load(f)
        else:
            check_dict = None
        fig = None
        if args.catalog:
            cat = read_events(args.catalog)             
            cat.events.sort(key=lambda ev: ev.origins[-1].time)
            print("There are {0} events in this file".format(len(cat)))
        else:
            assert os.path.isdir(args.plot_dir), f"{args.plot_dir} does not exist"
            cat = None
        check_dict, fig = manual_check(
            catalog=cat, client=client,
            checked_dict=check_dict, save_progress=True, fig=fig,
            plot_dir=args.plot_dir, length=args.length)
        with open("manual_check_complete.json", "w") as f:
            json.dump(check_dict, f)


if __name__ == "__main__":
    main()

