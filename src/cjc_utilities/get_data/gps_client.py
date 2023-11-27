"""
Simple client for GPS data from GeoNet.

"""

import requests
import numpy as np
import json
import copy

import matplotlib.pyplot as plt

from datetime import datetime as dt

from typing import List, Union, Iterable


GEONET_FITS = "http://fits.geonet.org.nz/observation"
GEONET_FITS_STATION = "http://fits.geonet.org.nz/site?"
DEFAULT_COLORS = {"u": "grey", "e": "lightcoral", "n": "teal"}


class GPSData():
    def __init__(
        self, 
        receiver: str, 
        component: str,
        times: np.ndarray,
        observations: np.ndarray, 
        errors: np.ndarray,
    ):
        """ Holder for GPS data. """
        assert times.shape == observations.shape, "Times is not shaped like observations"
        assert observations.shape == errors.shape, "Observations is not shaped like errors"
        self.receiver = receiver
        self.component = component
        # Sort everything
        sort_mask = np.argsort(times)
        self.times = times[sort_mask]
        self.observations = observations[sort_mask]
        self.errors = errors[sort_mask]

    def copy(self):
        return copy.deepcopy(self)

    def __repr__(self):
        return f"GPSData(receiver: {self.receiver}, component: {self.component}, ...)"

    def _time_mask(self, starttime: dt = None, endtime: dt = None) -> np.array:
        starttime = starttime or self.times[0]
        endtime = endtime or self.times[-1]
        mask = np.where(
            np.logical_and(self.times >= starttime, self.times <= endtime))
        return mask


    def trim(self, starttime: dt = None, endtime: dt = None):
        """Trim data between start and end times - works in place. """
        mask = self._time_mask(starttime, endtime)
        self.times = self.times[mask]
        self.observations = self.observations[mask]
        self.errors = self.errors[mask]
        return self

    def detrend(self, starttime: dt = None, endtime: dt = None):
        """ Detrend data - works in place. """
        mask = self._time_mask(starttime, endtime)
        data, times = self.observations[mask], self.times[mask]
        x1, x2 = data[0], data[-1]
        t1, t2 = times[0], times[-1]
        dx = x2 - x1
        dt = (t2 - t1).total_seconds()
        gradient = dx / dt
        seconds = np.array([(t - times[0]).total_seconds() for t in self.times])
        offsets = (gradient * seconds) + x1
        self.observations -= offsets
        return self

    def plot(
        self, 
        show: bool = False, 
        ax: plt.Axes = None,
        **kwargs
    ) -> plt.Figure:
        """Plot GPS data. Can be given an existing axes object to plot onto. """
        if ax:
            fig = ax.figure
        else:
            fig, ax = plt.subplots()

        internal_kwargs = kwargs.copy()
        color = internal_kwargs.pop("color", DEFAULT_COLORS[self.component])

        default_kwargs = {"markersize": 1.0, "linewidth": 0.5}
        for key, value in default_kwargs.items():
            if key not in internal_kwargs.keys():
                internal_kwargs.update({key: value})

        ax.errorbar(
            self.times, self.observations, yerr=self.errors,
            label=f"{self.receiver} {self.component}", 
            fmt="x", 
            color=color,
            **internal_kwargs)

        if show:
            plt.show()
        return fig


class GPSStation():
    _receiver = None
    def __init__(
        self,
        components: List[GPSData] = [],
        latitude: float = None,
        longitude: float = None,
        elevation: float = None,
    ):
        """ Holder for channels from one station. """
        self.components = []
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        if len(components):
            self._receiver = components[0].receiver
            for component in components:
                self.__iadd__(component)

    def __repr__(self):
        return f"GPSStation(<{len(self.components)} components for {self.receiver}>)"

    @property
    def receiver(self):
        return self._receiver

    def copy(self):
        return copy.deepcopy(self)

    def __iter__(self):
         return list(self.components).__iter__()

    def __iadd__(self, other: Union[GPSData, Iterable[GPSData]]):
        if isinstance(other, GPSData):
            other = [other]
        for component in other:
            if not isinstance(component, GPSData):
                raise TypeError(f"Addition of {type(component)} to GPSStation not supported")
            if not self.receiver:
                self._receiver = component.receiver
            assert component.receiver == self.receiver, "Receivers do not match"
            self.components.append(component)
        return self

    def __add__(self, other: Union[GPSData, Iterable[GPSData]]):
        new = self.copy() 
        new += other
        return new

    def __getitem__(self, index):
        return self.components[index]

    def detrend(self, starttime: dt = None, endtime: dt = None):
        for c in self:
            c.detrend(starttime, endtime)

    def trim(self, starttime: dt = None, endtime: dt = None):
        for c in self:
            c.trim(starttime, endtime)

    def plot(self, seperate_channels: bool = False, 
             show: bool = False, ax: plt.Axes = None, **kwargs) -> plt.Figure:
        if not ax:
            if seperate_channels:
                fig, ax = plt.subplots(len(self.components), 1, sharex=True)
                fig.subplots_adjust(hspace=0)
            else:
                fig, ax = plt.subplots()
                ax = [ax for _ in self.components]
        else:
            ax = [ax for _ in self.components]
            if seperate_channels:
                print("You have given an axis but asked for seperate channels. This is not supported")

        for component, _ax in zip(self.components, ax):
            component.plot(ax=_ax, show=False, **kwargs)
        fig.legend()
        if show:
            plt.show()
        return fig


def get_gps_location(receiver: str) -> dict:
    parameters = {"siteID": receiver}
    response = requests.get(GEONET_FITS_STATION, params=parameters)
    assert response.status_code == 200, \
        f"Bad request getting station location from {GEONET_FITS_STATION}. Response code {response.status_code}"
    payload = response.content.decode("utf-8")
    payload = json.loads(payload)
    return {
        "latitude": payload['features'][0]['geometry']['coordinates'][1],
        "longitude": payload['features'][0]['geometry']['coordinates'][0],
        "elevation": payload['features'][0]['properties']['height'],
    }



def get_gps_data(receiver: str, component: str = None) -> GPSStation:
    components = ["u", "n", "e"]
    if component:
        assert component in components, \
            f"Component ({component}) must be in {components})"
        components = [component]

    location = get_gps_location(receiver)

    station = GPSStation(**location)
    for component in components:
        parameters = {"typeID": component, "siteID": receiver}
        response = requests.get(GEONET_FITS, params=parameters)
        assert response.status_code == 200, \
            f"Bad request getting data from {GEONET_FITS}. Response code {response.status_code}"
        payload = response.content.decode("utf-8").split("\n")
        # payload is a csv with header
        payload = [p.split(',') for p in payload]
        # Check that this is what we expect
        assert payload[0][0] == 'date-time', "Unkown format"
        assert len(payload[0]) == 3, "Unknown format"
        times, displacements, errors = zip(*[
            (dt.strptime(p[0], '%Y-%m-%dT%H:%M:%S.%fZ'),
            float(p[1]), float(p[2])) for p in payload[1:-1]])
        station += GPSData(
            receiver=receiver, component=component, times=np.array(times),
            observations=np.array(displacements), 
            errors=np.array(errors))
    return station


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Get and plot GNSS time series")

    parser.add_argument("-r", "--receiver", type=str, required=True,
                        help="Station/receiver name to get data for")
    parser.add_argument("-d", "--detrend", action="store_true",
                        help="Set to detrend data")
    parser.add_argument("-s", "--seperate-channels", action="store_true",
                        help="Plot channels on their own axes")
    parser.add_argument("--starttime", type=str, required=False, default=None,
                        help="Time to plot data from (YYYY-MM-DD)")
    parser.add_argument("--endtime", type=str, required=False, default=None,
                        help="Time to plot data to (YYYY-MM-DD)")

    args = parser.parse_args()

    starttime, endtime = None, None
    if args.starttime:
        starttime = dt.strptime(args.starttime, "%Y-%m-%d")
    if args.endtime:
        endtime = dt.strptime(args.endtime, "%Y-%m-%d")

    station = get_gps_data(args.receiver)
    station.trim(starttime=starttime, endtime=endtime)
    title = f"GNSS position for {args.receiver}."
    if args.detrend:
        station.detrend()
        title += " Detrended."
    fig = station.plot(show=False, seperate_channels=args.seperate_channels)
    fig.suptitle(title)
    plt.show()
