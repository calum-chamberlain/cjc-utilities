"""
Script to get basic GeoNet info from quakesearch.
"""

import urllib.request
import csv
import logging
from datetime import timedelta
from datetime import datetime as dt
from obspy import UTCDateTime

Logger = logging.getLogger("GeoNet-downloader")

BASE_URL = ("https://quakesearch.geonet.org.nz/csv?bbox={0},{1},{2},{3}"
            "&startdate={4}&enddate={5}")


def quake_search(    
    min_latitude=-49.0, max_latitude=-30.0,
    min_longitude=164.0, max_longitude=182.0,
    min_magnitude=0.0, max_magnitude=9.0,
    min_depth=0.0, max_depth=100.0,
    start_time=UTCDateTime(1960, 1, 1),
    end_time=UTCDateTime(2020, 1, 1),
):
    """
    Get a dataframe of the eatrhquakes in the GeoNet catalogue.
    
    Parameters
    ----------
    min_latitude
        Minimum latitude in degrees for search
    max_latitude
        Maximum latitude in degrees for search
    min_longitude
        Minimum longitude in degrees for search
    max_longitude
        Maximum longitude in degrees for search
    min_depth
        Minimum depth in km for search
    max_depth
        Maximum depth in km for search
    min_magnitude
        Minimum magnitude for search
    max_magnitude
        Maximum magnitude for search
    start_time
        Start date and time for search
    end_time
        End date and time for search
        
    Returns
    -------
    pandas.DateFrame of resulting events
    """
    quakes = []
    max_chunk_size = 365 * 86400
    _starttime = start_time
    _endtime = _starttime + max_chunk_size
    kwargs = dict(min_latitude=min_latitude, min_longitude=min_longitude, 
            max_latitude=max_latitude, max_longitude=max_longitude,
            min_depth=min_depth, max_depth=max_depth, 
            min_magnitude=min_magnitude, max_magnitude=max_magnitude)
    while _endtime < end_time:
        quakes.append(_get_geonet_quakes(           
            start_time=_starttime, end_time=_endtime, **kwargs))
        _starttime += max_chunk_size
        _endtime += max_chunk_size
    quakes.append(_get_geonet_quakes(
        start_time=_starttime, end_time=end_time, **kwargs))

    earthquakes = quakes[0]
    for df in quakes[1:]:
        earthquakes = earthquakes.append(df, ignore_index=True)
    return earthquakes


def _get_geonet_quakes(
    min_latitude, max_latitude, min_longitude, max_longitude, min_magnitude, 
    max_magnitude, min_depth, max_depth, start_time, end_time):
    import requests
    import pandas as pd
    import os

    # Convert start_time and end_time to strings
    start_time = start_time.strftime("%Y-%m-%dT%H:%M:%S")
    end_time = end_time.strftime("%Y-%m-%dT%H:%M:%S")
    # Use the more efficient f-string formatting
    query_string = (
        "https://quakesearch.geonet.org.nz/csv?bbox="
        f"{min_longitude},{min_latitude},{max_longitude},"
        f"{max_latitude}&minmag={min_magnitude}"
        f"&maxmag={max_magnitude}&mindepth={min_depth}"
        f"&maxdepth={max_depth}&startdate={start_time}"
        f"&enddate={end_time}")
    print(f"Using query: {query_string}")
    response = requests.get(query_string)
    if not response.ok:
        print(response.content)
        return
    with open(".earthquakes.csv", "wb") as f:
        f.write(response.content)
    earthquakes = pd.read_csv(
        ".earthquakes.csv", 
        parse_dates=["origintime", "modificationtime"],
        dtype={"publicid": str})
    earthquakes = earthquakes.rename(
        columns={" magnitude": "magnitude",
                " latitude": "latitude",
                " depth": "depth"})
    earthquakes = earthquakes.sort_values(by=["origintime"], ignore_index=True)
    os.remove(".earthquakes.csv")
    return earthquakes



def get_geonet_events(startdate, enddate, bbox=(163.96, -49.18, 182.6, -32.3),
                      max_chunk=365, write_out=False):
    if enddate - startdate > timedelta(days=max_chunk):
        Logger.info("Requested window too large, will download in chunks")
        chunks = []
        date = startdate
        while date < enddate:
            _enddate = date + timedelta(days=max_chunk)
            if _enddate > enddate:
                _enddate = enddate
            chunks.append((date, _enddate))
            date += timedelta(days=max_chunk)
    else:
        chunks = [(startdate, enddate)]
    event_info = []
    for chunk in chunks:
        Logger.info(
            "Downloading between {0} and {1}".format(chunk[0], chunk[1]))
        response = urllib.request.urlopen(
            BASE_URL.format(bbox[0], bbox[1], bbox[2], bbox[3],
                            chunk[0].strftime("%Y-%m-%d"),
                            chunk[1].strftime("%Y-%m-%d")))
        Logger.debug(BASE_URL.format(bbox[0], bbox[1], bbox[2], bbox[3],
                                     chunk[0].strftime("%Y-%m-%d"),
                                     chunk[1].strftime("%Y-%m-%d")))
        reader = csv.DictReader(response.read().decode('utf-8').splitlines())
        reader = [line for line in reader]
        for line in reader:
            event_info.append(
                {"latitude": float(line[' latitude']),
                 "longitude": float(line['longitude']),
                 "depth": float(line[' depth']),
                 "magnitude": float(line[' magnitude']),
                 "origin-time": dt.strptime(line['origintime'],
                                            "%Y-%m-%dT%H:%M:%S.%fZ"),
                 "id": line['publicid']})
    if write_out:
        outfile = "get_geonet_events_{0}-{1}.csv".format(
           startdate.strftime("%Y%m%d"), enddate.strftime("%Y%m%d"))
        with open(outfile, "w") as f:
            f.write("PublicID, Latitude (deg), Longitude (deg), Depth (km),"
                    " Magnitude, Origin-time (UTC)\n")
            for event in event_info:
                f.write("{0:.2f}, {1:.2f}, {2:.2f}, {3:.2f}, {4}\n".format(
                    event['id'], event['latitude'], event['longitude'],
                    event['depth'], event['magnitude'], 
                    event['origin-time'].strftime("%Y-%m-%dT%H:%M:%S.%f")))
    else:
        return event_info


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Download lat, lon, depth, mag from GeoNet')
    parser.add_argument(
            '-s', '--startdate',
            help='Startdate as %Y/%m/%d', required=True)
    parser.add_argument(
            '-e', '--enddate',
            help='Enddate as %Y/%m/%d', required=True)
    args = vars(parser.parse_args())
    get_geonet_events(startdate=dt.strptime(args['startdate'], "%Y/%m/%d"),
                      enddate=dt.strptime(args['enddate'], "%Y/%m/%d"),
                      write_out=True)
