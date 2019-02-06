"""
Script to get basic GeoNet info from quakesearch.
"""

import urllib.request
import csv
import logging
from datetime import timedelta
from datetime import datetime as dt

Logger = logging.getLogger("GeoNet-downloader")

BASE_URL = ("https://quakesearch.geonet.org.nz/csv?bbox={0},{1},{2},{3}"
            "&startdate={4}&enddate={5}")


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
