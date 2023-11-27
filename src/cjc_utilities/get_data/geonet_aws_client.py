"""
Python client for AWS data access to GeoNet's open bucket.

Simulates get_waveforms and get_waveforms_bulk in FDSN clients.
"""

import subprocess
import os
import logging
import datetime as dt
import shutil
import tempfile
import glob
import fnmatch
import datetime as dt

from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from typing import Union, Iterable
from obspy import UTCDateTime, Stream, read

import boto3
from botocore import UNSIGNED
from botocore.config import Config


GEONET_AWS = "geonet-open-data"

DAY_STRUCT = "waveforms/miniseed/{date.year}/{date.year}.{date.julday:03d}"
CHAN_STRUCT = ("{station}.{network}/{date.year}.{date.julday:03d}."
              "{station}.{location}-{channel}.{network}.D")

Logger = logging.getLogger(__name__)


class AWSClient:
    def __init__(
        self, 
        bucket_name: str = GEONET_AWS,
        day_structure: str = DAY_STRUCT,
        nslc_structure: str = CHAN_STRUCT,
        file_length_seconds: float = 86400,
        threaded: bool = True,
    ):
        self.bucket_name = bucket_name
        self.day_structure = day_structure
        self.nslc_structure = nslc_structure
        self.file_length_seconds = file_length_seconds
        self.threaded = threaded

    @property
    def _s3(self):
        """ s3 access as property to enforce instantiation for individual threads. """
        # return boto3.client("s3", config=Config(signature_version=UNSIGNED))
        bob = boto3.resource('s3', config=Config(signature_version=UNSIGNED))
        return bob.Bucket(self.bucket_name)

    def get_waveforms(
        self,
        network: str, 
        station: str, 
        location: str, 
        channel: str, 
        starttime: Union[dt.datetime, UTCDateTime], 
        endtime: Union[dt.datetime, UTCDateTime], 
        **kwargs
    ):
        return self.get_waveforms_bulk(
            [(network, station, location, channel, starttime, endtime)], **kwargs)

    def get_waveforms_bulk(
        self,
        bulk: Iterable,
        **kwargs
    ):
        remote_paths = []
        for _bulk in bulk:
            remote_paths.extend(self._bulk_to_remote(*_bulk, **kwargs))

        # Work out the maximum extent of data to read in
        min_start = min(_bulk[4] for _bulk in bulk)
        max_end = max(_bulk[5] for _bulk in bulk)

        temp_dir = tempfile.mkdtemp(prefix="AWSClient_")
        Logger.info(f"Downloading files to {temp_dir}")
        self._download_remote_paths(remote_paths, temp_dir=temp_dir)
        st = Stream()
        for path, _, files in os.walk(temp_dir):
            for f in files:
                Logger.info(f"Reading from {f}")
                try:
                    st += read(os.path.join(path, f), starttime=min_start, endtime=max_end)
                except Exception as e:
                    Logger.error(f"Could not read from {f} due to {e}; skipping")
                    continue
        st.merge()

        trimmed_st = Stream()
        for _bulk in bulk:
            n, s, l, c, _start, _end = _bulk
            trimmed_st += st.select(n, s, l, c).slice(_start, _end).copy()

        shutil.rmtree(temp_dir)
        trimmed_st = trimmed_st.merge().split()

        return trimmed_st
    
    def _download_remote_paths(
        self, 
        remote_paths: Iterable, 
        temp_dir: str
    ):
        
        def download(remote, max_retries=999):
            local_path = os.path.join(temp_dir, remote)
            os.makedirs(os.path.dirname(local_path), exist_ok=True)

            attempt = 0
            while attempt <= max_retries:
                try:
                    Logger.info(f"Downloading {remote} attempt {attempt} of {max_retries}")
                    self._s3.download_file(remote, local_path)
                except Exception as e:
                    Logger.error(f"Could not download {remote} due to {e}")
                    attempt += 1
                    continue
                break
            return
            
        if self.threaded:
            with ThreadPoolExecutor() as executor:
                results = executor.map(download, remote_paths)
        else:
            results = [download(remote) for remote in remote_paths]
        
        for res in results:
            if isinstance(res, Exception):
                raise(res)
            elif isinstance(res, str):
                Logger.warning(res)
        return    

    def _bulk_to_remote(
        self,
        network: str,
        station: str,
        location: str,
        channel: str,
        starttime: Union[dt.datetime, UTCDateTime],
        endtime: Union[dt.datetime, UTCDateTime],
        **kwargs
    ):

        if not isinstance(starttime, UTCDateTime):
            starttime = UTCDateTime(starttime)
        if not isinstance(endtime, UTCDateTime):
            endtime = UTCDateTime(endtime)
        
        _endtime = endtime + self.file_length_seconds

        remote_paths = []
        date = starttime - self.file_length_seconds
        while date < _endtime:
            remote_paths.extend(self._make_remote_path(
                network=network, station=station, location=location,
                channel=channel, date=date, **kwargs))
            date += dt.timedelta(seconds=self.file_length_seconds)

        return remote_paths

    def _make_remote_path(
        self,
        network: str,
        station: str,
        location: str,
        channel: str,
        date: UTCDateTime,
        **kwargs
    ):
        nslc = ".".join([network, station, location, channel])
        day_contents = _list_day_files(
            bucket=self._s3, day_structure=self.day_structure, date=date.datetime)

        path_structure = '/'.join((self.day_structure, self.nslc_structure))
        nslc_path = path_structure.format(
            network=network, station=station, location=location, channel=channel,
            date=date)
        remote_paths = fnmatch.filter(day_contents, nslc_path)

        Logger.debug(f"{nslc} at {date}: {remote_paths}")
        return remote_paths


@lru_cache(10)
def _list_day_files(bucket, day_structure: str, date: dt.datetime):
    # Note we need to use datetime input to make this hashable.
    # List the contents for this date, then search in that list for matches
    Logger.info(f"Getting contents of {bucket} for {date}")
    day_contents = bucket.objects.filter(
        Prefix=day_structure.format(date=UTCDateTime(date)))
    files = {c.key for c in day_contents if not c.key.endswith('/')}
    return files
