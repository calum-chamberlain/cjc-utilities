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

from concurrent.futures import ThreadPoolExecutor
from typing import Union, Iterable
from obspy import UTCDateTime, Stream, read


GEONET_AWS = "s3://geonet-open-data/waveforms/miniseed/"

PATH_STRUCT = ("{date.year}/{date.year}.{date.julday:03d}/"
               "{station}.{network}/{date.year}.{date.julday:03d}."
               "{station}.{location}-{channel}.{network}.D")

Logger = logging.getLogger(__name__)


class AWSClient:
    def __init__(
        self, 
        base_url: str = GEONET_AWS, 
        path_structure: str = PATH_STRUCT,
        file_length_seconds: float = 86400
    ):
        self.base_url = base_url
        self.path_structure = path_structure
        self.file_length_seconds = file_length_seconds

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
        threaded: bool = True,
        **kwargs
    ):
        remote_paths = []
        for _bulk in bulk:
            remote_paths.extend(self._bulk_to_remote(*_bulk, **kwargs))

        temp_dir = tempfile.mkdtemp(prefix="AWSClient_")
        self._download_remote_paths(remote_paths, threaded=threaded, temp_dir=temp_dir)
        st = Stream()
        local_paths = glob.glob(temp_dir + "/*")
        for local_path in local_paths:
            st += read(local_path)
        st.merge()
        
        trimmed_st = Stream()
        for _bulk in bulk:
            n, s, l, c, _start, _end = _bulk
            trimmed_st += st.select(n, s, l, c).slice(_start, _end).copy()

        shutil.rmtree(temp_dir)

        return trimmed_st.merge()
    
    def _download_remote_paths(
        self, 
        remote_paths: Iterable, 
        threaded: bool,
        temp_dir: str
    ):
        
        def download(remote):
            # TODO: Cope with wildcards - use include and exclude? https://www.janbasktraining.com/community/aws/how-to-use-aws-s3-cp-wildcards-to-copy-group-of-files-in-aws-cli
            ret = subprocess.run(
                ["aws", "s3", "cp", "--no-sign-request", remote, temp_dir],
                capture_output=True)
            return ret
            
        if threaded:
            with ThreadPoolExecutor() as executor:
                results = executor.map(download, remote_paths)
        else:
            results = [download(remote) for remote in remote_paths]
        
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
            remote_paths.append(self._make_remote_path(
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
        if "*" in nslc or "?" in nslc:
            raise NotImplementedError(f"Wildcards found - not designed for this: {nslc}")
        remote_path = self.base_url + self.path_structure
        remote_path = remote_path.format(
            network=network, station=station, location=location, channel=channel,
            date=date)
        return remote_path

