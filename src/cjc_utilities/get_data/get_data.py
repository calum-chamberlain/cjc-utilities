"""
Script to get data from known networks on the IRIS DMC and write miniseed \
files in hour-long, single-channel files in year/julian day directories.

Author: Calum Chamberlain 2016
"""

from obspy.clients.fdsn import Client
from obspy.core.event import Event
from obspy import Stream, Inventory

from typing import List, Union, Tuple



def get_event_data(
    client: Client,
    eventid: Union[str, None] = None,
    event: Union[Event, None] = None,
    length: float = 60.0,
    all_channels: bool = False,
    ignore_rotated: bool = True,
    default_channels: List[str] = ["E??", "H??", "S??"],
    start_at_origin: bool = True,
    return_inv: bool = False,
) -> Union[Tuple[Event, Stream], Tuple[Event, Stream, Inventory]]:
    """
    Get data for an event

    """

    from obspy import Stream
    from obspy.clients.fdsn.client import FDSNNoDataException, FDSNException

    event = event or client.get_events(eventid=eventid)[0]

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
            channel = pick.waveform_id.channel_code or default_channels
        if not start_at_origin:
            # Start 10% of length before pick time
            t1 = pick.time - (0.1 * length)
            t2 = t1 + length
        chan_info = (pick.waveform_id.network_code or "*",
                     pick.waveform_id.station_code or "*",
                     pick.waveform_id.location_code or "*", channel,
                     t1, t2)
        if ignore_rotated and channel[-1] in ["R", "T"]:
            continue
        if isinstance(chan_info[3], list):
            for _chan in chan_info[3]:
                _chan_info = (
                    chan_info[0], chan_info[1], chan_info[2], _chan,
                    chan_info[4], chan_info[5])
                if _chan_info not in bulk:
                    bulk.append(_chan_info)
        elif chan_info not in bulk:
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
    if return_inv:
        inv = client.get_stations_bulk(bulk, level="response")
        return event, st, inv
    return event, st


def get_data(network, starttime, endtime, outdir='.'):
    """
    Function to download all data for a given network.

    Will download all stations available between start and end times given. \
    It will create day-long files, and will download data as day-long segments.

    :type network: str
    :param network: Network code
    :type starttime: UTCDateTime
    :param starttime: Time to begin donloading from, will use the date
    :type endtime: UTCDateTime
    :param endtime: Time to end download, will use the date
    :type outdir: str
    :param outdir: Path to write to, will write in Y????/R???.01 directories \
        within this head directory.

    .. note:: This function is slow and doesn't cope with errors, suggest \
        using the mass-downloader.
    """
    import obspy
    import os
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn import Client
    else:
        from obspy.fdsn import Client

    kdays = (endtime.datetime - starttime.datetime).days + 1
    for i in range(kdays):
        t1 = starttime + (86400 * i)
        t2 = t1 + 86400
        bulk_info = [(network, 'WZ01', '*', '*', t1, t2)]

        client = Client('IRIS')
        print('Downloading data for day ' + str(t1))
        print(bulk_info)
        st = client.get_waveforms_bulk(bulk_info)
        # st.plot()
        # Now you should split the data into traces by station and channel and
        # Save them into Y????/R???.01 directories, named something useful,
        # See the SAMBA_archive for details.
        for tr in st:
            f_name = '.'.join(tr.stats.station, tr.stats.network,
                              tr.stats.location, tr.stats.channel,
                              tr.stats.starttime.strftime('%Y.%j'))
            path = outdir + tr.stats.starttime.strftime('/Y%Y/R%j.01/')
            if not os.path.isdir(path):
                os.path.makedirs(path)
            tr.write(path + f_name, format='MSEED', encoding='STEIM2')


def get_mass_data(network, stations, starttime, endtime, outdir='.', lat=0, lon=0,
                  minrad=0, maxrad=180, providers=['IRIS']):
    """
    Use obspy's massdownloader to download lots of data, stores as day-long \
    miniseed files using IRIS DMC naming conventions.

    :type network: str
    :param network: Network code
    :type starttime: UTCDateTime
    :param starttime: Time to begin donloading from, will use the date
    :type endtime: UTCDateTime
    :param endtime: Time to end download, will use the date
    :type outdir: str
    :param outdir: Path to write to, will write in Y????/R???.01 directories \
        within this head directory.
    :type lat: float
    :param lat: Origin latitude
    :type lon: float
    :param lon: Origin longitude
    :type minrad: float
    :param minrad: Minumum radius in degrees for stations from lat/lon.
    :type maxrad: float
    :param maxrad: Maximum radius in degrees for stations from lat/lon
    :type providers: list
    :param providers: List of providers to query.  Default is IRIS. Can parse \
        an empty list and will query all providers, but slow.

    .. note:: Currently selects data using a circular domain, default is \
        set to entire globe so that a whole network can be downloaded. \
        Options left in function to allow repurposing.
    """
    def get_mseed_storage(network, station, location, channel,
                          starttime, endtime):
        """Function to define the naming for the stored file.

        .. note:: Can only have these arguments.  As such this needs to be a \
            function defined within the other function to have access to the \
            outdir variable.
        """
        import os
        # Returning True means that neither the data nor the StationXML file
        # will be downloaded.
        # If a string is returned the file will be saved in that location.
        path = os.path.join(outdir,
                            "%s/%s.%s.%s.%s.%s" % (starttime.
                                                   strftime('Y%Y/R%j.01'),
                                                   network, station, location,
                                                   channel, starttime.
                                                   strftime('%Y.%j')))
        if os.path.exists(path):
            return True
        return path

    import obspy
    if int(obspy.__version__.split('.')[0]) >= 1:
        from obspy.clients.fdsn.mass_downloader import CircularDomain, \
            Restrictions, MassDownloader
    else:
        from obspy.fdsn.mass_downloader import CircularDomain, Restrictions, \
            MassDownloader

    domain = CircularDomain(latitude=lat, longitude=lon,
                            minradius=minrad, maxradius=maxrad)
    for station in stations:
        restrictions = Restrictions(
            starttime=starttime,
            endtime=endtime,
            # Chunk it to have one file per day.
            chunklength_in_sec=86400,
            # Considering the enormous amount of data associated with continuous
            # requests, you might want to limit the data based on SEED identifiers.
            # If the location code is specified, the location priority list is not
            # used; the same is true for the channel argument and priority list.
            network=network, station=station, channel="*",
            # The typical use case for such a data set are noise correlations where
            # gaps are dealt with at a later stage.
            reject_channels_with_gaps=False,
            # Same is true with the minimum length. All data might be useful.
            minimum_length=0.0,
            # Guard against the same station having different names.
            minimum_interstation_distance_in_m=0.0,
            sanitize=True)

        # Restrict the number of providers if you know which serve the desired
        # data. If in doubt just don't specify - then all providers will be
        # queried.
        mdl = MassDownloader(providers=providers)
        mseed_storage = get_mseed_storage
        #  + "/{station}.{network}.{location}.{channel}.{starttime}")
        mdl.download(domain, restrictions, mseed_storage=mseed_storage,
                     stationxml_storage="stations", threads_per_client=3)


def main():
    import sys
    from obspy import UTCDateTime
    if not len(sys.argv) == 5:
        raise IOError('Insufficient arguments, need network code, ' +
                      'start-time and end-time, which need to be ' +
                      'convertable to UTCDateTime, and output directory.')
    network = sys.argv[1]
    starttime = UTCDateTime(sys.argv[2])
    endtime = UTCDateTime(sys.argv[3])
    outdir = sys.argv[4]
    get_mass_data(network, starttime, endtime, outdir)


if __name__ == "__main__":
    main()
