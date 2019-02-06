"""
Simple function for reading REST origins
"""


def read_hypoDD(reloc, mapper=None):
    """
    Read origins from a hypoDD reloc file.

    Note: Not all attributes are read in.

    :type reloc: str
    :param reloc: File to read from
    :type mapper: dict
    :param mapper: Dictionary of event-ids to map to as {hypodd_id: event_id}

    :returns: dict of origins keyed by event id.
    """
    origins = {}
    with open(reloc, 'rb') as f:
        origin_lines = f.read().decode("UTF8").splitlines()
    for origin_line in origin_lines:
        hypoDD_id, origin = _read_hypoDD_origin(origin_line)
        if mapper is not None:
            event_id = mapper.get(hypoDD_id, hypoDD_id)
        origins.update({event_id: origin})
    return origins


def _read_hypoDD_origin(origin_line):
    from obspy import UTCDateTime
    from obspy.core.event import Origin, ResourceIdentifier
    
    l = origin_line.split()
    # unpack
    (event_id, lat, lon, dep, x, y, z, ex, ey, ez, yr, mo, dy, hr, mi, sc, 
     mag, nccp, nccs, nctp, ncts, rcc, rct, cid) = l
    origin_time = UTCDateTime(
        int(yr), int(mo), int(dy), int(hr), int(mi), float(sc))

    origin = Origin(
        time=origin_time, latitude=float(lat), longitude=float(lon),
        depth=float(dep) * 1000, method=ResourceIdentifier(id="HypoDD"),
        origin_type="hypocenter")
    return event_id, origin

