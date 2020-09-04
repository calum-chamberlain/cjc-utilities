from obspy.core.event import (
    Origin, Magnitude, ResourceIdentifier, OriginUncertainty, OriginQuality)
from obspy import UTCDateTime

"""
2009  1 18  0 49 21.000         1 -42.23833  173.79283  19.900  2.85       1   12359       1     0     0     0  0.00  0.00  -1.000  -1.000  -1.000   -42.23833  173.79283  19.900
"""

FORMATTER = {
    "year": (0, int),
    "month": (1, int),
    "day": (2, int),
    "hour": (3, int),
    "minute": (4, int),
    "second": (5, float),
    "eventid": (6, int),
    "latitude": (7, float),
    "longitude": (8, float),
    "depth": (9, float),
    "magnitude": (10, float),
    "qID": (11, int),
    "cID": (12, int),
    "nbranch": (13, int),
    "qnpair": (14, int),
    "qndiffP": (15, int),
    "qndiffS": (16, int),
    "rmsP": (17, float),
    "rmsS": (18, float),
    "eh": (19, float),
    "ez": (20, float),
    "et": (21, float),
    "latitude_original": (22, float),
    "longitude_original": (23, float),
    "depth_origins": (24, float)}


def growclust_line_to_origin(line: str) -> [str, Origin]:
    line = line.split()
    deserialized = {key: val[1](line[val[0]]) for key, val in FORMATTER.items()}
    if deserialized["eh"] == -1.:
        # print(f"Event {deserialized['eventid']} not relocated")
        return deserialized["eventid"], None
    origin_time = UTCDateTime(
        year=deserialized["year"], month=deserialized["month"], 
        day=deserialized["day"], hour=deserialized["hour"], 
        minute=deserialized["minute"]) + deserialized["second"]
    try:
        p_standard_error = deserialized["rmsP"] / deserialized["qndiffP"]
    except ZeroDivisionError:
        p_standard_error = 0.0
    try:
        s_standard_error = deserialized["rmsS"] / deserialized["qndiffS"]
    except ZeroDivisionError:
        s_standard_error = 0.0
    origin = Origin(
        latitude=deserialized["latitude"], longitude=deserialized["longitude"],
        depth=deserialized["depth"] * 1000, time=origin_time, 
        method_id=ResourceIdentifier("GrowClust"), 
        time_errors={"uncertainty": deserialized["et"]},
        depth_errors={"uncertainty": deserialized["ez"] * 1000.0},
        time_fixed=False,
        origin_uncertainty=OriginUncertainty(
            horizontal_uncertainty=deserialized["eh"] * 1000.0),
        quality=OriginQuality(
            used_phase_count=deserialized["qndiffP"] + deserialized["qndiffS"],
            standard_error=(
                deserialized["qndiffP"] + deserialized["qndiffS"]) * 
                (p_standard_error + s_standard_error)))
    return deserialized["eventid"], origin


def read_growclust(fname: str, event_mapper: dict = None) -> dict:
    """
    Read growclust origins from a relocated file.

    Parameters
    ----------
    fname:
        File to read from - should be a growclust_cat file.
    event_mapper:
        Event id mapping of {growclust id: desired id}
    
    Returns
    -------
    Dictionary of origins keyed by event id.
    """
    with open(fname, "r") as f:
        lines = f.read().splitlines()
    
    origins = dict()
    for line in lines:
        event_id, growclust_origin = growclust_line_to_origin(line)
        if event_mapper:
            event_id = event_mapper.get(event_id, f"{event_id}_notmapped")
        origins.update({event_id: growclust_origin})
    return origins
