"""
Functions for writing to SIMUL PHS format.

"""
from copy import deepcopy
from typing import Tuple
from collections import namedtuple
from obspy.core.event import Catalog, Event

HEAD_STR = (
    "{year:02d}{month:02d}{day:02d} {hour:02d}{minute:02d} {seconds:05.2f}"
    "{lat:3d}{N_S:1s}{lat_min:5.2f} {lon:03d}{E_W:1s}{lon_min:5.2f} "
    "{depth:7.2f}{magnitude:6.1f}  {nobs:3d} {dunno:3d} "
    "{dunnoTwo:3d} {event_id:20s}")
PICK_STR = ("{station:5s} {phase:1s} {weight:1d}   {travel_time:.5f}")


def write_simul(catalog: Catalog, filename: str):
    event_str = "\n".join(
        [_get_simul_string(event) for event in catalog])
    with open(filename, "w") as f:
        f.write(event_str)
    return


def _get_simul_string(event: Event) -> str:
    out = []
    try:
        origin = event.preferred_origin() or event.origins[0]
    except IndexError:
        print("No origin found, skipping")
        return
    try:
        magnitude = event.preferred_magnitude() or event.magnitudes[0]
    except IndexError:
        print("No magnitude found, setting to 0.0")
        _magnitude = namedtuple("fake_magnitude", ["mag"])
        magnitude = _magnitude(0.0)
    lat, n_s, lat_min = _deg_to_deg_dec_min(origin.latitude, "latitude")
    lon, e_w, lon_min = _deg_to_deg_dec_min(origin.longitude, "longitude")
    depth = origin.depth / 1000.

    out.append(HEAD_STR.format(
        year=origin.time.year % 100, month=origin.time.month, 
        day=origin.time.day, hour=origin.time.hour, 
        minute=origin.time.minute,
        seconds=origin.time.second + (origin.time.microsecond / 1e6),
        lat=lat, N_S=n_s, lat_min=lat_min, lon=lon, E_W=e_w, 
        lon_min=lon_min, depth=depth, magnitude=magnitude.mag,
        nobs=len(event.picks), dunno=100, dunnoTwo=0, 
        event_id=event.resource_id.id.split('/')[-1]))
    
    for pick in event.picks:
        out.append(PICK_STR.format(
            station=pick.waveform_id.station_code,
            phase=pick.phase_hint,
            weight=1,
            travel_time=pick.time - origin.time))

    return "\n".join(out)


def _deg_to_deg_dec_min(
    value: float,
    lat_long: str
) -> Tuple[int, str, float]:
    assert lat_long in ["latitude", "longitude"]
    deg, orientation, dec_deg = None, None, None
    # Copy value
    working_value = deepcopy(value)

    if lat_long == "latitude":
        if value < 0.0:
            orientation = "S"
            working_value *= -1
        else:
            orientation = "N"
    else:
        if value < 0.0:
            orientation = "W"
            working_value *= -1
        else:
            orientation = "E"
    deg = int(working_value // 1)
    dec_deg = (working_value % 1) * 60.

    return deg, orientation, dec_deg

