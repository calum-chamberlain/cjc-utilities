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


def write_simul(catalog: Catalog, filename: str, min_stations: int = 5) -> int:
    event_str = [_get_simul_string(event, min_stations=min_stations)
                 for event in catalog]
    event_str = [_ for _ in event_str if _]
    # print(f"Writing {len(event_str)} events to file.")
    event_str = "\n    0\n".join(event_str)
    with open(filename, "w") as f:
        f.write(event_str)
        f.write("\n    0\n")
    return len(event_str)


def _get_simul_string(event: Event, min_stations: int) -> str:
    out = []
    try:
        origin = event.preferred_origin() or event.origins[0]
    except IndexError:
        print("No origin found, skipping")
        return
    try:
        magnitude = event.preferred_magnitude() or event.magnitudes[0]
    except IndexError:
        # print("No magnitude found, setting to 0.0")
        _magnitude = namedtuple("fake_magnitude", ["mag"])
        magnitude = _magnitude(0.0)
    lat, n_s, lat_min = _deg_to_deg_dec_min(origin.latitude, "latitude")
    lon, e_w, lon_min = _deg_to_deg_dec_min(origin.longitude, "longitude")
    depth = origin.depth / 1000.

    origin_time = origin.time
    for pick in event.picks:
        # SIMUL does not cope with negative pick times, but the origin time 
        # doesn't matter that much.
        if origin_time > pick.time:
            origin_time = pick.time - 5

    out.append(HEAD_STR.format(
        year=origin_time.year % 100, month=origin_time.month, 
        day=origin_time.day, hour=origin_time.hour, 
        minute=origin_time.minute,
        seconds=origin_time.second + (origin_time.microsecond / 1e6),
        lat=lat, N_S=n_s, lat_min=lat_min, lon=lon, E_W=e_w, 
        lon_min=lon_min, depth=depth, magnitude=magnitude.mag,
        nobs=len(event.picks), dunno=100, dunnoTwo=0, 
        event_id=event.resource_id.id.split('/')[-1]))
    
    # Only write out P and S with P - SIMUL uses S-P time, not S-time
    p_picked_stations = {p.waveform_id.station_code for p in event.picks 
                         if p.phase_hint.startswith("P")}
    if len(p_picked_stations) < min_stations or len(p_picked_stations) == 0:
        return None
    picks = [p for p in event.picks if p.waveform_id.station_code in p_picked_stations]
    for pick in picks:
        if not pick.phase_hint.startswith(("P", "S")):
            # SIMUL will silently crash on other types of pick :(
            continue
        out.append(PICK_STR.format(
            station=pick.waveform_id.station_code,
            phase=pick.phase_hint,
            weight=1,
            travel_time=pick.time - origin_time))

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

