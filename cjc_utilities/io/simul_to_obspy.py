""" Read SIMUL output to obspy Catalog. """

from typing import Tuple

from obspy import UTCDateTime
from obspy.core.event import (
    Event, Pick, Catalog, ResourceIdentifier, WaveformStreamID, Origin,
    Arrival, OriginQuality, Magnitude)
from obspy.geodetics import kilometers2degrees


def read_simul(phase_file: str, location_file: str) -> Catalog:
    """ 
    Read SIMUL output to Obspy.

    Parameters
    ----------
    phase_file:
        The input phase file for SIMUL
    location_file:
        The output location file from SIMUL
    
    Returns
    -------
    An obspy Catalog of events.
    """
    catalog = Catalog()

    with open(phase_file, "r") as f:
        phase_lines = f.read().split("\n")
    with open(location_file, "r") as f:
        location_lines = f.read().split("\n")
    
    # We want Event -IDs from the input file
    head = location_lines[0:2]
    events = []
    event = []
    for line in location_lines[3:]:
        if line.startswith("  DATE") and len(event) > 0:
            events.append(event)
            event = []
        elif len(line) > 0 and not line.startswith("  DATE"):
            event.append(line)
    
    input_origins = _get_input_origins(
        origin_lines=[line for line in phase_lines if len(line) > 30])

    for event_lines in events:
        event = _convert_event_lines(event_lines)
        # Find closest origin time.
        matched_origins = sorted(
            [(key, abs(event.origins[0].time - value))
             for key, value in input_origins.items()],
            key=lambda e: e[1])
        if matched_origins[0][1] > 5.0:
            raise ValueError("No origin found within 5s of original...")
        event_id = matched_origins[0][0]
        event.resource_id = ResourceIdentifier(event_id)
        catalog += event

    used_origin_ids = {ev.resource_id.id for ev in catalog}
    original_origin_ids = {key for key in input_origins.keys()}
    for origin_id in original_origin_ids.difference(used_origin_ids):
        print("Event origin {0} not matched to a location".format(origin_id))

    return catalog


def _get_input_origins(origin_lines: list) -> dict:
    input_origins = dict()
    for line in origin_lines:
        event_id = line.split()[-1]
        origin_time = UTCDateTime.strptime(
            line[0:11], "%y%m%d %H%M") + float(line[12:17])
        input_origins.update({event_id: origin_time})

    return input_origins


def _convert_event_lines(lines: list) -> Event:
    event = Event()

    origin_time = UTCDateTime(
        year=2000 + int(lines[0][1:3]), month=int(lines[0][3:5]),
        day=int(lines[0][5:7]), hour=int(lines[0][8:10]),
        minute=int(lines[0][10:12])) + float(lines[0][13:18])
    latitude= float(lines[0][19:21]) + float(lines[0][22:26]) / 60
    if lines[0][21] == "S":
        latitude *= -1
    longitude = float(lines[0][28:31]) + float(lines[0][32:37]) / 60
    if lines[0][31] == "W":
        longitude *= -1
    depth = float(lines[0][37:44]) * 1000.0  # Obspy depths are in m
    magnitude = float(lines[0][46:51])
    nobs = int(lines[0][51:54])
    rms = float(lines[0][62:])
    origin = Origin(
        time=origin_time, latitude=latitude, longitude=longitude,
        depth=depth, method_id=ResourceIdentifier("SIMUL"),
        origin_quality=OriginQuality(
            used_phase_count=nobs, standard_error=rms),
        arrivals=[])
    magnitude = Magnitude(mag=magnitude)

    picks = []
    arrivals = []

    p_lines = {line.split()[0]: line for line in lines[2:] if line[23] == "P"}
    sp_lines = {line.split()[0]: line 
                for line in lines[2:] if line[23:25] == "Sp"}
    for station, line in p_lines.items():
        pick, arrival = _extract_pick_info(line=line, origin_time=origin_time)
        picks.append(pick)
        arrivals.append(arrival)
        sp_line = sp_lines.get(station, None)
        if sp_line:
            pick, arrival = _extract_pick_info(
                line=sp_line, origin_time=origin_time, p_time=pick.time)
            picks.append(pick)
            arrivals.append(arrival)

    # Associate everything with the event!
    origin.arrivals = arrivals
    event.origins = [origin]
    event.magnitudes = [magnitude]
    event.picks = picks
    return event


def _extract_pick_info(
    line: str,
    origin_time: UTCDateTime,
    p_time: UTCDateTime = None
) -> Tuple[Pick, Arrival]:
    station = line[1:6]
    distance = kilometers2degrees(float(line[6:13]))
    azimuth = float(line[14:17])
    phase = line[23]
    if phase == "S" and p_time is None:
        raise NotImplementedError("S-P phase given, but no P-time")
    weight = int(line[25])
    if phase == "P":
        toa = int(line[18:21])
        hour = int(line[27:29])
        minute = int(line[29:31])
        seconds = float(line[33:38])
        pick_time = UTCDateTime(
            year=origin_time.year, month=origin_time.month, 
            day=origin_time.day, hour=hour, minute=minute) + seconds
    else:
        toa = None
        pick_time = p_time + float(line[40:45])
    time_residual = float(line[53:59])
    pwt = float(line[60:])  # I don't know what this is?

    pick = Pick(
        time=pick_time, waveform_id=WaveformStreamID(station_code=station),
        phase_hint=phase)
    arrival = Arrival(
        pick_id=pick.resource_id, phase=phase, time_residual=time_residual,
        distance=distance, takeoff_angle=toa, azimuth=azimuth)

    return (pick, arrival)


def print_line_numbers(line):
    """ Unused convenience func for fixed format working. """
    numbers = "".join([str(i % 10) for i in range(len(line))])
    print("\n".join([line, numbers]))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="SIMUL conversion tool")
    parser.add_argument(
        "-p", "--phase-file", help="Input phase file for SIMUL", required=True)
    parser.add_argument(
        "-l", "--location-file", help="SIMUL output file", required=True)
    parser.add_argument(
        "-o", "--output", help="Output for QuakeML writing", required=True)
    
    args = vars(parser.parse_args())

    cat = read_simul(
        phase_file=args["phase_file"], 
        location_file=args["location_file"])
    cat.write(args["output"], format="QUAKEML")