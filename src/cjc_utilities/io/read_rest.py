"""
Simple functions for reading REST formatted files.

"""

from obspy.core.event import (
    Event, Catalog, Pick, Origin, Magnitude, QuantityError,
    ResourceIdentifier, OriginQuality, WaveformStreamID,
    Arrival, Amplitude)
from obspy import UTCDateTime
from obspy.geodetics import kilometers2degrees


SPLITS = [0, 6, 10, 15, 18, 22, 28, 29, 41, 49, 56, 64, -1]


def rest_to_obspy(filename):
    """
    Read REST formatted event info to an ObsPy Event object.

    :param filename: File to read from
    :type filename: str

    :returns: :class:`obspy.core.event.Event`
    """
    catalog = Catalog()
    with open(filename, 'r') as f:
        full_str = [line for line in f]
    event_str = []
    for line in full_str:
        if len(line.rstrip(" \n\r")) != 0:
            event_str.append(line)
        else:
            event = read_origin(event_str)
            event = read_picks(event_str, event)
            catalog.events.append(event)
            event_str = []
    return catalog


def is_rest(event_str):
    """
    Check if the format is as expected.

    :param event_str: Contents of file as list of str
    :type event_str: list

    :returns: bool
    """
    if len(event_str[0].rstrip()) is not 141:
        return False
    if event_str[0][0:4] != 'YEAR':
        return False
    if event_str[3][0:3] != 'STA':
        return False
    return True


def read_origin(event_str):
    """
    Read the origin information from the REST file string

    :param event_str: Contents of file as list of str
    :type event_str: list

    :returns: :class:`obspy.core.event.Event`
    """
    head = event_str[0].split()
    event_id = head[8]
    if event_id == 'NOEVID':
        try:
            event_id = head[20]
        except IndexError:
            event_id = None
    if event_id:
        event = Event(resource_id=ResourceIdentifier(id=event_id))
    else:
        event = Event()
    try:
        gap = float(head[17])
    except IndexError:
        gap = None
    origin = Origin(
        time=UTCDateTime(
            year=int(head[0]), julday=int(head[1]), hour=int(head[2]),
            minute=int(head[3])) + float(head[4]),
        latitude=float(head[5]), longitude=float(head[6]),
        depth=float(head[7]) * 1000, origin_quality=OriginQuality(
            standard_error=float(head[9]),
            azimuthal_gap=gap,
            used_phase_count=int(head[17])),
        longitude_errors=QuantityError(
            uncertainty=kilometers2degrees(float(head[12]))),
        latitude_errors=QuantityError(
            uncertainty=kilometers2degrees(float(head[11]))),
        depth_errors=QuantityError(uncertainty=float(head[13]) * 1000),
        method_id=ResourceIdentifier("smi:local/REST"),
        evaluation_mode="automatic")
    event.origins.append(origin)
    try:
        event.magnitudes.append(Magnitude(
            mag=float(head[18]), magnitude_type="M"))
    except IndexError:
        pass
    return event


def read_picks(event_str, event):
    """
    Read the picks from the REST file string

    :param event_str: Contents of file as list of str
    :type event_str: list
    :param event:
        Event to assoaite the picks with. Note old picks will
        not be overwritten. Event should have only one origin.
    :type event: :class:`obspy.core.event.Event`

    :returns: :class:`obspy.core.event.Event`
    """
    for line in event_str[1:]:
        pick, arrival, amplitude = read_pick(line)
        event.picks.append(pick)
        event.origins[0].arrivals.append(arrival)
        if amplitude:
            event.amplitudes.append(amplitude)
    return event


def read_pick(line):
    """
    Convert REST pick string to ObsPy Pick object

    :param line: string containing pick information
    :type line: str

    :returns:
        :class: `obspy.core.event.Pick` and
        :class: `obspy.core.event.origin.Arrival`
        :class: `obspy.core.event.origin.Amplitude`
    """
    # line = line.split()  # Cannot just split the line :(
    amplitude = None
    _line = []
    for split in range(len(SPLITS) - 1):
        _line.append(line[SPLITS[split]: SPLITS[split + 1]].strip())
    line = _line
    if '*' in line[6] or 'X' in line[6]:
        pick = Pick(time=UTCDateTime(
            year=int(line[1]), julday=int(line[2]), hour=int(line[3]),
            minute=int(line[4])) + float(line[5]), phase_hint=line[7],
            evaluation_mode="automatic",
            method_id=ResourceIdentifier("smi:local/REST"),
            waveform_id=WaveformStreamID(station_code=line[0]),
            time_errors=QuantityError(uncertainty=float(line[8]), confidence_level=0))
    else:
        pick = Pick(time=UTCDateTime(
            year=int(line[1]), julday=int(line[2]), hour=int(line[3]),
            minute=int(line[4])) + float(line[5]), phase_hint=line[7],
            evaluation_mode="automatic",
            method_id=ResourceIdentifier("smi:local/REST"),
            waveform_id=WaveformStreamID(station_code=line[0]),
            time_errors=QuantityError(uncertainty=float(line[8])))
    onset = float(line[-3])
    if onset == 0:
        pick.polarity = "undecidable"
    else:
        amplitude = Amplitude(
            waveform_id=pick.waveform_id,  pick_id=pick.resource_id,
            generic_amplitude=onset, type='A', category='point')
        if onset < 0:
            pick.polarity = "negative"
        elif onset > 0:
            pick.polarity = "positive"
    arrival = Arrival(
        pick_id=pick.resource_id, time_residual=float(line[9]))
    return pick, arrival, amplitude