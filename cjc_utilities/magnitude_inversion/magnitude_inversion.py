"""
Python port of old Matlab code from Carolin Boese, updated by Emily Warren-Smith
and ported by Calum Chamberlain.

Calibrate amplitude picks from a network against "known" magnitudes for some
events.
"""

import numpy as np
import pandas as pd
import warnings

import matplotlib.pyplot as plt

from typing import Union, Tuple

from progressbar import ProgressBar

from scipy import linalg

from obspy.clients.fdsn import Client
from obspy import read_events, UTCDateTime, Inventory
from obspy.geodetics import gps2dist_azimuth, degrees2kilometers
from obspy.core.event import (
    Magnitude, StationMagnitude, StationMagnitudeContribution, 
    ResourceIdentifier, Catalog, Event, Amplitude, Arrival)


plt.style.use("ggplot")


MIN_CONSTRAINTS = 20  # Minimum number of matched events to constrain the inversion


def _get_origin_attrib(
    event: Event, 
    attribute: str
) -> Union[float, UTCDateTime]:
    """ Get an event origin attribute. """
    try:
        return (event.preferred_origin() or event.origins[-1]
                ).__getattribute__(attribute)
    except (IndexError, AttributeError):
        return None


def _get_magnitude_attrib(
    event: Event, 
    attribute: str, 
    magnitude_type: str = None
) -> Union[float, str]:
    """ Get a magnitude attribute. """
    if magnitude_type:
        magnitude = [magnitude for magnitude in event.magnitudes 
                     if magnitude.magnitude_type == magnitude_type]
        if len(magnitude) == 0:
            return None
        magnitude = magnitude[0]
    else:
        try:
            magnitude = event.preferred_magnitude() or event.magnitudes[-1]
        except IndexError:
            return None
    try:
        return magnitude.__getattribute__(attribute)
    except AttributeError:
        return None


def _get_arrival_for_amplitude(amplitude: Amplitude, event: Event) -> Arrival:
    ori = event.preferred_origin() or event.origins[-1]
    if amplitude.pick_id is None:
        print("Amplitude not matched to pick")
        return None
    pick = amplitude.pick_id.get_referred_object()
    if pick is None:
        print("Amplitude not matched to pick")
        return None
    # Need an arrival on this station to get the distance
    arrival = [arr for arr in ori.arrivals 
                if arr.pick_id.get_referred_object().waveform_id.station_code == pick.waveform_id.station_code
                and arr.distance]
    if len(arrival) == 0:
        print(f"No arrival found for {pick.waveform_id.station_code}, skipping")
        return None
    arrival = arrival[0]
    return arrival


def _get_amplitude_value(amplitude: Amplitude) -> float:
    # Convert to m
    amp = amplitude.generic_amplitude
    if amplitude.type == "IAML":
        # Apply Wood Anderson sensitivity - ISAPEI standard is removed
        amp *= 2080.0
    if amplitude.unit is None:
        warnings.warn(
            "No amplitude unit specified, assuming this is m - "
            "use with caution!!!")
    elif amplitude.unit != "m":
        raise NotImplementedError(
            "Only written to work with SI displacement units.")
    return amp


def get_event(
    catalog: Catalog, 
    event_id: Union[str, ResourceIdentifier]
) -> Event:
    """
    Get an event from a catalog given an event id.
    """
    if isinstance(event_id, ResourceIdentifier):
        filtered = [ev for ev in catalog if ev.resource_id == event_id]
    elif isinstance(event_id, str):
        filtered = [ev for ev in catalog if event_id in ev.resource_id.id]
    else:
        print("event id must be either string or ResourceIdentifier")
        return None
    if len(filtered) == 0:
        print("No match")
        return None
    elif len(filtered) > 1:
        print("Multiple matches")
        return filtered
    return filtered[0]


def summarize_catalog(catalog: Catalog, magnitude_type: str = None) -> pd.DataFrame:
    """
    Summarize a catalog into a sparse dataframe of origin information
    """
    if magnitude_type:
        events = [ev for ev in catalog.events 
                  if magnitude_type in [mag.magnitude_type 
                                        for mag in ev.magnitudes]]
    else:
        events = catalog.events
    event_ids = [ev.resource_id.id.split('/')[-1] for ev in events]
    origin_times = [_get_origin_attrib(ev, "time") for ev in events]
    latitudes = [_get_origin_attrib(ev, "latitude") for ev in events]
    longitudes = [_get_origin_attrib(ev, "longitude") for ev in events]
    depths = [_get_origin_attrib(ev, "depth") for ev in events]
    magnitudes = [
        _get_magnitude_attrib(ev, "mag", magnitude_type=magnitude_type)
        for ev in events]
    magnitude_types = [
        _get_magnitude_attrib(ev, "magnitude_type", 
                              magnitude_type=magnitude_type) 
        for ev in events]
    
    catalog_df = pd.DataFrame(data=dict(
        event_id=event_ids, origin_time=origin_times, latitude=latitudes,
        longitude=longitudes, depth=depths, magnitude=magnitudes,
        magnitude_type=magnitude_types))
    return catalog_df


def get_epicentral_distance(event: Event, inventory: Inventory, station: str) -> float:
    origin_time = _get_origin_attrib(event, "time")
    _station = inventory.select(station=station, starttime=origin_time - 10, 
                                endtime=origin_time + 10)
    if len(_station) > 1:
        print("Multiple stations found")
    try:
        _station = _station[0][0]
    except IndexError:
        print("No matching station")
        return None

    dist, _, _ = gps2dist_azimuth(
        lat1=_get_origin_attrib(event, "latitude"),
        lon1=_get_origin_attrib(event, "longitude"),
        lat2=_station.latitude, 
        lon2=_station.longitude)
    return dist / 1000.


def summarize_amplitudes(
    catalog: Catalog, 
) -> pd.DataFrame:
    """
    """
    total_observations = len([_ for ev in catalog for _ in ev.amplitudes 
                              if _.type in ["IAML", "AML", "ML"]])
    out = pd.DataFrame(
        index=np.arange(total_observations), 
        columns=("event_id", "depth", # "origin_time", "latitude", "longitude",
                 "epicentral_distance", "seed_id", "amplitude", 
                 "period", "seisan_magnitude")) # , "namps"))
    
    i = 0
    bar = ProgressBar(max_value=total_observations)
    for ev in catalog:
        ori = ev.preferred_origin() or ev.origins[-1]
        namps = len([amp for amp in ev.amplitudes 
                     if amp.type in ["IAML", "AML", "ML"]])
        rid, ori_time, depth, lat, lon = (
            ev.resource_id.id.split('/')[-1],
            _get_origin_attrib(ev, "time"),
            _get_origin_attrib(ev, "depth") / 1000.0,
            _get_origin_attrib(ev, "latitude"),
            _get_origin_attrib(ev, "longitude"))
        for amplitude in ev.amplitudes:
            if amplitude.type not in  ["IAML", "AML", "ML"]:
                continue         
            arrival = _get_arrival_for_amplitude(amplitude, ev)
            if arrival is None:
                continue
            # Convert to m
            amp = _get_amplitude_value(amplitude)

            epi_dist = degrees2kilometers(arrival.distance)

            out.event_id[i] = rid
            # out.origin_time[i] = ori_time
            out.depth[i] = depth
            # out.latitude[i] = lat
            # out.longitude[i] = lon 
            out.epicentral_distance[i] = epi_dist
            out.seed_id[i] = amplitude.waveform_id.get_seed_string()
            out.amplitude[i] = amp
            out.period[i] = amplitude.period
            out.seisan_magnitude[i] = (
                np.log10(amp * 1000) + np.log10(epi_dist) + 
                (0.0067 * 0.4343 * amp * 1000))
            # out.namps[i] = namps
            i += 1
            bar.update(i)
    bar.finish()
    # Set dtypes
    floating_columns = (
        "depth", "epicentral_distance", "amplitude", "period", 
        "seisan_magnitude")

    for column in floating_columns:
        out[column] = out[column].astype(np.float32)
    return out.dropna(how="all")
            

def find_matching_events(
    catalog_1: Catalog, 
    catalog_2: Catalog, 
    time_difference: float = 5.0,
    epicentral_difference: float = 20.0,
    depth_difference: float = 40.0,
    magnitude_type: str = None,
) -> dict:
    """
    Find matching events between two catalogs.

    Parameters
    ----------
    catalog_1
        A catalog to compare to catalog_2
    catalog_2
        A catalog to compare to catalog_1
    time_difference
        Maximum allowed difference in origin time in seconds between matching 
        events
    epicentral_difference
        Maximum allowed difference in epicentral distance in km between 
        matching events
    depth_difference
        Maximum allowed difference in depth in km between matching events.
    magnitude_type
        Magnitude type for comparison, will only return events in catalog_1
        with this magnitude

    Returns
    -------
    Dictionary of matching events ids. Keys will be from catalog_1.
    """
    df_1 = summarize_catalog(catalog_1, magnitude_type=magnitude_type)
    df_2 = summarize_catalog(catalog_2)
    if len(df_2) == 0:
        return None

    swapped = False
    if len(df_1) > len(df_2):
        # Flip for efficiency, will loop over the shorter of the two and use 
        # more efficient vectorized methods on the longer one
        df_1, df_2 = df_2, df_1
        swapped = True

    timestamp = min(min(df_1.origin_time), min(df_2.origin_time))
    comparison_times = np.array([t - timestamp for t in df_2.origin_time])

    print("Starting event comparison.")
    bar = ProgressBar(max_value=len(df_1))
    matched_quakes = dict()
    for i in range(len(df_1)):
        origin_time = UTCDateTime(df_1.origin_time[i])
        origin_seconds = origin_time - timestamp
        deltas = np.abs(comparison_times - origin_seconds)
        index = np.argmin(deltas)
        delta = deltas[index]
        if delta > time_difference:
            continue  # Too far away in time
        depth_sep = abs(df_1.depth[i] - df_2.depth[index]) / 1000.0  # to km
        if depth_sep > depth_difference:
            continue  # Too far away in depth
        # distance check
        dist, _, _ = gps2dist_azimuth(
            lat1=df_2.latitude[index],
            lon1=df_2.longitude[index],
            lat2=df_1.latitude[i],
            lon2=df_1.longitude[i])
        dist /= 1000.
        if dist > epicentral_difference:
            continue  # Too far away epicentrally
        matched_id = df_2.event_id[index]
        if matched_id in matched_quakes.keys():
            # Check whether this is a better match
            if delta > matched_quakes[matched_id]["delta"] or\
                 dist > matched_quakes[matched_id]["dist"] or\
                 depth_sep > matched_quakes[matched_id]["depth_sep"]:
                continue  # Already matched to another, better matched event
        matched_quakes.update(
            {matched_id: dict(
                delta=delta, dist=dist, depth_sep=depth_sep, 
                matched_id=df_1.event_id[i])})
        bar.update(i + 1)
    bar.finish()
    
    # We just want the event mapper
    if not swapped:
        return {key: value["matched_id"] 
                for key, value in matched_quakes.items()}
    else:
        return {value["matched_id"]: key 
                for key, value in matched_quakes.items()}


def get_comparison_catalog(
    catalog: Catalog, 
    client: Union[Client, str],
    region_tolerance: float = 0.1,
) -> Catalog:
    """
    Get a catalog spanning the same region as the given catalog.
    """
    if isinstance(client, str):
        client = Client(client)
    assert isinstance(client, Client)

    starttime = min(_get_origin_attrib(ev, "time") for ev in catalog) - 3600
    endtime = max(_get_origin_attrib(ev, "time") for ev in catalog) + 3600
    min_lat = min(_get_origin_attrib(ev, "latitude") 
                  for ev in catalog) - region_tolerance
    max_lat = max(_get_origin_attrib(ev, "latitude") 
                  for ev in catalog) + region_tolerance
    min_lon = min(_get_origin_attrib(ev, "longitude") 
                  for ev in catalog) - region_tolerance
    max_lon = max(_get_origin_attrib(ev, "longitude") 
                  for ev in catalog) + region_tolerance

    if (endtime - starttime) / (24 * 3600) > 30:
        comparison_cat = Catalog()
        _starttime, _endtime = starttime, starttime + (30 * 34 * 3600)
        while _endtime <= endtime:
            comparison_cat += client.get_events(
                starttime=_starttime, endtime=_endtime, minlatitude=min_lat, 
                maxlatitude=max_lat, minlongitude=min_lon, maxlongitude=max_lon)
            _starttime += (30 * 34 * 3600)
            _endtime += (30 * 34 * 3600)
        comparison_cat += client.get_events(
            starttime=_starttime, endtime=endtime, minlatitude=min_lat, 
            maxlatitude=max_lat, minlongitude=min_lon, maxlongitude=max_lon)
    else:
        comparison_cat = client.get_events(
            starttime=starttime, endtime=endtime, minlatitude=min_lat, 
            maxlatitude=max_lat, minlongitude=min_lon, maxlongitude=max_lon)
    return comparison_cat


def insert_magnitude(
    event: Event,
    magnitude: float,
    gamma: float,
    frequency_dependent: bool,
    station_corrections: dict, 
    geometric_parameter: float,
    use_mean_correction: bool = False,
) -> Event:
    """
    """
    event_out = event.copy()
    ori = event_out.preferred_origin() or event_out.origins[-1]
    mag_type = "ML"
    mean_sta_corr = np.mean([corr for corr in station_corrections.values()])

    station_magnitudes = []
    for amplitude in event_out.amplitudes:
        if amplitude.type not in  ["IAML", "AML", "ML"]:
            continue
        arrival = _get_arrival_for_amplitude(amplitude, event)
        if arrival is None:
            continue       
        amp = _get_amplitude_value(amplitude)
        hyp_dist = (
            degrees2kilometers(arrival.distance) ** 2 +
            (ori.depth / 1000) ** 2) ** .5  
        # Implictly assuming that stations are at the surface

        if frequency_dependent:
            gamma_freq = gamma / amplitude.period
        else:
            gamma_freq = gamma

        seed_id = amplitude.pick_id.get_referred_object(
            ).waveform_id.get_seed_string()

        sta_corr = station_corrections.get(seed_id, None)
        if sta_corr is None and not use_mean_correction:
            print(f"No station correction for {seed_id}, skipping")
            continue
        elif sta_corr is None:
            sta_corr = mean_sta_corr
            print(f"No station correction for {seed_id}, using mean: {sta_corr}")
        station_mag = (
                np.log10(amp * 1000) + 
                geometric_parameter * np.log10(hyp_dist) + 
                (gamma_freq * 0.4343 * hyp_dist) +  # This was wrong!
                sta_corr)

        station_magnitude = StationMagnitude(
            origin_id=ori.resource_id, mag=station_mag, 
            amplitude_id=amplitude.resource_id, 
            station_magnitude_type=mag_type,
            method_id=ResourceIdentifier("smi:local/magnitude_inversion_CJC"))
        station_magnitudes.append(station_magnitude)
    
    if len(station_magnitudes) == 0:
        print("No suitable amplitudes")
        return event_out
    
    if magnitude is None:
        magnitude = sum(
            sta_mag.mag for sta_mag in station_magnitudes) / len(station_magnitudes)
        method = "smi:local/station-mean"
    else:
        method = "smi:local/magnitude_inversion_CJC"
    final_magnitude = Magnitude(
        mag=magnitude, magnitude_type=mag_type, 
        station_count=len(station_magnitudes),
        method_id=ResourceIdentifier(method),
        station_magnitude_contributions=[
            StationMagnitudeContribution(
                sta_mag.resource_id, residual=sta_mag.mag - magnitude)
            for sta_mag in station_magnitudes])

    event_out.station_magnitudes.extend(station_magnitudes)
    event_out.magnitudes.append(final_magnitude)
    event_out.preferred_magnitude_id = final_magnitude.resource_id
    return event_out


def residual_plot(residuals: np.array) -> plt.Figure:
    """ Plot residuals from magnitude inversion. """
    fig, ax = plt.subplots()
    nbins = max(10, len(residuals) // 20)
    n, bins, _ = ax.hist(residuals, bins=nbins)
    ax.set_xlabel("Magnitude residual")
    ax.set_ylabel("Count")
    
    sigma = residuals.std()
    mu = residuals.mean()

    # Add a best fit normal distribution
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu)) ** 2))
    y /= y.max()
    y *= n.max()
    # ax1 = ax.twinx()
    ax.plot(bins, y, '--')
    
    ax.set_title(
        f"Inverted magnitude residuals, mean: {mu:.3f}, standard dev: {sigma:.3f}")
    return fig


def decay_plot(
    mag_data: pd.DataFrame, 
    gamma: float, 
    frequency_dependent: bool
) -> plt.Figure:
    """ 
    Replicate the period, amplitude and residual plot of Boese et al., 2012 (Fig 7). 
    """
    fig, axes = plt.subplots(nrows=3, sharex=True)

    period_ax = axes[0]
    distance_corr_ax = axes[1]
    residual_ax = axes[2]

    if "hypocentral_distance" not in mag_data.columns:
        assert "depth" in mag_data.columns, "Needs a depth column"
        assert "epicentral_distance" in mag_data.columns, "Needs an epicentral_distance column"
        hyp_dist = (mag_data.depth ** 2 + mag_data.epicentral_distance ** 2) ** .5
        mag_data["hypocentral_distance"] = hyp_dist

    required_columns = ("period", "residual", "hypocentral_distance")
    for col in required_columns:
        assert col in mag_data.columns, f"Needs a column: {col}"

    period_ax.plot(mag_data.hypocentral_distance, mag_data.period, 
                   marker="+", linestyle='None')
    period_ax.set_ylabel("Period at maximum amplitude (s)")

    if frequency_dependent:
        gammaf = gamma / mag_data.period
    else:
        gammaf = gamma
    dist_corr = np.log10(mag_data.hypocentral_distance) + (
        0.4343 * gammaf * mag_data.hypocentral_distance)
    distance_corr_ax.plot(mag_data.hypocentral_distance, dist_corr, 
                          marker="+", label="Calculated", linestyle='None')
    # Plot a line set to 1 Hz
    distances = np.arange(0.1, mag_data.hypocentral_distance.max(), 1.0)
    for freq in (1., 5., 10.):
        fit = np.log10(distances) + 0.4343 * gamma * freq * distances
        distance_corr_ax.plot(distances, fit, label=f"$\gamma({freq} Hz)={gamma * freq:.3f}$")
    distance_corr_ax.set_ylabel("Distance correction term")
    distance_corr_ax.legend()

    residual_ax.plot(mag_data.hypocentral_distance, mag_data.residual, 
                     marker="+", label="Residuals", linestyle='None')
    gradient, intercept = np.polyfit(
        mag_data.hypocentral_distance, mag_data.residual, 1)
    residual_ax.plot(distances, gradient * distances + intercept, 
                     label=f"Fit, ${gradient:.3f}\Delta + {intercept:.3f}$")
    residual_ax.legend()
    residual_ax.set_ylabel("Magnitude residual")
    residual_ax.set_xlabel("Hypocentral distance (km)")

    return fig


def magnitude_inversion(
    new_catalog: Catalog,
    callibration_catalog: Catalog,
    time_difference: float = 5.0,
    epicentral_difference: float = 20.0,
    depth_difference: float = 40.0,
    magnitude_type: str = "Mw",
    geometric_parameter: float = 1.0,
    magnitude_weight: float = 1.0,
    frequency_dependent: bool = True,
    frequency: float = 5.0,
    plot: bool = False,
    only_matched: bool = False,
    in_place: bool = False,
) -> Tuple[Catalog, float, dict]:
    """
    Compute magnitude scale for a given catalogue compared to another.

    Finds the least-squares solution to the attenuation parameter and station
    correction terms for the equation:

    $M_L = log_10(a) + \{alpha}log_10(\{delta}) + 0.4343\{gamma}\{delta} + S$

    Where a is amplitude in mm on a Wood-Anderson instrument, $\{alpha}$ is
    the geometric parameter for spherical attenuation (usually 1), $\{gamma} is
    the anelastic attenuation parameter, $\{delta}$ is the hypocentral distance
    and S is the station correction term.


    Parameters
    ----------
    new_catalog:
        Catalog of events to compute magnitudes for - catalog must have 
        locations and amplitude readings.
    callibration_catalog:
        Catalog to callibrate `new_catalog` against - needs to have magnitudes
        and locations.
    time_difference:
        Origin-time difference in seconds to match events together.
    epicentral_difference:
        Epicentral difference in km to match events together.
    depth_difference:
        Depth difference in km to match events together.
    magnitude_type:
        Magnitude type to callibrate against.
    geometric_parameter:
        Geometric spreading parameter ($\{alpha}$) in the above equation.
    magnitude_weight:
        Weight to apply to callibration magnitudes in the inversion.
    frequency_dependent:
        Whether to compute a frequency dependent $\{gamma} or not.
    frequency:
        Frequency to use for frequency independent inversion (only used if
        `frequency_dependent == False`)
    plot:
        Whether to make plots or not. Plots are shown to screen and not saved.
    only_matched:
        Whether to only use matched-events. Set to True to compute parameters
        for magnitude calculation for later use.
    in_place:
        Update magnitudes in-place in catalog (True), or return a new catalog.

    Returns
    -------
    catalog:
        Output catalog including magnitudes as calculated.
    gamma:
        anelastic attenuation term
    station_corrections:
        Dictionary of station corrections keyed by seed-id.
    """
    print("##########  Finding matching events ##########")

    comparison_events = find_matching_events(
        catalog_1=callibration_catalog, catalog_2=new_catalog,
        time_difference=time_difference, 
        epicentral_difference=epicentral_difference, 
        depth_difference=depth_difference, magnitude_type=magnitude_type)
    
    if len(comparison_events) < MIN_CONSTRAINTS:
        raise NotImplementedError(f"Only {len(comparison_events)} match")
    
    # Use a list here to maintain order.
    matched_ids = list(comparison_events.keys())
    callibration_ids = [comparison_events[matched_id] 
                        for matched_id in matched_ids]
    callibration_events = [
        get_event(callibration_catalog, callibration_id)
        for callibration_id in callibration_ids]
    callibration_magnitudes = np.array([
        _get_magnitude_attrib(event=_ev, attribute="mag",
                              magnitude_type=magnitude_type) 
        for _ev in callibration_events])

    if plot:
        Catalog(callibration_events).plot(projection="local", resolution="h")

    if only_matched:
        matched_events = Catalog(
            [ev for ev in new_catalog 
            if ev.resource_id.id.split('/')[-1] in matched_ids])
        mag_data = summarize_amplitudes(matched_events)
    else:
        mag_data = summarize_amplitudes(new_catalog)

    used_event_ids = list(set(mag_data.event_id))
    used_seed_ids = list(set(mag_data.seed_id))

    # Commonly used lengths
    n_new_events = len(used_event_ids)
    n_observations = len(mag_data)
    n_constraining_events = len(comparison_events)
    n_seed_ids = len(used_seed_ids)
    total_length = n_observations + n_constraining_events

    ########################## Construct Y vector ##############################
    print("########## Constructing matrices ##########")
    Y = np.zeros(total_length, dtype=np.float64)
    # Construct Y
    slantdist = (mag_data.epicentral_distance.to_numpy() ** 2 +
                 mag_data.depth.to_numpy() ** 2) ** 0.5
    slantdist = slantdist.astype(np.float32)
    # Fill top of Y with amplitude and geometric spreading terms
    # Amplitudes are in m
    amplitudes = mag_data.amplitude.to_numpy().astype(np.float32)
    Y[0:n_observations] = np.log10(amplitudes * 1000) + (
        geometric_parameter * np.log10(slantdist))
    # Fill remainder with weighted geonet magnitudes
    Y[n_observations:] = callibration_magnitudes * magnitude_weight

    ########################## Construct X matrix ##############################

    x_width = n_new_events + n_seed_ids + 1
    X = np.zeros((total_length, x_width), dtype=np.float64)

    # Fill with frequency * distance
    if frequency_dependent:
        X[0:n_observations, n_new_events] = (1 / mag_data.period) * slantdist
    else:
        X[0:n_observations, n_new_events] = frequency * slantdist

    # Fill with ones where an amplitude measurement matches an event
    for i, used_event_id in enumerate(used_event_ids):    
        X[0:n_observations, i] = mag_data.event_id == used_event_id


    # Fill with ones where a seed id is used for that event.
    for i, seed_id in enumerate(used_seed_ids):
        X[0:n_observations, n_new_events + i + 1] = mag_data.seed_id == seed_id

    # Fill with ones where events are matched
    for j, callibration_id in enumerate(callibration_ids):
        for k, used_event_id in enumerate(used_event_ids):
            if used_event_id in comparison_events.keys() and \
                comparison_events[used_event_id] == callibration_id:
                 X[n_observations + j, k] = magnitude_weight

    ################################# Solve ####################################
    print("########## Solving inversion ##########")
    conditionX = X.T @ X
    conditionY = X.T @ Y.T
    l, u = linalg.lu(conditionX, permute_l=True)

    params = linalg.lstsq(u, linalg.lstsq(l, conditionY)[0])[0]
    # Least squares solution to matlabs l \ (u \ conditionY)

    _id = np.identity(x_width)    
    condition_inv = linalg.lstsq(u, linalg.lstsq(l, _id)[0])[0] 
    # This is the inverse for the variances
   
    assert(np.allclose(condition_inv @ conditionX, _id, atol=1e-5))
    assert(np.allclose(
        linalg.lstsq(u, linalg.lstsq(l, conditionX)[0])[0], _id, atol=1e-5))

    print("########## Extracting parameters ##########")
    p_std_multiplier = np.sqrt(condition_inv[n_new_events + 1,
                               n_new_events + 1])
    # multiply this by the std dev residuals to get the std dev for p (etc)
    mag_std_multiplier = np.sqrt(np.diag(condition_inv))

    res = Y.T - X @ params
    new_magnitudes = params[0:n_new_events]
    calibration_magnitude_residuals = res[n_observations:total_length]
    p_param = params[n_new_events]
    gamma = p_param / np.log10(np.exp(1))
    gamma *= -1
    cvalue = params[(n_new_events + 1): x_width]

    Xconstraint = X[n_observations + 1:total_length, 1:n_new_events]

    ####################### Get station corrections ############################
    # These are just the average adjustments
    
    # TODO: The station corrections seem wrong: 
    #   the magnitudes calculated using them are not the same as the 
    #   inverted magnitudes.

    best_mags = np.zeros(n_observations)

    for j, event_id in enumerate(used_event_ids):
        filt = np.where(mag_data.event_id == event_id)
        best_mags[filt] = new_magnitudes[j]

    station_corrections = {
        seed_id: -1 * cvalue[i] for i, seed_id in enumerate(used_seed_ids)}

    magnitude_differences = np.empty((n_seed_ids, n_observations)) * np.nan

    for j, seed_id in enumerate(used_seed_ids):
        filt = np.where(mag_data.seed_id == seed_id)
        magnitude_differences[j, filt] = (
            mag_data.seisan_magnitude.to_numpy()[filt] - best_mags[filt])
    
    average_magnitude_differences = np.nanmean(magnitude_differences, axis=1)

    ########################## Add in magnitudes ###############################
    print("########## Building output catalog. ##########")

    if in_place:
        callibrated_events = new_catalog.events
    else:
        callibrated_events = [None for _ in range(len(new_catalog))]

    if only_matched:
        # Calculate average magnitude from stations.
        bar = ProgressBar(max_value=len(new_catalog))
        for i, event in enumerate(new_catalog):
            callibrated_events[i] = insert_magnitude(
                event=event, magnitude=None, gamma=gamma,
                frequency_dependent=frequency_dependent,
                station_corrections=station_corrections,
                geometric_parameter=geometric_parameter)
            bar.update(i)
        bar.finish()
    else:
        # Use the inverted magnitudes.
        bar = ProgressBar(max_value=len(used_event_ids))
        i = 0
        for event_id, magnitude in zip(used_event_ids, new_magnitudes):
            event = [ev for ev in new_catalog 
                    if ev.resource_id.id.split('/')[-1] == event_id][0]
            callibrated_events[i] = insert_magnitude(
                event=event, magnitude=magnitude, gamma=gamma, 
                frequency_dependent=frequency_dependent,
                station_corrections=station_corrections, 
                geometric_parameter=geometric_parameter)
            i += 1
            bar.update(i)
        bar.finish()

    callibrated_catalog = Catalog(callibrated_events)

    ####################### Make plots #########################################
    
    if plot:
        # Add residuals to mag_data
        res_col = res[0:n_observations]
        mag_data["residual"] = res_col
        # Plot of residuals for each observations
        res_plot = residual_plot(residuals=res_col)
        # Plot of residuals for inidivual callibration magnitudes
        callibration_res_plot = residual_plot(
            residuals=calibration_magnitude_residuals)
        callibration_res_plot.suptitle("Callibration events only")

        distance_decay_plot = decay_plot(mag_data, gamma, frequency_dependent)
        plt.show()


    return  callibrated_catalog, gamma, station_corrections


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description="Invert for magnitudes")

    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help=("Input catalogue in an ObsPy readable file. Must have amplitude "
              "readings"))
    parser.add_argument(
        "-c", "--comparison", type=str, required=False, default=None,
        help="Comparison catalogue in ObsPy readable format")
    parser.add_argument(
        "-C", "--client", type=str, required=False, default="GEONET",
        help="FDSN Client to get a comparison catalogue from")
    parser.add_argument(
        "-m", "--magnitude-type", type=str, required=False, default="ML",
        help="Magnitude type to callibrate against")
    parser.add_argument(
        "-o", "--output", type=str, default="magnitude_inverted.xml",
        help="Output catalogue file, written to QUAKEML format")

    args = parser.parse_args()

    new_catalog = read_events(args.input)
    print(f"Read in {len(new_catalog)} events.")

    if args.comparison:
        callibration_catalog = read_events(args.comparison)
    elif args.client:
        callibration_catalog = get_comparison_catalog(
            catalog=new_catalog, client=args.client)
    else:
        IOError("Requires either comparison or client")
    print(f"These will be compared with {len(callibration_catalog)} events")

    output, gamma, station_corrections = magnitude_inversion(
        new_catalog=new_catalog, callibration_catalog=callibration_catalog,
        magnitude_type=args.magnitude_type)

    output.write(args.output, format="QUAKEML")
    print(f"Written callibrated catalog to {args.output}")

    out_parameters = {
        "gamma": gamma,
        "station_corrections": station_corrections}
    
    print("######## RESULTS #########")
    print(f"gamma:\t{gamma}")
    print("STATION CORRECTIONS")
    for station, correction in station_corrections.items():
        print(f"{station}:\t{correction}")
    print("\n\n")

    with open("magnitude_inversion.json", "w") as f:
        json.dump(out_parameters, f)    
    print("Written gamma and station corrections to magnitude_inversion.json")
