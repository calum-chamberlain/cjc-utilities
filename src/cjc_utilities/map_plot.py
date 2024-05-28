"""
Simple map plotting for earthquake catalogues

"""

import pygmt
import numpy as np

from obspy.core.event import Catalog, Event


def _get_origin_attr(eq: Event, attr: str):
    try:
        ori = eq.preferred_origin() or eq.origins[-1]
    except IndexError:
        return None
    return ori[attr]


def _get_magnitude_attr(eq: Event, attr: str):
    try:
        mag = eq.preferred_magnitude() or eq.magnitudes[-1]
    except IndexError:
        return None
    return mag[attr]


def map_plot(
    catalog: Catalog,
    mag_fill_value: float = 3.0
) -> pygmt.Figure:

    latitudes = np.array([_get_origin_attr(eq, "latitude") or np.nan for eq in catalog])
    longitudes = np.array([_get_origin_attr(eq, "longitude") or np.nan for eq in catalog])
    depths = np.array([_get_origin_attr(eq, "depth") or np.nan for eq in catalog])
    times = np.array([_get_origin_attr(eq, "time") for eq in catalog])
    magnitudes = np.array([_get_magnitude_attr(eq, "mag") or np.nan for eq in catalog])

    depths /= 1000.0

    magnitudes = np.nan_to_num(magnitudes, mag_fill_value)

    lat_range = latitudes.max() - latitudes.min()
    lon_range = longitudes.max() - longitudes.min()
    region = [
        longitudes.min() - 0.1 * lon_range,
        longitudes.max() + 0.1 * lon_range,
        latitudes.min() - 0.1 * lat_range,
        latitudes.max() + 0.1 * lat_range]

    fig = pygmt.Figure()
    pygmt.makecpt(cmap='lajolla', series=[
              depths.min(), depths.max()])

    pygmt.config(MAP_FRAME_TYPE='plain', FORMAT_GEO_MAP='ddd.xx')
    fig.coast(region=region,
              shorelines=True,
              land='grey',
              water='lightblue',
              projection='M10c',
              frame=['WSne', 'xa2f1', 'ya2f1'])
    fig.plot(x=longitudes,
             y=latitudes,
             size=0.1 * (2 ** magnitudes),
             fill=depths,
             cmap=True,
             style='cc', pen='black')
    fig.colorbar(frame='af+l"Depth (km)""')
    return fig


def main():
    import argparse
    from obspy import read_events

    parser = argparse.ArgumentParser(description="Make a quick map plot of a catalog")

    parser.add_argument("-c", "--catalog", type=str, required=True,
                        help="Catalog to plot")

    args = parser.parse_args()

    fig = map_plot(read_events(args.catalog))
    fig.show()


if __name__ == "__main__":
    main()
