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
    mag_fill_value: float = 3.0,
    scalebymagnitude: bool = True,
    colorby: str = "depth",
    cmap: str = "lajolla",
    transparency: float = 50.,
) -> pygmt.Figure:

    latitudes = np.array(
        [_get_origin_attr(eq, "latitude") or np.nan for eq in catalog])
    longitudes = np.array(
        [_get_origin_attr(eq, "longitude") or np.nan for eq in catalog])
    depths = np.array(
        [_get_origin_attr(eq, "depth") or np.nan for eq in catalog]) / 1000.0
    times = np.array(
        [_get_origin_attr(eq, "time") for eq in catalog])
    if scalebymagnitude:
        magnitudes = np.array(
            [_get_magnitude_attr(eq, "mag") or np.nan for eq in catalog])
    else:
        magnitudes = mag_fill_value * np.ones_like(latitudes)

    magnitudes = np.nan_to_num(magnitudes, mag_fill_value)
    if scalebymagnitude:
        size = 0.1 * (2 ** magnitudes)
    else:
        size = magnitudes

    if colorby is None:
        colors = np.zeros_like(latitudes)
        pen = "0.1p,black"
    elif colorby.lower() == "depth":
        colors = depths
        pen = None
    elif colorby.lower() == "time":
        zerotime = min(times)
        pen = None
        colors = np.array([t - zerotime for t in times])
        if max(colors) - min(colors) < 3600:
            colors /= 60
            time_unit = "Minutes"
        elif max(colors) - min(colors) < 360000:
            colors /= 3600
            time_unit = "Hours"
        else:
            colors /= 86400
            time_unit = "Days"
    else:
        raise Exception(f"Don't know how to colorby {colorby}")

    order = np.argsort(colors)

    lat_range = latitudes.max() - latitudes.min()
    lon_range = longitudes.max() - longitudes.min()
    region = [
        longitudes.min() - 0.1 * lon_range,
        longitudes.max() + 0.1 * lon_range,
        latitudes.min() - 0.1 * lat_range,
        latitudes.max() + 0.1 * lat_range]

    fig = pygmt.Figure()
    if colorby:
        pygmt.makecpt(cmap=cmap, series=[
            colors.min(), colors.max()])

    pygmt.config(MAP_FRAME_TYPE='plain', FORMAT_GEO_MAP='ddd.xx')
    fig.coast(region=region,
              shorelines=True,
              land='grey',
              water='lightblue',
              projection='M10c',
              frame=['WSne', 'xa2f1', 'ya2f1'])
    plot_kwargs = dict(
        x=longitudes[order],
        y=latitudes[order],
        size=size[order],
        style='cc', 
        pen=pen,
        transparency=transparency)
    if colorby:
        plot_kwargs.update(dict(
            fill=colors[order],
            cmap=True))
    else:
        plot_kwargs.update(dict(fill="gray98"))
    fig.plot(**plot_kwargs)
    if colorby and colorby.lower() == "depth":
        fig.colorbar(frame='af+lDepth (km)')
    elif colorby and colorby.lower() == "time":
        fig.colorbar(frame=f"af+l{time_unit} since {zerotime}")
    return fig


def main():
    import argparse
    from obspy import read_events

    parser = argparse.ArgumentParser(
        description="Make a quick map plot of a catalog")

    parser.add_argument("-c", "--catalog", type=str, required=True,
                        help="Catalog to plot")

    args = parser.parse_args()

    fig = map_plot(read_events(args.catalog))
    fig.show()


if __name__ == "__main__":
    main()
