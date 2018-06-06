from obspy.imaging.maps import *
from obspy.imaging.maps import _plot_basemap_into_axes
from obspy.imaging.cm import obspy_sequential
from obspy import UTCDateTime
import numpy as np
import datetime
from matplotlib.dates import AutoDateFormatter, AutoDateLocator, date2num


def test_animated_plot(lons, lats, times, size, color, step, decay, projection,
                       resolution, continent_fill_color, water_fill_color,
                       colormap=None, colorbar=None,
                       colorbar_ticklabel_format=None, marker="o", title=None,
                       interval=200, show=True, **kwargs):
    """
    Make an animated plot of locations and associated times.

    :type lons: list/tuple of floats
    :param lons: Longitudes of the data points.
    :type lats: list/tuple of floats
    :param lats: Latitudes of the data points.
    :type times: list/tuple of UTCDateTime
    :param times: Origin times of data points.
    :type size: float or list/tuple of floats
    :param size: Size of the individual points in the scatter plot.
    :type colors: list/tuple of floats (or objects that can be
        converted to floats, like e.g.
        :class:`~obspy.core.utcdatetime.UTCDateTime`)
    :param colors: Color information of the individual data points to be
        used in the specified color map (e.g. origin depths,
        origin times).
    :type labels: list/tuple of str
    :param labels: Annotations for the individual data points.
    :type step: int
    :param step: Time in seconds for grouping and plotting earthquakes,
        will plot all earthquakes within this step size together.  This is
        also the frame-rate, e.g. one frame per step.
    :type decay: int
    :param decay: Time in seconds to leave locations plotted.  Plotted
            locations will fade (reduce opacity) linearly through the decay.
    :type projection: str
    :param projection: The map projection.
        Currently supported are:

            * ``"global"`` (Will plot the whole world.)
            * ``"ortho"`` (Will center around the mean lat/long.)
            * ``"local"`` (Will plot around local events)
    :type resolution: str
    :param resolution: Resolution of the boundary database to use. Will be
        based directly to the basemap module. Possible values are:

            * ``"c"`` (crude)
            * ``"l"`` (low)
            * ``"i"`` (intermediate)
            * ``"h"`` (high)
            * ``"f"`` (full)

    :type continent_fill_color: Valid matplotlib color
    :param continent_fill_color:  Color of the continents.
    :type water_fill_color: Valid matplotlib color
    :param water_fill_color: Color of all water bodies.
    :type colormap: str, any matplotlib colormap
    :param colormap: The colormap for color-coding the events as provided
        in `color` kwarg.
        The event with the smallest `color` property will have the
        color of one end of the colormap and the event with the highest
        `color` property the color of the other end with all other events
        in between.
    :type marker: str
    :param marker: Any valid matplotlib marker style descriptor.
    :type title: str
    :param title: Title above plot.
    :type interval: int
    :param interval: Interval between frames in ms.
    :type show: bool
    :param show: Whether to show the figure after plotting or not. Can be used
        to do further customization of the plot before showing it.
    """
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    if not HAS_BASEMAP:
        raise ImportError('Basemap cannot be imported but was implicitly '
                          'requested.')
    if any([isinstance(c, (datetime.datetime, UTCDateTime)) for c in color]):
        datetimeplot = True
        color = [
            (np.isfinite(float(t)) and
             date2num(getattr(t, 'datetime', t)) or
             np.nan)
            for t in color]
    else:
        datetimeplot = False
    fig = plt.figure()
    if projection == "local":
        ax_x0, ax_width = 0.10, 0.80
    elif projection == "global":
        ax_x0, ax_width = 0.01, 0.98
    else:
        ax_x0, ax_width = 0.05, 0.90

    if colorbar:
        map_ax = fig.add_axes([ax_x0, 0.13, ax_width, 0.77])
        cm_ax = fig.add_axes([ax_x0, 0.05, ax_width, 0.05])
    else:
        ax_y0, ax_height = 0.05, 0.85
        if projection == "local":
            ax_y0 += 0.05
            ax_height -= 0.05
        map_ax = fig.add_axes([ax_x0, ax_y0, ax_width, ax_height])

    bmap = _plot_basemap_into_axes(ax=map_ax, lons=lons, lats=lats, size=size,
                                   color=color, projection=projection,
                                   resolution=resolution,
                                   continent_fill_color=continent_fill_color,
                                   water_fill_color=water_fill_color,
                                   title=title, animate=True)
    groups = []
    min_time = min(times)
    n_groups = int(((max(times) - min(times)) / step) + (decay / step))
    alphas = np.linspace(1, 0, n_groups)
    min_color = min(color)
    max_color = max(color)

    scal_map = ScalarMappable(norm=Normalize(min_color, max_color),
                              cmap=colormap)
    ev_info = [[lat, lon, c, s, t]
               for lat, lon, c, s, t in zip(lats, lons, color, size, times)]
    ev_info.sort(key=lambda tup: tup[-1])
    for i in range(n_groups):
        group = []
        for ev in ev_info:
            for j, alpha in enumerate(alphas):
                if (i + j) * step < ev[4] - min_time < (i + j + 1) * step:
                    ev.append(alpha)
                    group.append(list(ev))
                else:
                    continue
        groups.append(group)
    groups.insert(0, [])
    x = [g[1] for g in groups[0]]
    y = [g[0] for g in groups[0]]
    s = [g[3] for g in groups[0]]
    c = [g[2] for g in groups[0]]
    a = [g[5] for g in groups[0]]
    x, y = bmap(x, y)
    points = bmap.scatter(x, y, s=s, lw=0.5, marker=marker, c=c,
                          alpha=a[0], cmap=colormap, zorder=10)
    print(points)
    if colorbar:
        if colorbar_ticklabel_format is not None:
            if isinstance(colorbar_ticklabel_format, (str, native_str)):
                formatter = FormatStrFormatter(colorbar_ticklabel_format)
            elif hasattr(colorbar_ticklabel_format, '__call__'):
                formatter = FuncFormatter(colorbar_ticklabel_format)
            elif isinstance(colorbar_ticklabel_format, Formatter):
                formatter = colorbar_ticklabel_format
            locator = MaxNLocator(5)
        else:
            if datetimeplot:
                locator = AutoDateLocator()
                formatter = AutoDateFormatter(locator)
                # Compat with old matplotlib versions.
                if hasattr(formatter, "scaled"):
                    formatter.scaled[1 / (24. * 60.)] = '%H:%M:%S'
            else:
                locator = None
                formatter = None

        # normal case: axes for colorbar was set up in this routine
        if "cm_ax" in locals():
            cb_kwargs = {"cax": cm_ax}
        # unusual case: reusing a plot that has no colorbar set up previously
        else:
            cb_kwargs = {"ax": map_ax}
        cb = fig.colorbar(
            mappable=points, cmap=colormap, orientation='horizontal',
            ticks=locator, format=formatter, **cb_kwargs)

    def animate(_i):
        xy = [bmap(g[1], g[0]) for g in groups[_i]]
        s = [g[3] for g in groups[_i]]
        c = scal_map.to_rgba([g[2] for g in groups[_i]])
        a = [g[5] for g in groups[_i]]
        if len(a) == 0:
            a = 0
        else:
            a = a[0]
        points.set_facecolors(c)
        points.set_sizes(s)
        points.set_offsets(xy)
        points.set_alpha(a)
        return points

    anim = animation.FuncAnimation(plt.gcf(), animate, frames=n_groups,
                                   interval=interval)

    if show:
        plt.show()
    else:
        return anim


if __name__ == '__main__':
    lons = [173.8327766, 174.3092149, 173.4703075, 172.914321, 173.6952679,
            173.7258489, 173.8848195, 174.0049249, 174.2375215, 173.8310381,
            173.2884741, 173.4638208, 174.2928826, 174.2526478, 173.5961806,
            173.9787327, 173.8085293, 173.5475915, 174.0427485, 174.3182552]
    lats = [-42.38910312, -41.83240977, -42.39025399, -42.64619075,
            -42.07473242, -42.2358171, -42.29550808, -42.32569712, -41.6901737,
            -42.24054879, -42.47008968, -42.39173439, -41.69179934,
            -41.71313951, -42.46233008, -42.0095946, -42.26492239,
            -42.53100965, -41.87465439, -41.72924518]
    times = [UTCDateTime('2016-11-14T23:52:25.306668Z'),
             UTCDateTime('2016-11-14T23:47:19.725126Z'),
             UTCDateTime('2016-11-14T23:45:57.964739Z'),
             UTCDateTime('2016-11-14T23:44:52.573135Z'),
             UTCDateTime('2016-11-14T23:41:54.645464Z'),
             UTCDateTime('2016-11-14T23:37:51.389411Z'),
             UTCDateTime('2016-11-14T23:36:48.798458Z'),
             UTCDateTime('2016-11-14T23:36:01.110241Z'),
             UTCDateTime('2016-11-14T23:33:20.548524Z'),
             UTCDateTime('2016-11-14T23:31:01.385975Z'),
             UTCDateTime('2016-11-14T23:30:29.523533Z'),
             UTCDateTime('2016-11-14T23:29:01.861621Z'),
             UTCDateTime('2016-11-14T23:25:39.571550Z'),
             UTCDateTime('2016-11-14T23:20:17.766446Z'),
             UTCDateTime('2016-11-14T23:17:43.296723Z'),
             UTCDateTime('2016-11-14T23:13:40.581251Z'),
             UTCDateTime('2016-11-14T23:11:17.491870Z'),
             UTCDateTime('2016-11-14T23:09:36.462545Z'),
             UTCDateTime('2016-11-14T23:07:29.087792Z'),
             UTCDateTime('2016-11-14T23:04:11.732366Z')]
    colors = [35.46875, 16.71875, 19.0625, 18.59375, 22.34375, 27.03125,
              22.8125, 19.0625, 14.84375, 30.3125, 22.34375, 5.46875,
              5.46875, 15.546875, 15.78125, 24.21875, 25.15625, 7.34375,
              28.90625, 17.65625]
    fig = test_animated_plot(lons=lons, lats=lats, times=times,
                             size=np.arange(200, 220), color=colors, step=600,
                             decay=1200, projection='local', resolution='h',
                             continent_fill_color='0.9',
                             water_fill_color='1.0',
                             colormap=obspy_sequential, colorbar=True,
                             marker="o", title='test animation',
                             interval=200, show=True)
