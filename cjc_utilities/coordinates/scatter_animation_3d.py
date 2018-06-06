"""
3D scatter animation for earthquakes. Uses the coordinates module.
"""

from coordinates import Location, Geographic
from extract_cross_section import get_plane

import matplotlib.pyplot as plt
import matplotlib.animation
from mpl_toolkits.mplot3d import Axes3D

from datetime import datetime as dt
from datetime import timedelta


def scatter_animation_3d(locations, plane=None, step_size=86400.0, show=False):
    """
    Make a 3D scatter plot animation.

    :type locations: List
    :param locations: List of Location or Geographic objects
    :type plane: List
    :param plane: List of Location of Geographic defining a plane
    :type step_size: float
    :param step_size: Time in seconds for each frame.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(azim=-100)
    # Set up the writer
    Writer = matplotlib.animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    if plane:
        xs, ys, zs = ([0] * len(plane), [0] * len(plane), [0] * len(plane))
        for i, corner in enumerate(plane):
            if isinstance(corner, Geographic):
                xs[i] = corner.longitude
                ys[i] = corner.latitude
                zs[i] = corner.depth
            elif isinstance(corner, Location):
                xs[i] = corner.x
                ys[i] = corner.y
                zs[i] = corner.z
            else:
                raise TypeError("Plane is neither Location nor Geographic")
        ax.plot(xs, ys, zs, color='b')

    x, y, z, times, magnitudes = ([0] * len(locations), [0] * len(locations),
                                  [0] * len(locations), [0] * len(locations),
                                  [0] * len(locations))
    # Sort locations
    locations.sort(key=lambda l: l.time)
    for i, location in enumerate(locations):
        if isinstance(location, Geographic):
            x[i] = location.longitude
            y[i] = location.latitude
            z[i] = location.depth
            times[i] = location.time
            magnitudes[i] = location.magnitude
        elif isinstance(location, Location):
            x[i] = location.x
            y[i] = location.y
            z[i] = location.z
            times[i] = location.time
            magnitudes[i] = location.magnitude
        else:
            raise TypeError("Location is neither Geographic nor Location")
    time_label = ax.text(min(x), max(y), 0, "")

    # Make time-bins
    n_bins = (times[-1] - times[0]).total_seconds() / step_size
    n_bins = int(n_bins) + 1
    bins = []
    prev_stop = 0
    for bin_number in range(n_bins):
        _binstart = times[0] + timedelta(seconds=(step_size * bin_number))
        _bin = {"x": [], "y": [], "z": [], "times": [], "magnitudes": [],
                "bin_start": _binstart.strftime("%Y-%m-%dT%H:%M:%S.%f")}
        _binend = _binstart + timedelta(seconds=step_size)
        for i in range(prev_stop, len(locations)):
            if _binstart < times[i] < _binend:
                _bin["x"].append(x[i])
                _bin["y"].append(y[i])
                _bin["z"].append(z[i])
                _bin['times'].append(times[i])
                _bin['magnitudes'].append(magnitudes[i])
                prev_stop += 1
            elif times[i] > _binend:
                break
        bins.append(_bin)

    # Add background scatter
    # bg = ax.scatter(x, y, z, s=40, facecolors='white', edgecolors='black',
    #                 linewidth=0.5, alpha=0.5)
    # Set-up scatter objects
    sc = ax.scatter([],[],[], c='darkblue', s=40, alpha=1.0)
    sc1 = ax.scatter([], [], [], c='darkblue', s=35, alpha=0.9)
    sc2 = ax.scatter([], [], [], c='darkblue', s=30, alpha=0.8)
    sc3 = ax.scatter([], [], [], c='darkblue', s=25, alpha=0.7)
    sc4 = ax.scatter([], [], [], c='darkblue', s=20, alpha=0.6)
    sc5 = ax.scatter([], [], [], c='darkblue', s=10, alpha=0.5)
    sc6 = ax.scatter([], [], [], c='darkblue', s=15, alpha=0.4)
    sc7 = ax.scatter([], [], [], c='darkblue', s=5, alpha=0.3)
    sc8 = ax.scatter([], [], [], c='darkblue', s=2, alpha=0.2)
    scatters = [sc, sc1, sc2, sc3, sc4, sc5, sc6, sc7, sc8]

    def update(i):
        time_label.set_text(bins[i]["bin_start"])
        for j, scatter in enumerate(scatters):
            try:
                scatter._offsets3d = (
                    bins[i - j]["x"], bins[i - j]["y"], bins[i - j]["z"])
            except IndexError():
                pass

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ani = matplotlib.animation.FuncAnimation(
        fig, update, frames=len(bins), interval=10)

    plt.tight_layout()
    if show:
        plt.show()
    else:
        ani.save("scatter_animation_3d.mp4", writer=writer)

if __name__ == "__main__":
    locations = []
    with open("locations.csv", "r") as f:
        header = f.readline()
        for line in f:
            line = line.rstrip().split(', ')
            locations.append(Geographic(
                latitude=float(line[0]), longitude=float(line[1]),
                depth=-1 * float(line[2]), magnitude=float(line[3]),
                time=dt.strptime(line[4], "%Y-%m-%dT%H:%M:%S.%fZ")))
    origin = Geographic(latitude=-44.056691, longitude=168.723146,
                        depth=-0.015) # Location of Whataroa valley
    strike = 57.26  # Rotation of Alpine Fault clockwise from North
    dip = 50.0  # Dip of Alpine Fault from horizontal.
    # Make a plane to plot on the figures
    alpine_fault = get_plane(origin, strike, dip, length=300)

    scatter_animation_3d(locations, plane=alpine_fault)