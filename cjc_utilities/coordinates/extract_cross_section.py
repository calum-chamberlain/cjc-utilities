"""
Code to extract a subset of a catalog adjacent to an arbitrary cross-section.

:author: Calum Chamberlain
:date: 23/03/2018
"""

from cjc_utilities.coordinates.coordinates import Location, Geographic

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_plane(origin, strike, dip, length=100, height=-50):
    """
    Get four corners of a plane around an origin.

    :type origin: Geographic
    :type strike: float
    :param strike: Strike of plane in degrees
    :type dip: float
    :param dip: Dip of plane in degrees
    :type length: float
    :param length:
        Total length along strike, will add on the North side of the origin
    type height: float
    :param height: Total length down-dip, will be one-sided from origin.

    :return: List of Geographic
    """
    points = [Location(0, 0, 0, origin, strike, dip),
              Location(0, length, 0, origin, strike, dip),
              Location(0, length, height, origin, strike, dip),
              Location(0, 0, height, origin, strike, dip),
              Location(0, 0, 0, origin, strike, dip)]
    point_geographics = [point.to_geographic() for point in points]
    return point_geographics

def plot_xyz(locations, plane=None):
    """
    Make a 3D scatter plot of locations, either list of Location, or list of
    tuples of (x, y, z).
    """
    x, y, z = ([0] * len(locations), [0] * len(locations),
               [0] * len(locations))
    for i, location in enumerate(locations):
        if isinstance(location, Geographic):
            x[i] = location.longitude
            y[i] = location.latitude
            z[i] = location.depth
        elif isinstance(location, Location):
            x[i] = location.x
            y[i] = location.y
            z[i] = location.z
        else:
            raise TypeError("Location is neither Geographic nor Location")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='k', marker='o')
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
    if isinstance(locations[0], Geographic):
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
        ax.set_zlabel("Depth (km) -ve down")
    else:
        ax.set_xlabel('Distance normal to plane (km)')
        ax.set_ylabel('Distance along strike (km)')
        ax.set_zlabel('Distance down-dip (km) -ve down')
    return fig


def plot_xy(locations):
    """
    Plot 2D scatter plot coloured by depth.
    """
    x, y, z, s = ([0] * len(locations), [0] * len(locations),
                  [0] * len(locations), [1] * len(locations))
    for i, location in enumerate(locations):
        if isinstance(location, Geographic):
            x[i] = location.longitude
            y[i] = location.latitude
            z[i] = location.depth
            if location.magnitude is not None:
                s[i] = location.magnitude
        elif isinstance(location, Location):
            x[i] = location.x
            y[i] = location.y
            z[i] = location.z
            if location.magnitude is not None:
                s[i] = location.magnitude
        else:
            raise TypeError("Location is neither Geographic nor Location")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y, s=s, c=z)
    if isinstance(locations[0], Geographic):
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
    else:
        ax.set_xlabel('Distance normal to plane (km)')
        ax.set_ylabel('Distance along strike (km)')

    return fig


if __name__ == '__main__':
    locations = []
    with open("data/locations.csv", "r") as f:
        header = f.readline()
        for line in f:
            line = line.split(',')
            locations.append(Geographic(
                latitude=float(line[0]), longitude=float(line[1]),
                depth=-1 * float(line[2])))
    origin = Geographic(latitude=-44.056691, longitude=168.723146,
                        depth=-0.015) # Location of Whataroa valley
    strike = 57.26  # Rotation of Alpine Fault clockwise from North
    dip = 50.0  # Dip of Alpine Fault from horizontal.
    # Make a plane to plot on the figures
    alpine_fault = get_plane(origin, strike, dip, length=300)
    projected = [location.to_xyz(origin, strike=strike, dip=dip)
                 for location in locations]
    projected_alpine_fault = [corner.to_xyz(origin, strike=strike, dip=dip)
                              for corner in alpine_fault]
    fig = plot_xyz(locations, alpine_fault)
    plt.title("Input locations")
    fig.show()
    fig = plot_xyz(projected, projected_alpine_fault)
    plt.title("Projected locations")
    fig.show()
    # Plot just the +/- 2km swath that V wants
    swath = [location for location in projected if abs(location.x) < 2.0]
    fig = plot_xyz(swath)
    plt.title("Swath +/-km around the fault.")
    fig.show()
    # Rotate swath back to lat-lon-depth and plot on a map
    swath_latlondepth = [location.to_geographic() for location in swath]
    fig = plot_xyz(swath_latlondepth)
    plt.title("Swath in lat-lon")
    plt.show(block=True)

