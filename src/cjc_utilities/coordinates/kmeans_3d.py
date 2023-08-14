"""
Example clustering in 3D using scikit-learn kmeans clustering
"""
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D

from coordinates import Location, Geographic

from sklearn.cluster import KMeans

def kmeans_3d(locations, n_clusters):
    """
    Cluster a list of locations using k_means.

    :type locations: list
    :param locations: List of Location
    :type n_clusters: int
    :param n_clusters: Number of cluster centroids. 
    
    :returns: list of lists of clusters
    """
    estimator = KMeans(n_clusters=n_clusters)
    _locations = np.array([np.array([l.x, l.y, l.z]) for l in locations])
    estimator.fit(_locations)
    clusters = []
    for i in range(n_clusters):
        cluster = [l for j, l in enumerate(locations) 
                   if estimator.labels_[j] == i]
        clusters.append(cluster)
    return clusters
    
def plot_clusters(clusters):
    """
    Plot clusters in different colours.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colours = iter(cm.rainbow(np.linspace(0, 1, len(clusters))))
    for cluster in clusters:
        x, y, z = ([0] * len(cluster), [0] * len(cluster),
                   [0] * len(cluster))
        for i, location in enumerate(cluster):
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
        ax.scatter(x, y, z, c=next(colours), marker='o')
    if isinstance(clusters[0][0], Geographic):
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
        ax.set_zlabel("Depth (km) -ve down")
    else:
        ax.set_xlabel('Distance normal to plane (km)')
        ax.set_ylabel('Distance along strike (km)')
        ax.set_zlabel('Distance down-dip (km) -ve down')
    plt.show()


if __name__ == '__main__':
    # Read in the konstantinos dataset for testing.
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
    projected = [location.to_xyz(origin, strike=strike, dip=dip)
                 for location in locations]
    clusters = kmeans_3d(projected, 8)
    _clusters = []
    for cluster in clusters:
        cluster = [l.to_geographic() for l in cluster]
        _clusters.append(cluster)
    clusters = _clusters
    plot_clusters(clusters)