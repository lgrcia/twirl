from itertools import combinations

import numpy as np

from twirl.geometry import triangle_angles


def order_points(triangles):
    """
    Orders the vertices of each triangle in a consistent manner.

    Vertices are ordered by their distance to the centroid of the triangle.

    Parameters
    ----------
    triangles : ndarray
        An array of shape (n_triangles, 3, 2) representing the vertices of each triangle.

    Returns
    -------
    ndarray
        An array of shape (n_triangles, 3, 2) with the vertices of each triangle ordered consistently.
    """

    # Compute the centroid of the triangles
    centroids = np.mean(triangles, axis=1, keepdims=True)

    # Compute the distances from the centroid to the vertices
    distances = np.linalg.norm(triangles - centroids, axis=-1)

    # Get the indices that would sort the distances
    sort_indices = np.argsort(distances, axis=1)

    # Use numpy's advanced indexing to sort the triangles
    ordered_triangles = np.take_along_axis(triangles, sort_indices[:, :, None], axis=1)

    return ordered_triangles


def hashes(xy, min_angle=np.deg2rad(30)):
    """
    Computes the hashes of the triangles formed by the points in xy.

    Parameters
    ----------
    xy : ndarray
        An array of shape (n_points, 2) representing the x and y coordinates of each point.
    min_angle : float, optional
        The minimum angle (in radians) that a triangle must have to be included in the hashes. Default is 30 degrees.

    Returns
    -------
    hashes : ndarray
        An array of shape (n_hashes, 2) representing the hashes of the triangles.
    triangles : ndarray
        An array of shape (n_triangles, 3, 2) representing the vertices of each triangle.
    """

    triangles_idxs = np.array(list(combinations(np.arange(xy.shape[0]), 3)))
    triangles = xy[triangles_idxs]
    triangles = order_points(triangles)
    angles = triangle_angles(triangles)
    # keep only triangles with any angle > min_angle
    mask = np.all(np.abs(angles) > min_angle, axis=1)
    triangles = triangles[mask]
    angles = angles[mask]
    hashes = np.sort(angles, axis=1)[:, 0:2]
    return hashes, triangles
