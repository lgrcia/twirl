from typing import Optional

import numpy as np
from scipy.spatial import cKDTree

from twirl.geometry import get_transform_matrix, pad
from twirl.quads import hashes as hash4
from twirl.triangles import hashes as hash3


def count_cross_match(coords1, coords2, tol=1e-3):
    """
    Counts the number of cross-matches between two sets of 2D points.

    Parameters
    ----------
    coords1 : np.ndarray
        First set of points coordinates, shape (n, 2).
    coords2 : np.ndarray
        Second set of points coordinates, shape (m, 2).
    tol : float, optional
        Tolerance of the match, by default 1e-3.

    Returns
    -------
    int
        Number of cross-matches between coords1 and coords2.
    """
    return np.count_nonzero(
        np.min(np.linalg.norm(coords1[:, None, :] - coords2[None, :, :], axis=-1), 0)
        < tol
    )


def cross_match(coords1, coords2, tolerance=10):
    """
    Finds the closest matches between two sets of 2D points.

    Parameters
    ----------
    coords1 : np.ndarray
        First set of points coordinates, shape (n, 2).
    coords2 : np.ndarray
        Second set of points coordinates, shape (m, 2).
    tolerance : float, optional
        Tolerance of the match, given in coords1 points units, by default 10.

    Returns
    -------
    np.ndarray
        Array of matched indices, where each row contains the indices of the matched points in coords1 and coords2.
    """

    matches = []

    for i, s in enumerate(coords1):
        distances = np.linalg.norm(s - coords2, axis=1)
        closest = np.argmin(distances)
        if distances[closest] < tolerance:
            matches.append([i, closest])

    return np.array(matches)


def find_transform(
    radecs: np.ndarray,
    pixels: np.ndarray,
    min_match: float = 0.7,
    asterism: int = 4,
    quads_tolerance: float = 0.02,
    tolerance: float = 12,
) -> np.ndarray:
    """
    Finds the transformation matrix from `radecs` to `pixels`.

    Parameters
    ----------
    radecs : np.ndarray
        The coordinates to be transformed, shape (n, 2).
    pixels : np.ndarray
        The target coordinates, shape (m, 2).
    min_match : float, optional
        The fraction of `pixels` coordinates that must be matched to stop the search.
        I.e., if the number of matched points is `>= min_match * len(pixels)`, the
        search stops and return the found transform. By default 0.7.
    asterism : int, optional
        The asterism to use for hashing, either 3 or 4, by default 4.
    quads_tolerance : float, optional
        The minimum euclidean distance between two quads to be matched and tested.
        By default 0.02.
    tolerance : float, optional
        The minimum distance between two coordinates to be considered cross-matched
        (in `pixels` units). This serves to compute the number of coordinates being
        matched between `radecs` and `pixels` for a given transform.
        By default 12.

    Returns
    -------
    np.ndarray
        The transformation matrix from `radecs` to `pixels`.
    """

    if asterism == 3:
        asterism_function = hash3
    elif asterism == 4:
        asterism_function = hash4
    else:
        raise ValueError("available asterisms are 3 and 4")

    hashes_pixels, asterism_pixels = asterism_function(pixels)
    hashes_radecs, asterism_radecs = asterism_function(radecs)

    tree_pixels = cKDTree(hashes_pixels)
    tree_radecs = cKDTree(hashes_radecs)
    pairs = []

    ball_query = tree_pixels.query_ball_tree(tree_radecs, r=quads_tolerance)
    matches = []

    for i, j in enumerate(ball_query):
        if len(j) > 0:
            pairs += [[i, k] for k in j]

    for i, j in pairs:
        M = get_transform_matrix(asterism_radecs[j], asterism_pixels[i])
        test = (M @ pad(radecs).T)[0:2].T
        match = count_cross_match(pixels, test, tolerance)
        matches.append(match)

        if min_match is not None:
            if match >= min_match * len(pixels):
                break

    if len(matches) == 0:
        return None
    else:
        i, j = pairs[np.argmax(matches)]
        M = get_transform_matrix(asterism_radecs[j], asterism_pixels[i])
        return M
