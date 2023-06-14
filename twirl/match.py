from typing import Optional

import numpy as np

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
    tolerance: int = 12,
    min_match: Optional[int] = None,
    asterism: int = 4,
    max_search: int = 10000,
) -> np.ndarray:
    """
    Finds the transformation matrix that maps the coordinates in `coords2` to the coordinates in `coords1`.

    Parameters
    ----------
    radecs : np.ndarray
        The coordinates to be transformed, shape (n, 2).
    pixels : np.ndarray
        The target coordinates, shape (m, 2).
    tolerance : int, optional
        The tolerance of the match, given in `coords1` points units, by default 12.
    min_match : Optional[int], optional
        The minimum number of matches required to stop the search, by default None.
    asterism : int, optional
        The asterism to use for hashing, either 3 or 4, by default 4.

    Returns
    -------
    np.ndarray
        The transformation matrix that maps `radecs` to `pixels`.
    """

    if asterism == 3:
        asterism_function = hash3
    elif asterism == 4:
        asterism_function = hash4
    else:
        raise ValueError("available asterisms are 3 and 4")

    hashes_pixels, asterism_pixels = asterism_function(pixels)
    hashes_radecs, asterism_radecs = asterism_function(radecs)
    distances = np.linalg.norm(
        hashes_pixels[:, None, :] - hashes_radecs[None, :, :], axis=2
    )
    best_couples = np.array(
        np.unravel_index(np.argsort(distances.flatten())[0:max_search], distances.shape)
    ).T
    ns = []

    for i, j in best_couples:
        M = get_transform_matrix(asterism_radecs[j], asterism_pixels[i])
        test = (M @ pad(radecs).T)[0:2].T
        n = count_cross_match(pixels, test, tolerance)
        ns.append(n)

        if min_match is not None:
            if n >= min_match:
                break

    i, j = best_couples[np.argmax(ns)]

    M = get_transform_matrix(asterism_radecs[j], asterism_pixels[i])

    return M
