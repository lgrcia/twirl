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


def _match_counts_numpy(
    asterism_radecs, asterism_pixels, i, j, radecs, pixels, tolerance, stop
):
    """Pure NumPy version of `twirl._numba.match_counts`."""
    padded_radecs = pad(radecs)
    counts = []
    for p in range(len(i)):
        M = get_transform_matrix(asterism_radecs[j[p]], asterism_pixels[i[p]])
        test = (M @ padded_radecs.T)[0:2].T
        counts.append(count_cross_match(pixels, test, tolerance))
        if stop >= 0 and counts[-1] >= stop:
            break
    return np.array(counts, dtype=int)


def find_transform(
    radecs: np.ndarray,
    pixels: np.ndarray,
    min_match: float = 0.7,
    asterism: int = 4,
    quads_tolerance: float = 0.02,
    tolerance: float = 12,
    backend: str = "auto",
) -> Optional[np.ndarray]:
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
    backend : str, optional
        How candidate asterism pairs are tested: "numba" (compiled, much faster,
        requires numba), "numpy", or "auto" to use numba when it can be imported.
        All backends return the same transform. The first numba call in a process
        compiles the loop (about 1 s, then cached on disk). By default "auto".

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

    if backend == "auto":
        # numba can be installed but fail to import (e.g. unsupported NumPy)
        try:
            from twirl._numba import match_counts
        except ImportError:
            match_counts = _match_counts_numpy
    elif backend == "numba":
        from twirl._numba import match_counts
    elif backend == "numpy":
        match_counts = _match_counts_numpy
    else:
        raise ValueError("available backends are 'auto', 'numba' and 'numpy'")

    hashes_pixels, asterism_pixels = asterism_function(pixels)
    hashes_radecs, asterism_radecs = asterism_function(radecs)

    tree_pixels = cKDTree(hashes_pixels)
    tree_radecs = cKDTree(hashes_radecs)
    ball_query = tree_pixels.query_ball_tree(tree_radecs, r=quads_tolerance)

    # candidate pairs (pixels asterism i, radecs asterism j), as arrays
    lengths = np.fromiter(map(len, ball_query), dtype=np.int64, count=len(ball_query))
    i = np.repeat(np.arange(len(ball_query), dtype=np.int64), lengths)
    j = np.fromiter(
        (k for js in ball_query for k in js), dtype=np.int64, count=int(lengths.sum())
    )

    if len(i) == 0:
        return None

    stop = -1.0 if min_match is None else float(min_match * len(pixels))
    matches = match_counts(
        np.ascontiguousarray(asterism_radecs, dtype=float),
        np.ascontiguousarray(asterism_pixels, dtype=float),
        i,
        j,
        np.ascontiguousarray(radecs, dtype=float),
        np.ascontiguousarray(pixels, dtype=float),
        float(tolerance),
        stop,
    )

    best = np.argmax(matches)
    return get_transform_matrix(asterism_radecs[j[best]], asterism_pixels[i[best]])
