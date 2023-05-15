import numpy as np

from twirl.geometry import get_transform_matrix, pad
from twirl.quads import hashes as hash4
from twirl.triangles import hashes as hash3


def count_cross_match(xy1, xy2, tol=1e-3):
    return np.count_nonzero(
        np.min(np.linalg.norm(xy1[:, None, :] - xy2[None, :, :], axis=-1), 0) < tol
    )


def cross_match(xy1, xy2, tolerance=10):
    matches = []

    for i, s in enumerate(xy1):
        distances = np.linalg.norm(s - xy2, axis=1)
        closest = np.argmin(distances)
        if distances[closest] < tolerance:
            matches.append([i, closest])

    return np.array(matches)


def find_transform(
    xy2: np.ndarray,
    xy1: np.ndarray,
    tolerance: int = 12,
    min_match: int = None,
    asterism: int = 4,
) -> np.ndarray:
    """_summary_

    Parameters
    ----------
    xy2 : np.ndarray
        first set of points coordinates, shape (n, 2)
    xy1 : np.ndarray
        second set of points coordinates, shape (m, 2)
    tolerance : int, optional
        tolerance of the match, given in xy1 points units, by default 12
    min_match : _type_, optional
        if None, stops when the number of crossed matched points
        is min_match, by default None

    Returns
    -------
    _type_
        _description_
    """

    if asterism == 3:
        asterism = hash3
    elif asterism == 4:
        asterism = hash4
    else:
        raise ValueError("available asterisms are 3 and 4")

    hash1, points1 = asterism(xy1)
    hash2, points2 = asterism(xy2)
    distances = np.linalg.norm(hash1[:, None, :] - hash2[None, :, :], axis=2)
    shortest_hash = np.argmin(distances, 1)
    ns = []

    for i, j in enumerate(shortest_hash):
        M = get_transform_matrix(points2[j], points1[i])
        test = (M @ pad(xy2).T)[0:2].T
        n = count_cross_match(xy1, test, tolerance)
        ns.append(n)
        if min_match is not None:
            if n >= min_match:
                break

        i = np.argmax(ns)
        M = get_transform_matrix(points2[np.argmin(distances, 1)[i]], points1[i])

    return M
