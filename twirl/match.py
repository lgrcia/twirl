import numpy as np
from .geometry import get_transform_matrix, pad
from .quads import hashes


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
    xy2: np.ndarray, xy1: np.ndarray, tolerance: int = 12, min_match: int = None
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

    hash1, quads1 = hashes(xy1)
    hash2, quads2 = hashes(xy2)
    distances = np.linalg.norm(hash1[:, None, :] - hash2[None, :, :], axis=2)
    shortest_hash = np.argmin(distances, 1)
    ns = []

    for i, j in enumerate(shortest_hash):
        M = get_transform_matrix(quads2[j], quads1[i])
        test = (M @ pad(xy2).T)[0:2].T
        n = count_cross_match(xy1, test, tolerance)
        if min_match is not None:
            if n >= min_match:
                break
        else:
            ns.append(n)

        i = np.argmax(ns)
        M = get_transform_matrix(quads2[np.argmin(distances, 1)[i]], quads1[i])

    return M
