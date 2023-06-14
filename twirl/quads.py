from itertools import combinations

import numpy as np

from twirl.geometry import proj, u1u2


def reorder(quads):
    distances = np.linalg.norm(
        quads[:, :, :, None] - np.rollaxis(quads.T, 2)[:, None, :, :], axis=2
    )
    i = np.argmax(np.max(distances, 1), 1)
    idxs = np.roll(
        np.argsort(distances[range(len(quads)), i], axis=1)[:, ::-1], 1, axis=1
    )
    return np.array([q[i] for q, i in zip(quads, idxs)])


def good_quads(quads, circletol=0.01):
    """
    whether all points are contained in a circle (see Lang2009)
    """
    q0, q1 = quads[:, 0], quads[:, 1]
    r = np.linalg.norm(q1 - q0, axis=1) / 2
    center = q0 + (q1 - q0) / 2
    in_circle = np.linalg.norm(quads - center[:, None, :], axis=2) <= r[:, None] * (
        1 + circletol
    )
    return np.all(in_circle, axis=1)


def quad_hash(quads):
    a, b, c, d = np.rollaxis(quads, 1)

    norm = np.linalg.norm(b - a, axis=1)

    def u1u2(a, b):
        """compute x, y basis as defined in Lang2009"""
        co_m = np.cos(-np.pi / 4)
        si_m = np.sin(-np.pi / 4)
        co_p = np.cos(np.pi / 4)
        si_p = np.sin(np.pi / 4)
        rm = np.array([[co_m, -si_m], [si_m, co_m]])
        rp = np.array([[co_p, -si_p], [si_p, co_p]])
        x = (rm @ (b - a).T).T + a
        y = (rp @ (b - a).T).T + a
        return x, y

    u1, u2 = u1u2(a, b)

    def proj(p, origin, axe):
        """projection of a point p on a segment from origin to axe"""
        n = axe - origin
        n /= np.linalg.norm(n, axis=1)[:, None]
        return np.sum((p - origin) * n, 1)

    return (
        np.array(
            [
                proj(c, a, u1) / norm,
                proj(d, a, u1) / norm,
                proj(c, a, u2) / norm,
                proj(d, a, u2) / norm,
            ]
        ).T,
        np.rollaxis(np.array([a, b]), 1),
    )


def clean_quads(xy):
    assert xy.shape[1] == 2
    quads_idxs = np.array(list(combinations(np.arange(xy.shape[0]), 4)))
    quads = xy[quads_idxs]
    ordered_quads = reorder(quads)
    return ordered_quads[good_quads(ordered_quads)]


def hashes(xy):
    assert xy.shape[1] == 2
    quads_idxs = np.array(list(combinations(np.arange(xy.shape[0]), 4)))
    quads = xy[quads_idxs]
    ordered_quads = reorder(quads)
    good_ordered_quads = ordered_quads[good_quads(ordered_quads)]
    h, u = quad_hash(good_ordered_quads)
    # we sort hashes from larger AB (see Lang 2008)
    idxs = np.argsort(np.linalg.norm(u[:, 1] - u[:, 0], axis=1))
    return h[idxs[::-1]], good_ordered_quads[idxs[::-1]]
