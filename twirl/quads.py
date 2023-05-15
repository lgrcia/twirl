from itertools import combinations

import numpy as np

from twirl.geometry import proj, u1u2


def _reorder(q):
    """
    order coordinates from closest to first one
    """
    distances = np.linalg.norm(q[:, :, None] - q.T[None, :, :], axis=1)
    q0i = np.min(np.unravel_index(np.argmax(distances), distances.shape))
    q0 = q[q0i]
    q123 = np.delete(q.copy(), q0i, axis=0)
    distance_from_q0 = np.linalg.norm(q123 - q0, axis=1)
    idxs = np.argsort(1 / distance_from_q0)
    return np.array([q0, *q123[idxs]])


reorder = np.vectorize(_reorder, signature="(n, m)->(n, m)")


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


def _quad_hash(quad):
    """
    from 4 coordinates produce the quad hash code
    """
    a, b, c, d = quad
    u1, u2 = u1u2(a, b)
    h = np.linalg.norm(b - a)
    return np.array(
        [
            proj(c, a, u1, norm=True) / h,
            proj(d, a, u1, norm=True) / h,
            proj(c, a, u2, norm=True) / h,
            proj(d, a, u2, norm=True) / h,
        ]
    ), np.array([a, b])


quad_hash = np.vectorize(_quad_hash, signature="(n,m)->(k), (l, 2)")


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
