import numpy as np


def pad(x):
    return np.hstack([x, np.ones((x.shape[0], 1))])


def transform_matrix(scale=1.0, rotation=0.0, translation=(0.0, 0.0)):
    """Affine transformation matix

    Parameters
    ----------
    scale : float, optional
        scaling factor, by default 1.
    rotation : float, optional
        rotation angle in radians, by default 0.
    translation : tuple, optional
        translation vector, by default (0., 0.)

    Returns
    -------
    np.ndarray
        transformation matrix
    """

    R = np.array(
        [
            [np.cos(rotation), np.sin(-rotation), 0],
            [np.sin(rotation), np.cos(rotation), 0],
            [0, 0, 1],
        ]
    )
    T = np.array([[1, 0, translation[0]], [0, 1, translation[1]], [0, 0, 1]])
    S = np.array([[scale, 0, 0], [0, scale, 0], [0, 0, 1]])
    return T @ S @ R


def proj(p, origin, axe, norm=False):
    """
    projection of a point p on a segment from origin to axe
    """
    n = axe - origin
    n /= np.linalg.norm(n, 2)
    return np.where(norm, np.dot(p - origin, n), origin + n * np.dot(p - origin, n))[0]


def rotation(point, pivot, angle, norm=False):
    """
    rotate point around pivot of certain angle
    """
    co = np.cos(angle)
    si = np.sin(angle)
    r = np.array([[co, -si], [si, co]])
    x = r @ (point - pivot)
    x = np.where(norm, (x / np.linalg.norm(x)) * np.linalg.norm(point - pivot) * co, x)
    return x + pivot


def u1u2(a, b, norm=False):
    """
    x, y basis as defined in Lang2009
    """
    norm = np.where(norm, np.linalg.norm(b - a), False)
    x = rotation(b, a, -np.pi / 4, norm=norm)
    y = rotation(b, a, np.pi / 4, norm=norm)
    return x, y


def _get_transform_matrix(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Compute the transformation matrix between two vectors defined by pairs of points coordinates

    Parameters
    ----------
    v1 : np.array
        pair of points coordinates, shape (2, 2)
    v2 :  np.array
        pair of points coordinates, shape (2, 2)

    Returns
    -------
    np.ndarray
        transformation matrix of shape (3, 3)
    """
    u1 = v1[1] - v1[0]
    u2 = v2[1] - v2[0]
    n1 = np.linalg.norm(u1)
    n2 = np.linalg.norm(u2)
    theta = np.arccos(np.clip(np.dot(u1 / n1, u2 / n2), -1.0, 1.0))
    scale = n2 / n1
    A = transform_matrix(scale=scale, rotation=theta)
    a = np.pad(v1[0], (0, 1), constant_values=1)
    b = np.pad(v2[0], (0, 1), constant_values=1)
    t = (b - A @ a)[0:2]
    return transform_matrix(scale=scale, rotation=theta, translation=t)


def get_transform_matrix(xy1, xy2):
    XY1 = pad(xy1)
    XY2 = pad(xy2)
    M, _, _, _ = np.linalg.lstsq(XY1, XY2, rcond=None)
    return M.T


def triangle_angles(trios):
    if trios.shape[1:] != (3, 2):
        raise ValueError("The input array must have shape (n, 3, 2)")

    # Calculate the vectors between the points
    vec1 = trios[:, 1, :] - trios[:, 0, :]
    vec2 = trios[:, 2, :] - trios[:, 1, :]
    vec3 = trios[:, 0, :] - trios[:, 2, :]

    # Calculate the lengths of the vectors (distances between the points)
    a = np.linalg.norm(vec2, axis=1)
    b = np.linalg.norm(vec3, axis=1)
    c = np.linalg.norm(vec1, axis=1)

    # Use the law of cosines to find the angles
    angle_A = np.arccos((b**2 + c**2 - a**2) / (2 * b * c))
    angle_B = np.arccos((c**2 + a**2 - b**2) / (2 * c * a))
    angle_C = np.arccos((a**2 + b**2 - c**2) / (2 * a * b))

    # Convert to degrees and return
    return np.array([angle_A, angle_B, angle_C]).T


def sparsify(coords: np.ndarray, radius: float) -> np.ndarray:
    """
    Returns a sparsified version of the input coordinates array, where only the points that are at least `radius` distance
    apart from each other are kept.

    Parameters:
    -----------
    coords : np.ndarray
        The input array of shape (n, 2) containing the coordinates of the points.
    radius : float
        The minimum distance between two points to be kept in the output array.

    Returns:
    --------
    np.ndarray
        The sparsified array of shape (m, 2), where `m` is the number of points that are at least `radius` distance apart
        from each other.
    """
    _coords = coords.copy()
    deleted_coords = np.zeros([], dtype=int)
    sparse_coords = []

    for i, s in enumerate(_coords):
        if not i in deleted_coords:
            distances = np.linalg.norm(_coords - s, axis=1)
            idxs = np.flatnonzero(distances < radius)
            sparse_coords.append(s)
            deleted_coords = np.hstack([deleted_coords, idxs])

    return np.array(sparse_coords)
