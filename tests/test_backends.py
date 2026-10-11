import numpy as np
import pytest
from astropy.wcs import WCS

from twirl.geometry import pad, transform_matrix
from twirl.match import find_transform
from twirl.quads import reorder
from twirl.utils import compute_wcs

pytest.importorskip("numba")


def reorder_loop(quads):
    """Reference implementation of `twirl.quads.reorder` (before vectorization)."""
    distances = np.linalg.norm(
        quads[:, :, :, None] - np.rollaxis(quads.T, 2)[:, None, :, :], axis=2
    )
    i = np.argmax(np.max(distances, 1), 1)
    idxs = np.roll(
        np.argsort(distances[range(len(quads)), i], axis=1)[:, ::-1], 1, axis=1
    )
    return np.array([q[i] for q, i in zip(quads, idxs)])


def simulated_points(n, seed, extra=5, noise=0.0):
    rng = np.random.default_rng(seed)
    xy1 = rng.random((n, 2))
    true_M = transform_matrix(scale=8.0, rotation=np.pi, translation=(0.3, 0.1))
    xy2 = (true_M @ pad(np.vstack([xy1, rng.random((extra, 2))])).T)[0:2].T
    rng.shuffle(xy2)
    xy2 += noise * rng.random(xy2.shape)
    return xy1, xy2


def test_reorder_matches_loop():
    rng = np.random.default_rng(0)
    quads = rng.random((500, 4, 2))
    np.testing.assert_array_equal(reorder(quads), reorder_loop(quads))


@pytest.mark.parametrize("asterism", [3, 4])
@pytest.mark.parametrize("n", [6, 12, 25])
@pytest.mark.parametrize("seed", [3, 5])
@pytest.mark.parametrize("min_match", [0.7, None])
@pytest.mark.parametrize("noise", [0.0, 0.01])
def test_backends_same_transform(asterism, n, seed, min_match, noise):
    xy1, xy2 = simulated_points(n, seed, noise=noise)
    kwargs = dict(tolerance=0.02, asterism=asterism, min_match=min_match)
    M_numpy = find_transform(xy1, xy2, backend="numpy", **kwargs)
    M_numba = find_transform(xy1, xy2, backend="numba", **kwargs)
    np.testing.assert_allclose(M_numba, M_numpy, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("asterism", [3, 4])
def test_backends_same_wcs(asterism):
    rng = np.random.default_rng(1)
    true_wcs = WCS(naxis=2)
    true_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    true_wcs.wcs.crval = [150.0, -30.0]
    true_wcs.wcs.crpix = [1024, 1024]
    angle = 0.7
    true_wcs.wcs.cd = (0.35 / 3600) * np.array(
        [[-np.cos(angle), np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    xy = rng.uniform(0, 2048, (30, 2))
    radecs = np.array(true_wcs.pixel_to_world_values(*xy.T)).T
    pixels = xy[:20] + rng.normal(0, 0.5, (20, 2))

    wcs_numpy = compute_wcs(pixels, radecs, asterism=asterism, backend="numpy")
    wcs_numba = compute_wcs(pixels, radecs, asterism=asterism, backend="numba")
    np.testing.assert_allclose(
        wcs_numba.pixel_to_world_values(*xy.T),
        wcs_numpy.pixel_to_world_values(*xy.T),
        rtol=0,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        wcs_numba.world_to_pixel_values(*radecs.T), xy.T, atol=1.0
    )


def test_unknown_backend():
    xy1, xy2 = simulated_points(8, 0)
    with pytest.raises(ValueError):
        find_transform(xy1, xy2, backend="cython")
