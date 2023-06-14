import numpy as np
import pytest

from twirl.geometry import pad, transform_matrix
from twirl.match import count_cross_match, find_transform


@pytest.mark.parametrize("asterism", [3, 4])
@pytest.mark.parametrize("n", [12, 7, 6])
@pytest.mark.parametrize("seed", [5, 3])
def test_true_geometric_match(asterism, n, seed, extra=5):
    """Test perfect match between transformation matrix
    from two sets of simulated points

    Testing :code:`M = find_transform(xy1, xy2)`

    Parameters
    ----------
    n : int, optional
        number of simulated points xy1, by default 15
    extra : int, optional
        number of extra points in xy2 not in xy1, by default 5
    """
    np.random.seed(seed)
    xy1 = np.random.rand(n, 2)
    true_M = transform_matrix(scale=8.0, rotation=np.pi, translation=(0.3, 0.1))
    xy2 = (true_M @ pad(np.array([*xy1, *np.random.rand(extra, 2)])).T)[0:2].T
    np.random.shuffle(xy2)

    M = find_transform(xy1, xy2, tolerance=0.01, asterism=asterism)
    np.testing.assert_allclose(M, true_M, atol=1e-6)


@pytest.mark.parametrize("asterism", [3, 4])
@pytest.mark.parametrize("n", [12, 7, 6])
@pytest.mark.parametrize("seed", [5, 3])
def test_fuzzy_match(asterism, n, seed):
    """Test perfect match between transformation matrix from
    two sets of simulated points, with positions errors

    Testing :code:`M = find_transform(xy1, xy2)

    Parameters
    ----------@pytest.mark.parametrize("asterism", [3, 4])

    n : int, optional
        number of simulated points xy1, by default 15
    extra : int, optional
        number of extra points in xy2 not in xy1, by default 5
    """
    np.random.seed(seed * 2)
    extra = 5
    xy1 = np.random.rand(n, 2)
    true_M = transform_matrix(scale=8.0, rotation=np.pi, translation=(0.3, 0.1))
    xy2 = (true_M @ pad(np.array([*xy1, *np.random.rand(extra, 2)])).T)[0:2].T
    np.random.shuffle(xy2)
    xy2 += 0.01 * np.random.rand(len(xy2), 2)

    M = find_transform(xy1, xy2, tolerance=0.02, asterism=asterism)
    cn = count_cross_match((M @ pad(xy1).T)[0:2].T, xy2, tol=0.02)
    assert cn > 0.8 * n


@pytest.mark.parametrize("n", [25, 30])
@pytest.mark.parametrize("seed", [5, 6])
def test_large_fuzzy_match(n, seed, asterism=3):
    """Test perfect match between transformation matrix from
    two sets of simulated points, with positions errors

    Testing :code:`M = find_transform(xy1, xy2)

    Parameters
    ----------@pytest.mark.parametrize("asterism", [3, 4])

    n : int, optional
        number of simulated points xy1, by default 15
    extra : int, optional
        number of extra points in xy2 not in xy1, by default 5
    """
    np.random.seed(seed)
    extra = 10
    xy1 = np.random.rand(n, 2)
    true_M = transform_matrix(scale=8.0, rotation=np.pi, translation=(0.3, 0.1))
    xy2 = (true_M @ pad(np.array([*xy1, *np.random.rand(extra, 2)])).T)[0:2].T
    np.random.shuffle(xy2)
    xy2 += 0.01 * np.random.rand(len(xy2), 2)

    M = find_transform(xy1, xy2, tolerance=0.02, asterism=asterism)
    cn = count_cross_match((M @ pad(xy1).T)[0:2].T, xy2, tol=0.02)
    assert cn > 0.8 * n


def test_realistic_match():
    pixels = np.array(
        [
            [
                1874.23682135,
                906.47045569,
                744.72233708,
                862.92482627,
                1554.26186158,
                377.22820759,
                31.92397294,
                66.8977319,
                1705.18939397,
                386.0952344,
                1435.64933862,
                1353.8251716,
                1534.62432011,
                425.8644209,
                386.17091984,
            ],
            [
                946.4496088,
                460.92037899,
                647.68930592,
                867.8975847,
                22.10989145,
                699.83604197,
                1521.24951425,
                1924.84590848,
                764.02973769,
                396.82036351,
                525.75953215,
                1201.6207975,
                493.05218508,
                1869.49881979,
                114.5057253,
            ],
        ]
    ).T

    radecs = np.array(
        [
            [
                274.73148949,
                274.86490329,
                274.84561483,
                274.85110499,
                274.90849086,
                274.90721201,
                274.75127985,
                274.78820829,
                274.77117239,
                274.79333541,
                274.75682744,
                274.87694318,
                274.76899521,
                274.80959496,
                274.92856186,
            ],
            [
                -68.14639379,
                -68.15990274,
                -68.16806503,
                -68.15018564,
                -68.15770553,
                -68.17102645,
                -68.15448174,
                -68.12906399,
                -68.16645804,
                -68.13535781,
                -68.13438334,
                -68.14672887,
                -68.16049091,
                -68.17644741,
                -68.12362587,
            ],
        ]
    ).T

    M = find_transform(radecs, pixels, tolerance=5, asterism=4)
    assert count_cross_match(pixels, (M @ pad(radecs).T)[0:2].T, tol=5) == 9
