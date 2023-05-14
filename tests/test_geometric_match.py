import numpy as np
import pytest

from twirl.geometry import pad, transform_matrix
from twirl.match import count_cross_match, find_transform


@pytest.mark.parametrize("asterism", [3, 4])
@pytest.mark.parametrize("n", [12, 10, 7, 6])
@pytest.mark.parametrize("seed", [10, 7, 5, 4, 3, 1, 2])
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
@pytest.mark.parametrize("n", [12, 10, 7, 6])
@pytest.mark.parametrize("seed", [10, 7, 5, 4, 3, 1, 2])
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


@pytest.mark.parametrize("n", [25, 30, 20, 22])
@pytest.mark.parametrize("seed", [10, 7, 5, 4, 3, 1, 2])
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
