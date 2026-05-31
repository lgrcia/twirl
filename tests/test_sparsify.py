import numpy as np

from twirl import sparsify


def test_sparsify_does_not_drop_first():
    # Regression test for issue lgarci/twirl#41

    # Note that these are 0.1 degree from each other
    input_coords = np.array([[1, 0], [1.1, 0], [1.2, 0]])

    # All coordinate separations are larger than 0.01 deg
    # so none should get dropped.
    output_coords = sparsify(input_coords, 0.01)

    # make sure none got dropped based on length
    assert output_coords.shape == input_coords.shape

    # Make sure the values were not changed
    assert np.all(output_coords == input_coords)
