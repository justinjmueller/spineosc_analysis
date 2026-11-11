import numpy as np
import pytest

from fdgen import dials

S = np.array([-3., -2., -1., 0., 1., 2., 3.])


def test_spline_passes_through_every_knot():
    rng = np.random.default_rng(0)
    y = 1 + 0.3 * rng.standard_normal((5, 7))
    y[:, 3] = 1.0
    for k, x in enumerate(S):
        assert np.allclose(dials.cafana_spline(x, S, y), y[:, k])


def test_exact_for_a_quadratic_on_even_knots():
    y = (1 + 0.2 * S + 0.05 * S ** 2)[None, :]
    for x in (-2.5, -0.4, 0.3, 1.7, 2.9):
        assert np.allclose(dials.cafana_spline(x, S, y), 1 + 0.2 * x + 0.05 * x ** 2)


def test_two_knots_are_linear_and_out_of_range_is_refused():
    y = np.array([[1.0, 1.4]])
    assert np.allclose(dials.cafana_spline(0.5, np.array([0., 1.]), y), 1.2)
    with pytest.raises(ValueError):
        dials.cafana_spline(1.5, np.array([0., 1.]), y)
