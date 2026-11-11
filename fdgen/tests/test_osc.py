import numpy as np

from fdgen import osc

EDGES = osc.le_bin_edges(0.0, 3.0, 300)


def test_bin_edges_follow_profit_convention():
    # lower edge belongs to its bin, the last edge and anything above or below is dropped, NaN is dropped
    idx = osc.le_bin_index(np.array([0.0, 0.0099, 0.01, 2.9999, 3.0, 3.5, -0.1, np.nan], dtype=np.float32), EDGES)
    assert idx.tolist() == [0, 0, 1, 299, -1, -1, -1, -1]


def test_probability_is_evaluated_at_the_bin_midpoint():
    le = np.array([0.1000, 0.1049, 0.1099], dtype=np.float32)       # three events in the same bin [0.10, 0.11)
    p = osc.pmumu_profit(le, EDGES, 3.0, 0.4)
    mid = np.float32(0.105)
    want = 1 - 0.4 * np.sin(1.266932679 * 3.0 * float(mid)) ** 2
    assert np.allclose(p, want, rtol=1e-6) and np.ptp(p) == 0


def test_no_oscillation_and_clamping():
    le = np.array([0.5, 1.5], dtype=np.float32)
    assert np.all(osc.pmumu_profit(le, EDGES, 3.0, 0.0) == 1.0)
    assert np.allclose(osc.pmumu_profit(le, EDGES, 3.0, 2.0), osc.pmumu_profit(le, EDGES, 3.0, 1.0))   # s clamped to 1
    assert np.allclose(osc.pmumu_profit(le, EDGES, 3.0, -1.0), 1.0)                                   # and to 0


def test_outside_grid_is_nan():
    p = osc.pmumu_profit(np.array([3.2, np.nan], dtype=np.float32), EDGES, 3.0, 0.4)
    assert np.isnan(p).all()
