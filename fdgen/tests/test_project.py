import numpy as np

from fdgen import osc
from fdgen.generate import project_true_e

TE = osc.le_bin_edges(0.3, 3.0, 20)


def test_projection_is_a_weighted_mean_per_branch_and_true_e_bin():
    ev = dict(aw=np.array([1., 1, 2, 2, 0, 1]),
              true_e=np.array([0.5, 0.5, 0.5, 5.0, 0.5, 0.31], dtype=np.float32),
              branch=np.array([0, 0, 1, 0, 0, 0]))
    w = np.array([2., 4, 3, 9, 7, 1.])
    # branch 0, bin of 0.5: (1*2 + 1*4) / 2 = 3; branch 1 is its own bin; E = 5 is outside the bins; aw = 0 stays 1;
    # E = 0.31 sits in another bin alone with ratio w = 1
    assert project_true_e(w, ev, TE).tolist() == [3., 3., 3., 1., 1., 1.]


def test_projection_preserves_the_weighted_total_inside_the_bins():
    rng = np.random.default_rng(1)
    n = 5000
    ev = dict(aw=rng.uniform(0.5, 1.5, n), true_e=rng.uniform(0.3, 3.0, n).astype(np.float32), branch=rng.integers(0, 3, n))
    w = rng.uniform(0.5, 2.0, n)
    p = project_true_e(w, ev, TE)
    assert np.isclose((ev["aw"] * p).sum(), (ev["aw"] * w).sum())
