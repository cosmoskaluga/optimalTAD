import numpy as np

from optimalTAD.optimization.staircaller import deltaH_statistic, get_bootstraps


def test_bootstrap_is_reproducible_with_same_seed():
    samples = (
        np.array([0.1, 0.4, 0.9, 1.2, 1.8, 2.3]),
        np.array([0.2, 0.3, 0.7, 1.0, 1.1, 1.5]),
    )

    first = get_bootstraps(samples, deltaH_statistic, rng=np.random.default_rng(1234))
    second = get_bootstraps(samples, deltaH_statistic, rng=np.random.default_rng(1234))

    assert first.confidence_interval.low == second.confidence_interval.low
    assert first.confidence_interval.high == second.confidence_interval.high
