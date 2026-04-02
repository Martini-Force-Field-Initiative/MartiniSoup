import pytest
import numpy as np
from martinisoup.utils import binned_histogram, bootstrap_ci


class TestBinnedHistogram:

    def test_linear_bins_shape(self):
        data = np.linspace(1, 100, 50)
        result = binned_histogram(data, nbins=10, log_bins=False)
        assert len(result['counts']) == 10
        assert len(result['bin_centers']) == 10
        assert len(result['edges']) == 11

    def test_log_bins_shape(self):
        data = np.logspace(0, 2, 50)
        result = binned_histogram(data, nbins=10, log_bins=True)
        assert len(result['counts']) == 10
        assert len(result['bin_centers']) == 10
        assert len(result['edges']) == 11

    def test_linear_pdf_integrates_to_one(self):
        np.random.seed(0)
        data = np.random.exponential(scale=10, size=500)
        result = binned_histogram(data, nbins=20, log_bins=False)
        area = np.sum(result['pdf'] * np.diff(result['edges']))
        assert abs(area - 1.0) < 0.01

    def test_log_pdf_integrates_to_one(self):
        np.random.seed(0)
        data = np.random.exponential(scale=10, size=500)
        result = binned_histogram(data, nbins=20, log_bins=True)
        area = np.sum(result['pdf'] * np.diff(result['edges']))
        assert abs(area - 1.0) < 0.01

    def test_all_zero_data_log_bins_returns_empty(self):
        result = binned_histogram([0, 0, 0], nbins=10, log_bins=True)
        assert result['counts'] == []

    def test_linear_bin_centers_are_midpoints(self):
        data = np.linspace(0, 10, 100)
        result = binned_histogram(data, nbins=5, log_bins=False)
        expected = 0.5 * (result['edges'][:-1] + result['edges'][1:])
        np.testing.assert_allclose(result['bin_centers'], expected)

    def test_total_count_matches_input_size(self):
        data = np.arange(1, 101, dtype=float)
        result = binned_histogram(data, nbins=10, log_bins=False)
        assert result['counts'].sum() == 100


class TestBootstrapCI:

    def test_returns_two_percentiles(self):
        ci = bootstrap_ci([1, 2, 3, 4, 5], n_boot=100)
        assert ci.shape == (2,)
        assert ci[0] < ci[1]

    def test_ci_contains_true_mean(self):
        np.random.seed(42)
        data = np.random.normal(loc=5.0, scale=1.0, size=300)
        ci = bootstrap_ci(data, statistic='mean', n_boot=500)
        assert ci[0] < 5.0 < ci[1]

    def test_ci_contains_true_median(self):
        np.random.seed(42)
        data = np.random.normal(loc=5.0, scale=1.0, size=300)
        ci = bootstrap_ci(data, statistic='median', n_boot=500)
        assert ci[0] < 5.0 < ci[1]

    def test_unknown_statistic_raises_value_error(self):
        with pytest.raises(ValueError):
            bootstrap_ci([1, 2, 3], statistic='variance', n_boot=10)
