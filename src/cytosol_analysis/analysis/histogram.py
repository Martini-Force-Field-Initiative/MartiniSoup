import numpy as np

class HistogramAnalysis:
    """
    Compute linear or log-binned histograms of residence times.
    """

    def __init__(self, durations):
        self.durations = np.asarray(durations, dtype=float)

    def compute(self, nbins=30, log_bins=False):
        """
        Parameters
        ----------
        nbins : int
        log_bins : bool

        Returns
        -------
        dict with keys ["bin_edges", "bin_centers", "counts", "pdf"]
        """
        data = self.durations[self.durations > 0]

        if len(data) == 0:
            return {
                "bin_edges": np.array([]),
                "bin_centers": np.array([]),
                "counts": np.array([]),
                "pdf": np.array([]),
            }

        if log_bins:
            edges = np.logspace(
                np.log10(data.min()),
                np.log10(data.max()),
                nbins + 1
            )
        else:
            edges = np.linspace(data.min(), data.max(), nbins + 1)

        counts, edges = np.histogram(data, bins=edges)

        # PDF normalized
        widths = np.diff(edges)
        pdf = counts / (counts.sum() * widths)

        centers = 0.5 * (edges[:-1] + edges[-1:])

        return {
            "bin_edges": edges,
            "bin_centers": centers,
            "counts": counts,
            "pdf": pdf,
        }
