import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from toad.shifts import ShiftsMethod
from toad.shifts import ASDETECT

class FilteredASDETECT(ShiftsMethod):
    def __init__(
        self,
        lmin: int = 5,
        lmax: int = None,
        segmentation: str = "two_sided",
        window_indices: int = 30,
        min_rel_change: float = 0.25,
        min_abs_change: float = 0.0
    ):
        # Store ASDETECT parameters for cloning
        self.lmin = lmin
        self.lmax = lmax
        self.segmentation = segmentation

        # Parameters for the 2nd-tier filtering
        # Assuming annual data, window_indices corresponds to years
        self.window_indices = window_indices
        self.min_rel_change = min_rel_change
        self.min_abs_change = min_abs_change

    @property
    def asdetect(self):
        """Lazy-load ASDETECT instance (not stored in __dict__ to avoid clone issues).

        Creates the ASDETECT instance on-demand instead of storing it.
        This prevents TOAD's _clone_method from trying to pass 'asdetect'
        as an __init__ parameter during parallelization.
        """
        if not hasattr(self, '_asdetect_cache'):
            self._asdetect_cache = ASDETECT(lmin=self.lmin, lmax=self.lmax, segmentation=self.segmentation)
        return self._asdetect_cache

    def fit_predict(
        self,
        values_1d: np.ndarray,
        times_1d: np.ndarray,
    ) -> np.ndarray:
        # 1. Get raw shift scores using the base ASDETECT (3 MAD criteria)
        raw_shifts = self.asdetect.fit_predict(values_1d, times_1d)
        filtered_shifts = np.zeros_like(raw_shifts)

        # 2. Find indices where ASDETECT detected a potential shift (non-zero values)
        candidate_indices = np.where(raw_shifts != 0)[0]

        # 3. Apply the 2nd-tier physical/relative threshold filtering
        for idx in candidate_indices:
            start_idx = max(0, idx - self.window_indices)
            end_idx = min(len(values_1d), idx + self.window_indices)

            # Skip if there is not enough data before or after to compare
            if idx == start_idx or idx == end_idx:
                continue

            # Calculate the mean values before and after the detected shift
            mean_before = np.mean(values_1d[start_idx:idx])
            mean_after = np.mean(values_1d[idx:end_idx])

            # Calculate absolute and relative changes
            abs_change = np.abs(mean_after - mean_before)
            
            # Prevent division by zero
            if mean_before != 0:
                rel_change = abs_change / np.abs(mean_before)
            else:
                rel_change = 0.0

            # 4. Filter condition: Must meet BOTH absolute and relative minimum changes
            if (abs_change >= self.min_abs_change) and (rel_change >= self.min_rel_change):
                # Retain the original shift score if it passes the filter
                filtered_shifts[idx] = raw_shifts[idx]

        return filtered_shifts