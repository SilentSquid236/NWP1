import os
import re
import glob
import torch
from torch.utils.data import Dataset

# Matches live_hrrr_f00.pt, live_hrrr_f6.pt, live_hrrr_f018.pt ...
_FXX = re.compile(r"_f(\d+)\.pt$")


def _forecast_hour(path: str):
    """Extract the integer forecast hour from a tensor filename, or None."""
    m = _FXX.search(os.path.basename(path))
    return int(m.group(1)) if m else None


class AutoregressiveDataset3D(Dataset):
    """
    Maps the atmospheric state from Time T to Time T+1.

    Expects one or more directories of sequential hourly .pt tensors named
    live_hrrr_f<HH>.pt. Only consecutive hours (t+1 - t == 1) become pairs,
    so gaps in a run are skipped rather than silently learned as a jump.
    """

    def __init__(self, run_folders, verbose=True):
        self.valid_pairs = []

        for folder in run_folders:
            files = glob.glob(os.path.join(str(folder), "live_hrrr_f*.pt"))

            # Sort numerically by forecast hour, NOT lexicographically --
            # a plain sort puts f10 before f2 whenever hours aren't padded.
            hours = [(h, p) for p in files if (h := _forecast_hour(p)) is not None]
            hours.sort(key=lambda x: x[0])

            for (h_t, p_t), (h_t1, p_t1) in zip(hours, hours[1:]):
                if h_t1 - h_t == 1:
                    self.valid_pairs.append((p_t, p_t1))

        if verbose:
            print(f"Autoregressive dataset: {len(self.valid_pairs)} hourly transitions "
                  f"from {len(run_folders)} run folder(s).")
            if not self.valid_pairs:
                print("  WARNING: no consecutive-hour pairs found. Check that the "
                      "ingestion step has written live_hrrr_f*.pt files.")

    def __len__(self):
        return len(self.valid_pairs)

    def __getitem__(self, idx):
        path_t, path_t1 = self.valid_pairs[idx]

        # Shape: [C, L, Y, X] -- channels, vertical levels, grid
        data_t = torch.load(path_t, weights_only=True)["hrrr_features"]
        data_t1 = torch.load(path_t1, weights_only=True)["hrrr_features"]

        return data_t.float(), data_t1.float()
