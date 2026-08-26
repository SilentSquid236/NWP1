"""
Central path + hyperparameter configuration.

Everything is driven by environment variables so the exact same code runs
unchanged in WSL and on the Xeon server. Set NWP_DATA_ROOT per machine:

    WSL:    export NWP_DATA_ROOT=/mnt/c/Users/Epier/Desktop/NWP/NWP_Deployment_Package/data
    Server: export NWP_DATA_ROOT=/data5/pierce/Data5/NWP/data

If unset, it defaults to ./data relative to this file.
"""

import os
from pathlib import Path

# --- Paths -----------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent

DATA_ROOT = Path(os.environ.get("NWP_DATA_ROOT", PROJECT_ROOT / "data"))
TENSOR_DIR = DATA_ROOT / "tensors_3d"
CHECKPOINT_DIR = Path(os.environ.get("NWP_CHECKPOINT_DIR", DATA_ROOT / "checkpoints"))

# --- Atmospheric state definition -----------------------------------------

# Channel order is fixed and must match the ingestion script's output.
CHANNELS = ["TMP", "RH", "UGRD", "VGRD", "HGT"]
N_CHANNELS = len(CHANNELS)

# Pressure levels in hPa, ordered surface -> top. Spacing is deliberately
# uneven: 25 hPa through the boundary layer and lower troposphere where
# gradients are sharp, coarsening aloft where fields are smoother.
# HRRR's prs product carries every 25 hPa, so any subset of these is free.
PRESSURE_LEVELS = [
    1000, 975, 950, 925, 900, 875, 850, 825, 800, 750,
     700, 650, 600, 550, 500, 450, 400, 300, 250, 200,
]
N_LEVELS = len(PRESSURE_LEVELS)

# Regional domain -- approximates the HRRR Northeast sector.
DOMAIN = {
    "name": "northeast",
    "lat_min": 37.0, "lat_max": 47.5,
    "lon_min": -82.0, "lon_max": -66.0,
}

# Index of temperature, used by the control-variable loss weighting.
TEMP_CHANNEL_IDX = 0

# --- Training defaults ------------------------------------------------------

BATCH_SIZE = int(os.environ.get("NWP_BATCH_SIZE", 1))
LEARNING_RATE = float(os.environ.get("NWP_LR", 1e-4))
EPOCHS = int(os.environ.get("NWP_EPOCHS", 50))
# Loss weights, aligned to CHANNELS. Temperature is the strict control
# variable, so it carries a heavy multiplier.
LOSS_WEIGHTS = [5.0, 1.0, 1.0, 1.0, 1.0]

# Thread/worker limits live in resources.py, which must be imported before
# torch. Default policy: never use more than 50% of the server's cores.


def ensure_dirs():
    """Create output directories if they don't exist."""
    TENSOR_DIR.mkdir(parents=True, exist_ok=True)
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)


def describe() -> str:
    return (
        f"  DATA_ROOT      : {DATA_ROOT}\n"
        f"  TENSOR_DIR     : {TENSOR_DIR}\n"
        f"  CHECKPOINT_DIR : {CHECKPOINT_DIR}\n"
        f"  state          : {N_CHANNELS} vars x {N_LEVELS} levels {CHANNELS}\n"
        f"  levels (hPa)   : {PRESSURE_LEVELS[0]} -> {PRESSURE_LEVELS[-1]}\n"
        f"  domain         : {DOMAIN['name']} "
        f"{DOMAIN['lat_min']}-{DOMAIN['lat_max']}N "
        f"{DOMAIN['lon_min']}-{DOMAIN['lon_max']}E\n"
        f"  batch/lr/epochs: {BATCH_SIZE} / {LEARNING_RATE} / {EPOCHS}"
    )


if __name__ == "__main__":
    print("NWP configuration:")
    print(describe())
