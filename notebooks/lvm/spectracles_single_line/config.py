"""Shared configuration for the single emission line tutorial notebooks.

Change LINE below to switch which emission line both notebooks use.
Look at the NOTE(s) in different places below for other things you might want to consider changing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from lvm_tools import DataConfig, FitDataBuilder, LVMTile, LVMTileCollection
from lvm_tools.fit_data.filtering import BAD_FLUX_THRESHOLD

# ======================================================================
# Pick which line to fit / plot — change this one variable
# ======================================================================

# NOTE: Choose what line to fit and plot here.
LINE = "nii"  # "halpha", "nii", or "oiii"

# ======================================================================
# Per-line settings
# ======================================================================


@dataclass(frozen=True)
class LineConfig:
    name: str
    lambda_centre: float
    label: str
    norm_scale: float
    model_file: str
    v_syst_init: float = 25.0  # initial guess for bulk systemic velocity [km/s]


# NOTE: ADD NEW LINES HERE IF YOU WANT TO FIT / PLOT THEM IN THE NOTEBOOKS
LINES = {
    "halpha": LineConfig("halpha", 6562.85, r"H$\alpha$", 1e-12, "halpha.model"),
    "nii": LineConfig("nii", 6583.45, r"[NII] 6583", 0.3e-12, "nii.model"),
    "oiii": LineConfig("oiii", 5006.84, r"[OIII] 5007", 0.01e-12, "oiii.model"),
}

CFG = LINES[LINE]

# ======================================================================
# Paths and data files — change these for a different dataset
# ======================================================================

# NOTE: Change this stuff to change what data you load.
# The example here is 3 tiles covering the flame nebula.
DATA_DIR = Path("../../../data/work/lvm/1.2.0")
TILE_IDS = [1028921, 1028922, 1028892]
DRP_FILES = [list(DATA_DIR.glob(f"*/{tile_id}/60*/*.fits"))[0] for tile_id in TILE_IDS]
MODEL_DIR = Path("models")
MODEL_PATH = MODEL_DIR / CFG.model_file

# ======================================================================
# Data loading
# ======================================================================

# NOTE: LINE_WINDOW_HALF controls how much spectrum to include around the line centre.
# You might want to change this to avoid close lines in some cases, for example.
LINE_WINDOW_HALF = 8.0  # Angstroms either side of line centre

# NOTE: You need more modes if your domain gets really big (in terms of sky footprint)
# 401 x 401 is sufficient for problems up to the size of 19 contiguous tiles, the example
# here is only 3, so it's way overkill. You could go down to 201 x 201 for a speed-up.
# NOTE: Both numbers MUST BE ODD, for annoying technical reasons.
N_MODES = (401, 401)


def load_flame(cfg: LineConfig = CFG):
    """Load LVM data for the configured line."""
    λ_range = (
        cfg.lambda_centre - LINE_WINDOW_HALF,
        cfg.lambda_centre + LINE_WINDOW_HALF,
    )

    for p in DRP_FILES:
        assert p.exists(), f"DRP file not found: {p}"
    tiles = LVMTileCollection.from_tiles([LVMTile.from_file(p) for p in DRP_FILES])

    data_config = DataConfig.from_tiles(
        tiles,
        λ_range,
        normalise_F_scale=cfg.norm_scale,
        normalise_F_offset=0.0,
        F_range=(BAD_FLUX_THRESHOLD, 1e-12),
        F_bad_strategy="spaxel",
    )

    builder = FitDataBuilder(tiles, data_config)
    return builder.build()
