"""FLOWR path configuration.

Fill in the absolute paths to your FLOWR installation here so the app
can launch jobs without requiring environment variables. These values
are only used when the corresponding keyword argument or environment
variable is missing.
"""
from __future__ import annotations

FLOWR_SETTINGS = {
    # Absolute path to the FLOWR repository (the folder that contains
    # the `flowr` python package and `setup.cfg`).
    "repo_path": "",
    # Python interpreter that has FLOWR + RDKit dependencies installed.
    "python_path": "",
    # Default checkpoint to sample from (.ckpt file).
    "checkpoint_path": "",
    # Optional paths below – leave blank to skip.
    "data_path": "",
    "dataset": "",
    "arch": "pocket",
    "pocket_type": "holo",
    "batch_cost": 100,
    "num_workers": 8,
    "gpus": 1,
}
