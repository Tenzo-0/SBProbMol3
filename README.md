# MoleAG

A lightweight web app for structure-conditioned molecular generation and 3D visualization. It integrates DiffSBDD for ligand generation, uses RDKit to convert SDF to PDB, and visualizes proteins/ligands with NGL.

## Highlights
- Target selection: choose curated proteins or upload your own PDB.
- Pocket definition: simple "Reference Ligand / Pocket" field (e.g., `A:300`).
- Sampling: configure DiffSBDD parameters (n_samples, timesteps, etc.).
- Job lifecycle: start, live progress, logs, and result artifacts.
- Results: list generated molecules; click to visualize with the selected protein.
- Visualization: NGL viewer with cartoon/stick/surface, quality and opacity controls.

## Requirements
- OS: Ubuntu 20.04/22.04 recommended (Windows for local dev is OK; production on Linux preferred).
- NVIDIA GPU with CUDA 11.8 (optional, recommended). CPU works but slower.
- Conda (Miniconda/Anaconda) for DiffSBDD environment.
- Python 3.10+ for the web app (can reuse DiffSBDD env or a separate venv).

### DiffSBDD environment (conda)
Use the environment.yaml from DiffSBDD (example pins):
- python=3.10.4, pytorch=2.0.1 (CUDA 11.8), cudatoolkit=11.8
- pytorch-scatter=2.1.2, pytorch-lightning=1.8.4, torchmetrics=1.4.2
- rdkit=2022.03.2, openbabel=3.1.1, biopython, scipy, numpy, pandas, networkx, seaborn, tqdm, protobuf 3.20.*
- channels: pyg, pytorch, nvidia, anaconda, conda-forge, defaults

## Setup

1) Clone
```bash
git clone https://github.com/arneschneuing/DiffSBDD.git
git clone https://github.com/Tenzo-0/MoleAG.git
```

2) Create the DiffSBDD env
```bash
cd DiffSBDD
conda env create -f DiffSBDD -n diffsbdd
conda activate diffsbdd
```

5) Download model checkpoints
- You need one or more DiffSBDD model checkpoints (e.g., conditional and/or unconditional).
```bash
wget -P checkpoints/ https://zenodo.org/record/8183747/files/crossdocked_fullatom_cond.ckpt
```

3) Install app dependencies (Flask, etc.)
```bash
cd
cd MoleAG
pip install -U pip setuptools wheel
pip install -e .  # if your pyproject lists deps
```

4) Configure environment variables
- DIFFSBDD_REPO: absolute path to DiffSBDD repo
- DIFFSBDD_PYTHON: full path to the python inside the `diffsbdd` conda env
- DIFFSBDD_CHECKPOINT: path to your DiffSBDD model checkpoint file
```bash
#example
export DIFFSBDD_REPO=/home/user/DiffSBDD
export DIFFSBDD_PYTHON=/home/user/anaconda3/envs/diffsbdd/bin/python
export DIFFSBDD_CHECKPOINT=/home/user/DiffSBDD/checkpoints/cond.ckpt
```

6) Run (development)
```bash
flask run #if you want to run it in debug mode, add [--debug]
```

## Usage
1. Pick a protein (or upload .pdb). The viewer loads the structure in NGL.
2. Define the pocket with the "Reference Ligand / Pocket" input (e.g., chain:resid).
3. Set sampling parameters and click Run Generation.
4. Watch progress. When done, Step 4 lists generated molecules; click to visualize.
   - The protein remains; the selected ligand replaces any previously shown ligand.
   - Download the full SDF from the Results panel.

## API (selected)
- POST `/upload_pdb` — upload a PDB file.
- POST `/save_sampling_config` — save UI parameters before starting a job.
- POST `/diffsbdd/start` — start a DiffSBDD job.
- GET `/diffsbdd/status/<job_id>` — poll job status.
- GET `/diffsbdd/result/<job_id>/sdf` — download generated SDF.
- GET `/diffsbdd/result/<job_id>/list` — list molecules in SDF.
- GET `/diffsbdd/result/<job_id>/pdb/<index>` — serve PDB for a molecule index.
- GET `/database/<filename>` — serve user-uploaded PDBs.

## Project structure
```
app.py
routes.py
utils/
  diffsbdd.py        # DiffSBDD runner (jobs, status, SDF counting, RDKit conversions)
static/
  css/style.css
  js/app.js          # NGL viewer + workflow and results UI
  img/logo.png
templates/
  base.html
  index.html
database/            # user uploads and job artifacts
```

## Deployment (production)
- Use Gunicorn (or uvicorn+ASGI wrapper) behind NGINX.
- Point NGINX to serve `/static` directly and reverse-proxy to the Flask app.
- Ensure `database/` is writable by the app service user.
- Set the DIFFSBDD_* env vars in your systemd unit or process manager.

## Troubleshooting
- RDKit not found: make sure `DIFFSBDD_PYTHON` points to the conda env that has RDKit.
- PyTorch/torch-scatter mismatch: install versions matching your CUDA (11.8 suggested).
- No results listed: check `database/diffsbdd/<job_id>/run.log` and the configured checkpoint path.
- 3D viewer blank: verify network access to `/database/<file>.pdb` or fallback to `rcsb://<pdbid>`.

## Acknowledgements
- DiffSBDD: Diffusion-based SBDD framework for conditional ligand generation.
- NGL: WebGL-based molecular visualization library.

## License
Add your license here (e.g., MIT).
# SBProbMol3
