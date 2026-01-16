# SBProbMol3

SBProbMol3 is a Flask web app for structure-based drug design workflows:

- Browse/select a protein target (curated list or upload your own PDB)
- Select a reference ligand (upload SDF/MOL/PDB)
- Generate candidate ligands via the Flowr/“DiffSBDD” generation pipeline (runs as a subprocess)
- Inspect results in-browser with an NGL viewer and download artifacts
- (Optional) run quick MD with FairChem (ASE + FAIRChemCalculator) and view the trajectory as PDB
- (Optional) annotate generated ligands with SA score/QED (RDKit) and docking scores (AutoDock-GPU toolchain)

This repo contains the web UI + orchestration layer; it expects external model code/checkpoints for generation.

## Quickstart

### 1) Create the environment

The simplest option is to use the provided conda environment file:

```bash
cd SBProbMol3_1
conda env create -f environment.yaml
conda activate SB
```

Notes:

- GPU support is optional but recommended for ligand generation.
- The MD feature additionally requires `ase` and `fairchem` (see “MD module” below).

### 2) Configure the generation backend (Flowr)

The server starts generation jobs by invoking a Python interpreter that has the Flowr code available.
Set these environment variables before running the app:

- `DIFFSBDD_REPO`: path to the Flowr repo (contains the `flowr` Python package)
- `DIFFSBDD_PYTHON`: path to the Python executable in the Flowr environment
- `DIFFSBDD_CHECKPOINT`: path to a model checkpoint file

Example:

```bash
export DIFFSBDD_REPO=/abs/path/to/flowr_root
export DIFFSBDD_PYTHON=/abs/path/to/conda/env/bin/python
export DIFFSBDD_CHECKPOINT=/abs/path/to/flowr_root_v2.ckpt
```

If you do not set them, the app falls back to local developer defaults defined in [utils/diffsbdd.py](utils/diffsbdd.py).

### 3) Run the web app

```bash
python app.py
```

Then open `http://localhost:5000`.

## How it works

### Web workflow

The UI guides you through four steps:

1. **Target selection**: pick from curated targets or upload a `.pdb`
2. **Pocket/reference**: select/upload a reference ligand (used as the “ligand_file” input)
3. **Sampling parameters**: choose `n_samples`, `timesteps`, etc.
4. **Generation & results**: monitor progress, list molecules, view them in NGL, download the SDF

### Outputs and storage

Uploaded files and generated artifacts are stored under [database/](database/):

- `database/proteins/`: uploaded PDBs
- `database/ligands/`: uploaded ligand files
- `database/diffsbdd/<job_id>/`: generation outputs (SDF, logs, derived PDBs)
- `database/MD/` and `database/MD/temp/`: MD inputs and cached trajectory conversions

Make sure this folder is writable by the process running Flask/Gunicorn.

## Key routes (API)

- `POST /upload_pdb` – upload a protein PDB
- `POST /upload_ligand` – upload a ligand file (PDB/MOL/SDF)
- `POST /save_sampling_config` – store sampling settings in the session
- `POST /diffsbdd/start` – start a generation job
- `GET /diffsbdd/status/<job_id>` – poll job status
- `GET /diffsbdd/result/<job_id>/sdf` – download job SDF
- `GET /diffsbdd/result/<job_id>/list` – list molecules in the SDF
- `GET /diffsbdd/result/<job_id>/pdb/<index>` – convert and serve a single molecule as PDB
- `GET /md` – MD viewer page
- `POST /md/upload` – upload MD input files
- `POST /md/run` – run FairChem MD (SDF only), convert to PDB and return a URL
- `POST /md/import_generated` – split a generated job SDF into per-ligand SDFs for MD

## Optional modules

### MD module (FairChem + ASE)

The MD page (`/md`) can run a short Langevin dynamics simulation using FairChem models.
This requires:

- `ase` (and the `ase` CLI, for `ase convert`)
- `fairchem` (imported as `fairchem.core`)

If your environment doesn’t include these, you can add them (example):

```bash
pip install ase fairchem-core
```

### Docking score annotation (AutoDock-GPU toolchain)

After generation completes, the app can optionally annotate the output SDF with docking scores.
This is best-effort and requires external binaries on `PATH`, including:

- `mk_prepare_receptor.py`
- `obabel`
- `autogrid4`
- `autodock_gpu_128wi`

If unavailable, the app will still generate ligands; it just won’t add `VINA_SCORE` properties.

## Deployment

For production, run behind a reverse proxy:

```bash
gunicorn -w 2 -b 0.0.0.0:5000 app:app
```

Set `SESSION_SECRET` to a strong value in production:

```bash
export SESSION_SECRET='change-me'
```

## Project structure

```text
app.py                # Flask app entrypoint
routes.py             # HTTP routes (UI + API)
data/                 # curated proteins + helpers
templates/            # Jinja2 templates
static/               # JS/CSS/assets (NGL-based viewer UI)
utils/                # Flowr runner, MD script, Vina/AutoDock-GPU helpers
database/             # uploads + job artifacts (runtime data)
```

## License

Add a license file if you plan to publish/distribute this project.
