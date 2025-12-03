import os
import uuid
import time
import threading
import subprocess
import logging
from dataclasses import dataclass, field
from typing import Dict, Optional, List

from urllib.request import urlretrieve


logger = logging.getLogger(__name__)


@dataclass
class DiffSBDDJob:
    id: str
    state: str = "queued"  # queued | running | success | failed
    message: str = ""
    started_at: float = field(default_factory=time.time)
    finished_at: Optional[float] = None
    n_samples: int = 0
    outfile: str = ""
    log_file: str = ""
    proc: Optional[subprocess.Popen] = None
    receptor_pdb: Optional[str] = None


@dataclass
class MoleculeScore:
    index: int
    title: str
    vina_score: Optional[float] = None
    qed: Optional[float] = None
    sa: Optional[float] = None


class DiffSBDDRunner:
    """
    Thin wrapper to run DiffSBDD's generate_ligands.py as a subprocess.
    Requires environment variables:
      - DIFFSBDD_PYTHON: path to Python interpreter in the DiffSBDD conda env (python.exe)
      - DIFFSBDD_REPO: local path to the DiffSBDD repo root
      - DIFFSBDD_CHECKPOINT: path to a .ckpt model file
    """

    def __init__(self,
                 uploads_dir: str,
                 repo_path: Optional[str] = None,
                 python_path: Optional[str] = None,
                 checkpoint_path: Optional[str] = None):
        self.uploads_dir = os.path.abspath(uploads_dir)
        self.repo_path = repo_path or os.environ.get('DIFFSBDD_REPO', '')
        self.python_path = python_path or os.environ.get('DIFFSBDD_PYTHON', '')
        self.checkpoint_path = checkpoint_path or os.environ.get('DIFFSBDD_CHECKPOINT', '')
        self.jobs: Dict[str, DiffSBDDJob] = {}

    def is_configured(self) -> bool:
        return all([
            self.repo_path and os.path.isdir(self.repo_path),
            self.python_path and os.path.isfile(self.python_path),
            self.checkpoint_path and os.path.isfile(self.checkpoint_path),
        ])

    def ensure_pdb_local(self, pdb_id: Optional[str], filename: Optional[str]) -> str:
        """Return absolute path to a local PDB file. If not present, try to download by pdb_id."""
        # 1) If an uploaded filename is provided and exists, use it
        if filename:
            local = os.path.join(self.uploads_dir, filename)
            if os.path.isfile(local):
                return local
        # 2) If a file with the PDB ID exists in uploads, use it
        if pdb_id:
            candidate = os.path.join(self.uploads_dir, f"{pdb_id}.pdb")
            if os.path.isfile(candidate):
                return candidate
        # 3) Try to download from RCSB using pdb_id
        if pdb_id:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            dest = os.path.join(self.uploads_dir, f"{pdb_id}.pdb")
            try:
                os.makedirs(self.uploads_dir, exist_ok=True)
                logger.info("Downloading PDB %s -> %s", url, dest)
                urlretrieve(url, dest)
                if os.path.isfile(dest):
                    return dest
            except Exception as e:
                logger.error("Failed to download PDB %s: %s", pdb_id, e)
                raise
        raise FileNotFoundError("No local PDB file found and download failed. Please upload a PDB file.")

    def start_job(self,
                  pdb_id: Optional[str],
                  filename: Optional[str],
                  ref_ligand: str,
                  n_samples: int,
                  num_nodes_lig: Optional[int] = None,
                  timesteps: Optional[int] = None,
                  resamplings: Optional[int] = None,
                  jump_length: Optional[int] = None,
                  keep_all_fragments: bool = False,
                  sanitize: bool = False,
                  relax: bool = False) -> DiffSBDDJob:
        if not self.is_configured():
            raise RuntimeError("DiffSBDD not configured. Set DIFFSBDD_REPO, DIFFSBDD_PYTHON, DIFFSBDD_CHECKPOINT.")

        pdb_path = self.ensure_pdb_local(pdb_id, filename)
        job_id = str(uuid.uuid4())
        outdir = os.path.join(self.uploads_dir, 'diffsbdd', job_id)
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, 'output.sdf')
        log_file = os.path.join(outdir, 'run.log')

        cmd = [
            self.python_path,
            os.path.join(self.repo_path, 'generate_ligands.py'),
            self.checkpoint_path,
            '--pdbfile', pdb_path,
            '--outfile', outfile,
            '--ref_ligand', ref_ligand,
            '--n_samples', str(int(n_samples)),
        ]
        if num_nodes_lig is not None:
            cmd += ['--num_nodes_lig', str(int(num_nodes_lig))]
        if timesteps is not None:
            cmd += ['--timesteps', str(int(timesteps))]
        if resamplings:
            cmd += ['--resamplings', str(int(resamplings))]
        if jump_length:
            cmd += ['--jump_length', str(int(jump_length))]
        if keep_all_fragments:
            cmd += ['--all_frags']
        if sanitize:
            cmd += ['--sanitize']
        if relax:
            cmd += ['--relax']

        job = DiffSBDDJob(id=job_id, state='running', message='Starting DiffSBDD...', n_samples=n_samples,
                   outfile=outfile, log_file=log_file, receptor_pdb=pdb_path)
        self.jobs[job_id] = job

        def _run():
            try:
                with open(log_file, 'w', encoding='utf-8', errors='ignore') as lf:
                    lf.write('Command: ' + ' '.join([f'"{c}"' if ' ' in str(c) else str(c) for c in cmd]) + '\n')
                    lf.flush()
                    proc = subprocess.Popen(cmd, cwd=self.repo_path, stdout=lf, stderr=subprocess.STDOUT, shell=False)
                    job.proc = proc
                    ret = proc.wait()
                job.finished_at = time.time()
                if ret == 0 and os.path.isfile(outfile):
                    job.state = 'success'
                    job.message = 'Completed'
                else:
                    job.state = 'failed'
                    job.message = f'Process exited with code {ret}'
            except Exception as e:
                job.finished_at = time.time()
                job.state = 'failed'
                job.message = f'Error: {e}'

        t = threading.Thread(target=_run, daemon=True)
        t.start()
        return job

    def count_generated_molecules(self, job: DiffSBDDJob) -> int:
        try:
            if not os.path.isfile(job.outfile):
                return 0
            count = 0
            with open(job.outfile, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.strip() == '$$$$':
                        count += 1
            return count
        except Exception:
            return 0

    def get_status(self, job_id: str) -> Dict:
        job = self.jobs.get(job_id)
        if not job:
            return {'exists': False}
        # Estimate progress based on molecules count relative to target
        mols = self.count_generated_molecules(job)
        progress = 5
        if job.n_samples > 0:
            progress = min(95, int(5 + 90 * (mols / max(1, job.n_samples))))
        if job.state == 'success':
            progress = 100
        if job.state == 'failed':
            progress = min(progress, 99)
        # Read last log lines (tail)
        tail = ''
        try:
            if os.path.isfile(job.log_file):
                with open(job.log_file, 'r', encoding='utf-8', errors='ignore') as f:
                    lines = f.readlines()
                    tail = ''.join(lines[-20:])
        except Exception:
            pass
        # Build a URL for download relative to /uploads route
        rel_outfile = os.path.relpath(job.outfile, self.uploads_dir).replace('\\', '/')
        return {
            'exists': True,
            'job_id': job.id,
            'state': job.state,
            'message': job.message,
            'progress': progress,
            'molecules': mols,
            'outfile_url': f"/uploads/{rel_outfile}" if os.path.isfile(job.outfile) else None,
            'log_tail': tail,
        }

    def get_job(self, job_id: str) -> Optional[DiffSBDDJob]:
        return self.jobs.get(job_id)

    def get_job_dir(self, job_id: str) -> str:
        return os.path.join(self.uploads_dir, 'diffsbdd', job_id)

    # --- AutoDock / PDBQT helpers -------------------------------------------------

    @staticmethod
    def _is_float(val: str) -> bool:
        try:
            float(val)
            return True
        except Exception:
            return False

    def _prepare_receptor_pdbqt(self, job: DiffSBDDJob) -> Optional[str]:
        """Convert receptor PDB to PDBQT using OpenBabel if needed.

        Requires the `obabel` executable in PATH. This function is best-effort:
        if anything fails it returns None and VINAScore will be skipped.
        """
        if not job.receptor_pdb or not os.path.isfile(job.receptor_pdb):
            return None
        workdir = self.get_job_dir(job.id)
        os.makedirs(workdir, exist_ok=True)
        receptor_pdbqt = os.path.join(workdir, 'receptor.pdbqt')
        if os.path.isfile(receptor_pdbqt):
            return receptor_pdbqt
        cmd = ['obabel', '-ipdb', job.receptor_pdb, '-xp', '-opdbqt', '-O', receptor_pdbqt]
        try:
            ret = subprocess.run(cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, check=False)
            if ret.returncode == 0 and os.path.isfile(receptor_pdbqt):
                return receptor_pdbqt
            logger.error('Failed to prepare receptor PDBQT (rc=%s): %s', ret.returncode, ret.stderr)
        except Exception as e:
            logger.exception('Error preparing receptor PDBQT: %s', e)
        return None

    def _prepare_ligand_pdbqt(self, mol, index: int, job: DiffSBDDJob) -> Optional[str]:
        """Write RDKit mol to PDB and convert to PDBQT using OpenBabel.

        Returns absolute path to ligand PDBQT or None on failure.
        """
        try:
            from rdkit import Chem
        except ImportError:
            logger.error('RDKit not available; cannot export ligand PDB for docking.')
            return None

        workdir = self.get_job_dir(job.id)
        os.makedirs(workdir, exist_ok=True)
        pdb_path = os.path.join(workdir, f'lig_{index+1}.pdb')
        lig_pdbqt = os.path.join(workdir, f'lig_{index+1}.pdbqt')
        if os.path.isfile(lig_pdbqt):
            return lig_pdbqt
        try:
            Chem.MolToPDBFile(mol, pdb_path)
        except Exception as e:
            logger.exception('Failed to write ligand PDB: %s', e)
            return None

        cmd = ['obabel', '-ipdb', pdb_path, '-xp', '-opdbqt', '-O', lig_pdbqt]
        try:
            ret = subprocess.run(cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, check=False)
            if ret.returncode == 0 and os.path.isfile(lig_pdbqt):
                return lig_pdbqt
            logger.error('Failed to prepare ligand PDBQT (rc=%s): %s', ret.returncode, ret.stderr)
        except Exception as e:
            logger.exception('Error preparing ligand PDBQT: %s', e)
        return None

    def _run_autodock_gpu(self, receptor_pdbqt: str, ligand_pdbqt: str, job: DiffSBDDJob, index: int) -> Optional[float]:
        """Run AutoDock-GPU for a single ligand and parse the best affinity.

        Requires AUTODOCK_GPU_BIN env var pointing to the executable. This
        uses a generic invocation and attempts to parse a Vina-like RESULT
        line from the log file. If anything fails, returns None.
        """
        autodock_bin = os.environ.get('AUTODOCK_GPU_BIN')
        if not autodock_bin or not os.path.isfile(autodock_bin):
            logger.error('AUTODOCK_GPU_BIN is not configured or not found; skipping VINAScore.')
            return None

        workdir = self.get_job_dir(job.id)
        os.makedirs(workdir, exist_ok=True)
        out_pdbqt = os.path.join(workdir, f'out_{index+1}.pdbqt')
        log_path = os.path.join(workdir, f'log_{index+1}.txt')

        # NOTE: Adjust arguments according to your AutoDock-GPU build/CLI.
        cmd = [
            autodock_bin,
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--out', out_pdbqt,
            '--log', log_path,
        ]

        try:
            ret = subprocess.run(cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, check=False)
            if ret.returncode != 0:
                logger.error('AutoDock-GPU failed (rc=%s): %s', ret.returncode, ret.stderr)
                return None
        except Exception as e:
            logger.exception('AutoDock-GPU execution error: %s', e)
            return None

        if not os.path.isfile(log_path):
            return None

        best_affinity: Optional[float] = None
        try:
            with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    # Try to parse Vina-style result lines
                    if 'VINA RESULT' in line or 'Affinity' in line:
                        parts = line.replace(',', ' ').split()
                        vals = [p for p in parts if self._is_float(p)]
                        if vals:
                            best_affinity = float(vals[0])
                            break
        except Exception as e:
            logger.exception('Failed to parse AutoDock-GPU log: %s', e)
            return None
        return best_affinity

    def convert_sdf_to_pdb(self, job_id: str, select: str = 'first') -> Optional[str]:
        """Convert a job's SDF to PDB using RDKit in the DiffSBDD env. Returns output path or None."""
        job = self.get_job(job_id)
        if not job:
            return None

    def list_sdf_molecules(self, job_id: str):
        """Return a list of dicts {index, title} for molecules in the job's SDF using RDKit."""
        job = self.get_job(job_id)
        if not job or not os.path.isfile(job.outfile):
            return []
        script = (
            "import sys,json; from rdkit import Chem; "
            "mols=[m for m in Chem.SDMolSupplier(sys.argv[1], removeHs=False) if m]; "
            "out=[{'index':i,'title':(m.GetProp('_Name') if m.HasProp('_Name') else f'mol_{i+1}')} for i,m in enumerate(mols)]; "
            "print(json.dumps(out))"
        )
        try:
            ret = subprocess.run([self.python_path, '-c', script, job.outfile], cwd=self.repo_path,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
            if ret.returncode == 0 and ret.stdout:
                import json as _json
                return _json.loads(ret.stdout.strip())
            logger.error('List SDF failed: rc=%s err=%s', ret.returncode, ret.stderr)
            return []
        except Exception as e:
            logger.exception('List SDF error: %s', e)
            return []

    def score_sdf_molecules(self, job_id: str) -> List[MoleculeScore]:
        """Return per-molecule VINAScore (placeholder), QED and SA scores using RDKit.

        VINAScore here is a simple placeholder based on heavy atom count so that
        the UI always has a numeric column even if an external docking engine is
        not configured. Replace with a proper docking call if available.
        """
        job = self.get_job(job_id)
        if not job or not os.path.isfile(job.outfile):
            return []

        try:
            from rdkit import Chem
            from rdkit.Chem import QED
            try:
                # SA score is in rdkit Contrib (sascorer)
                from rdkit.Chem import rdMolDescriptors
                have_sa = hasattr(rdMolDescriptors, 'CalcSyntheticAccessibilityScore')
            except Exception:  # pragma: no cover - contrib may be missing
                rdMolDescriptors = None
                have_sa = False
        except ImportError:
            logger.error('RDKit is not available in this environment; cannot compute scores.')
            return []

        scores: List[MoleculeScore] = []
        suppl = Chem.SDMolSupplier(job.outfile, removeHs=False)

        # Prepare receptor PDBQT once per job
        receptor_pdbqt = self._prepare_receptor_pdbqt(job)

        for idx, mol in enumerate(m for m in suppl if m):
            title = mol.GetProp('_Name') if mol.HasProp('_Name') else f'mol_{idx+1}'
            # Compute VINAScore via AutoDock-GPU if possible
            vina = None
            if receptor_pdbqt:
                lig_pdbqt = self._prepare_ligand_pdbqt(mol, idx, job)
                if lig_pdbqt:
                    vina = self._run_autodock_gpu(receptor_pdbqt, lig_pdbqt, job, idx)
            try:
                qed_val = float(QED.qed(mol))
            except Exception:
                qed_val = None
            sa_val = None
            if have_sa:
                try:
                    sa_val = float(rdMolDescriptors.CalcSyntheticAccessibilityScore(mol))
                except Exception:
                    sa_val = None
            scores.append(MoleculeScore(index=idx, title=title, vina_score=vina, qed=qed_val, sa=sa_val))
        return scores

    def convert_sdf_to_pdb_index(self, job_id: str, index: int) -> Optional[str]:
        """Convert one specific molecule (0-based index) from SDF to a PDB file and return its path."""
        job = self.get_job(job_id)
        if not job or not os.path.isfile(job.outfile):
            return None
        outdir = os.path.dirname(job.outfile)
        pdb = os.path.join(outdir, f'output_{index+1}.pdb')
        if os.path.isfile(pdb):
            return pdb
        script = (
            "import sys; from rdkit import Chem; "
            "idx=int(sys.argv[2]); mols=[m for m in Chem.SDMolSupplier(sys.argv[1], removeHs=False) if m]; "
            "assert 0 <= idx < len(mols), 'Index out of range'; Chem.MolToPDBFile(mols[idx], sys.argv[3])"
        )
        try:
            ret = subprocess.run([self.python_path, '-c', script, job.outfile, str(index), pdb], cwd=self.repo_path,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
            if ret.returncode == 0 and os.path.isfile(pdb):
                return pdb
            logger.error('Index conversion failed: rc=%s err=%s', ret.returncode, ret.stderr)
            return None
        except Exception as e:
            logger.exception('Index conversion error: %s', e)
            return None
        sdf = job.outfile
        if not os.path.isfile(sdf):
            return None
        outdir = os.path.dirname(sdf)
        pdb = os.path.join(outdir, 'output.pdb')
        if os.path.isfile(pdb):
            return pdb
        # Prepare a tiny RDKit script
        script = (
            "import sys; from rdkit import Chem; "
            "mols=[m for m in Chem.SDMolSupplier(sys.argv[1], removeHs=False) if m]; "
            "assert mols, 'No molecules in SDF'; "
            "mol = mols[0]; "
            "Chem.MolToPDBFile(mol, sys.argv[2])"
        )
        try:
            ret = subprocess.run([self.python_path, '-c', script, sdf, pdb], cwd=self.repo_path,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
            if ret.returncode == 0 and os.path.isfile(pdb):
                return pdb
            logger.error('RDKit conversion failed: rc=%s, err=%s', ret.returncode, ret.stderr)
            return None
        except Exception as e:
            logger.exception('RDKit conversion error: %s', e)
            return None


# Singleton helper
_runner_instance: Optional[DiffSBDDRunner] = None


def get_runner(uploads_dir: str) -> DiffSBDDRunner:
    global _runner_instance
    if _runner_instance is None:
        _runner_instance = DiffSBDDRunner(uploads_dir=uploads_dir)
    return _runner_instance
