import os
import uuid
import time
import threading
import subprocess
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, List

from urllib.request import urlretrieve

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem  # type: ignore
    from rdkit.Chem import SDMolSupplier, SDWriter  # type: ignore
    from rdkit.Chem import QED  # type: ignore
    from rdkit.Contrib.SA_Score import sascorer  # type: ignore
    _HAVE_SASCORE = True
except Exception:  # pragma: no cover
    Chem = None
    SDMolSupplier = None
    SDWriter = None
    QED = None
    sascorer = None
    _HAVE_SASCORE = False

try:
    # Local AutoDock-GPU pipeline helper
    from . import vina  # type: ignore
    _HAVE_VINA = True
except Exception:  # pragma: no cover
    vina = None
    _HAVE_VINA = False


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


class DiffSBDDRunner:
    """
    Thin wrapper to run Flowr's flowr.gen.generate_from_pdb as a subprocess.
    Expects a Flowr-ready Python interpreter, repository path, and checkpoint path.
    """

    def __init__(self,
                 uploads_dir: str,
                 repo_path: Optional[str] = None,
                 python_path: Optional[str] = None,
                 checkpoint_path: Optional[str] = None):
        self.uploads_dir = os.path.abspath(uploads_dir)
        DEFAULT_REPO_PATH = "/home/user/flowr_root"
        DEFAULT_PYTHON_PATH = "/home/user/anaconda3/envs/flowr_root/bin/python"
        DEFAULT_CHECKPOINT_PATH = "/home/user/flowr_root_v2.ckpt"

        self.repo_path = repo_path or DEFAULT_REPO_PATH
        self.python_path = python_path or DEFAULT_PYTHON_PATH
        self.checkpoint_path = checkpoint_path or DEFAULT_CHECKPOINT_PATH
        self.jobs: Dict[str, DiffSBDDJob] = {}

    def is_configured(self) -> bool:
        return all([
            self.repo_path and os.path.isdir(self.repo_path),
            self.python_path and os.path.isfile(self.python_path),
            self.checkpoint_path and os.path.isfile(self.checkpoint_path),
        ])

    def ensure_pdb_local(self, pdb_id: Optional[str], filename: Optional[str]) -> str:
        """Return absolute path to a local PDB file.

        Proteins are stored under database/proteins. If not present, try to download
        by pdb_id into that folder.
        """
        proteins_dir = os.path.join(self.uploads_dir, 'proteins')
        os.makedirs(proteins_dir, exist_ok=True)

        # 1) If an uploaded filename is provided and exists under database/proteins, use it
        if filename:
            local = os.path.join(proteins_dir, filename)
            if os.path.isfile(local):
                return local
        # 2) If a file with the PDB ID exists in database/proteins, use it
        if pdb_id:
            candidate = os.path.join(proteins_dir, f"{pdb_id}.pdb")
            if os.path.isfile(candidate):
                return candidate
        # 3) Try to download from RCSB using pdb_id into database/proteins
        if pdb_id:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            dest = os.path.join(proteins_dir, f"{pdb_id}.pdb")
            try:
                os.makedirs(proteins_dir, exist_ok=True)
                logger.info("Downloading PDB %s -> %s", url, dest)
                urlretrieve(url, dest)
                if os.path.isfile(dest):
                    return dest
            except Exception as e:
                logger.error("Failed to download PDB %s: %s", pdb_id, e)
                raise
        raise FileNotFoundError("No local PDB file found and download failed. Please upload a PDB file.")

    # --- SA score helpers -------------------------------------------------

    def _annotate_sascore(self, job: DiffSBDDJob) -> None:
        """Annotate each SDF molecule with SA_SCORE and QED in-place.

        Uses RDKit's official SA_Score and QED implementations if available.
        No-op if RDKit/SA_Score is unavailable or outfile is missing.
        """
        if not _HAVE_SASCORE or not os.path.isfile(job.outfile):
            return
        try:
            suppl = SDMolSupplier(job.outfile, removeHs=False)
            mols = [m for m in suppl if m]
            if not mols:
                return
            for m in mols:
                try:
                    sa = sascorer.calculateScore(m)
                    m.SetProp('SA_SCORE', f"{sa:.3f}")
                    if QED is not None:
                        qed = QED.qed(m)
                        m.SetProp('QED', f"{qed:.3f}")
                except Exception:
                    continue
            writer = SDWriter(job.outfile)
            for m in mols:
                writer.write(m)
            writer.close()
        except Exception:
            logger.exception("Failed to annotate SA scores for job %s", job.id)

    # --- Vina scoring helpers -------------------------------------------------

    def _run_vina_for_job(self, job: DiffSBDDJob, pdb_path: str) -> None:
        """Run AutoDock-GPU pipeline via utils.vina for each ligand and
        write Vina-like scores back into the SDF as VINA_SCORE.

        This is a best-effort post-processing step:
          - requires utils.vina and external tools (mk_prepare_receptor.py,
            obabel, autogrid4, autodock_gpu_128wi)
          - ignores all errors so DiffSBDD results are still usable.
        """
        if not _HAVE_VINA or vina is None:
            logger.info("Vina module not available; skipping Vina scoring for job %s", job.id)
            return
        if not os.path.isfile(job.outfile):
            return

        try:
            suppl = SDMolSupplier(job.outfile, removeHs=False)
            mols: List[Chem.Mol] = [m for m in suppl if m]
            if not mols:
                return

            # Write each molecule to a temporary SDF and score it inside a per-job Vina temp folder
            tmp_dir = os.path.join(os.path.dirname(job.outfile), "vina_tmp")
            os.makedirs(tmp_dir, exist_ok=True)
            scored: List[Chem.Mol] = []
            for idx, mol in enumerate(mols):
                try:
                    tmp_sdf = os.path.join(tmp_dir, f"_vina_tmp_{idx+1}.sdf")
                    w = SDWriter(tmp_sdf)
                    w.write(mol)
                    w.close()

                    try:
                        best, _xml = vina.run_full_pipeline(
                            protein_pdb=pdb_path,
                            ligand_sdf=tmp_sdf,
                            workdir=tmp_dir,
                            out_prefix=f"{Path(job.outfile).stem}_vina",
                        )
                    except Exception as e:  # pragma: no cover - external tools
                        logger.error("Vina pipeline failed for mol %d in job %s: %s", idx, job.id, e)
                        best = None

                    if best is not None:
                        mol.SetProp("VINA_SCORE", f"{best:.3f}")
                except Exception:
                    logger.exception("Error scoring molecule %d in job %s", idx, job.id)
                finally:
                    # Clean up temp ligand file if it exists
                    try:
                        if os.path.isfile(tmp_sdf):
                            os.remove(tmp_sdf)
                    except Exception:
                        pass
                scored.append(mol)

            # Rewrite SDF with VINA_SCORE properties
            try:
                writer = SDWriter(job.outfile)
                for m in scored:
                    writer.write(m)
                writer.close()
            except Exception:
                logger.exception("Failed to write Vina-scored SDF for job %s", job.id)
        except Exception:
            logger.exception("Unexpected error during Vina scoring for job %s", job.id)

    def start_job(self,
                  pdb_id: Optional[str],
                  filename: Optional[str],
                  ref_ligand: str,
                  n_samples: int,
                  batch_cost: int = 20,
                  num_workers: int = 12,
                  coord_noise_scale: float = 0.1,
                  timesteps: Optional[int] = None,
                  resamplings: Optional[int] = None,
                  jump_length: Optional[int] = None,
                  keep_all_fragments: bool = False,
                  sanitize: bool = False,
                  relax: bool = False) -> DiffSBDDJob:
        if not self.is_configured():
            raise RuntimeError("Flowr not configured. Set repo, python, and checkpoint paths.")

        pdb_path = self.ensure_pdb_local(pdb_id, filename)
        job_id = str(uuid.uuid4())
        outdir = os.path.join(self.uploads_dir, 'diffsbdd', job_id)
        os.makedirs(outdir, exist_ok=True)
        log_file = os.path.join(outdir, 'run.log')

        # Resolve reference ligand path: if it's just a filename, assume database/ligands
        lig_path = ref_ligand
        if lig_path and not os.path.isabs(lig_path) and not os.path.dirname(lig_path):
            lig_dir = os.path.join(self.uploads_dir, 'ligands')
            cand = os.path.join(lig_dir, lig_path)
            if os.path.isfile(cand):
                lig_path = cand

        target_name = Path(pdb_path).stem
        outfile = os.path.join(outdir, f"samples_{target_name}.sdf")
        os.makedirs(outdir, exist_ok=True)

        cmd = [
            self.python_path,
            '-m', 'flowr.gen.generate_from_pdb',
            '--pdb_file', pdb_path,
            '--ligand_file', lig_path,
            '--arch', 'pocket',
            '--pocket_type', 'holo',
            '--cut_pocket',
            '--pocket_cutoff', '7',
            '--gpus', '1',
            '--num_workers', str(int(num_workers)),
            '--batch_cost', str(int(batch_cost)),
            '--ckpt_path', self.checkpoint_path,
            '--save_dir', outdir,
            '--max_sample_iter', '20',
            '--coord_noise_scale', str(coord_noise_scale),
            '--sample_n_molecules_per_target', str(int(n_samples)),
            '--categorical_strategy', 'uniform-sample',
            '--filter_valid_unique',
            '--filter_diversity',
            '--diversity_threshold', '0.7',
        ]
        if timesteps is not None:
            cmd += ['--integration_steps', str(int(timesteps))]
        # Corrector iterations are disabled: the available checkpoint expects
        # unmasked categorical interpolation and will assert otherwise.

        job = DiffSBDDJob(id=job_id, state='running', message='Starting Flowr...', n_samples=n_samples,
                           outfile=outfile, log_file=log_file)
        self.jobs[job_id] = job

        env = os.environ.copy()
        env['PYTHONPATH'] = (
            f"{self.repo_path}:{env.get('PYTHONPATH', '')}"
            if env.get('PYTHONPATH')
            else str(self.repo_path)
        )

        def _run():
            try:
                with open(log_file, 'w', encoding='utf-8', errors='ignore') as lf:
                    lf.write('Command: ' + ' '.join([f'"{c}"' if ' ' in str(c) else str(c) for c in cmd]) + '\n')
                    lf.flush()
                    proc = subprocess.Popen(cmd, cwd=self.repo_path, stdout=lf, stderr=subprocess.STDOUT, shell=False, env=env)
                    job.proc = proc
                    ret = proc.wait()
                job.finished_at = time.time()
                if ret == 0 and os.path.isfile(outfile):
                    # Post-process: add SA_SCORE and VINA_SCORE to each generated ligand if possible
                    self._annotate_sascore(job)
                    try:
                        self._run_vina_for_job(job, pdb_path=pdb_path)
                    except Exception:
                        logger.exception("Vina scoring failed for job %s", job.id)
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
        # Build a URL for download relative to /database route
        rel_outfile = os.path.relpath(job.outfile, self.uploads_dir).replace('\\', '/')
        return {
            'exists': True,
            'job_id': job.id,
            'state': job.state,
            'message': job.message,
            'progress': progress,
            'molecules': mols,
            'outfile_url': f"/database/{rel_outfile}" if os.path.isfile(job.outfile) else None,
            'log_tail': tail,
        }

    def get_job(self, job_id: str) -> Optional[DiffSBDDJob]:
        return self.jobs.get(job_id)

    def get_job_dir(self, job_id: str) -> str:
        return os.path.join(self.uploads_dir, 'diffsbdd', job_id)

    def convert_sdf_to_pdb(self, job_id: str, select: str = 'first') -> Optional[str]:
        """Convert a job's SDF to PDB using RDKit in the DiffSBDD env. Returns output path or None."""
        job = self.get_job(job_id)
        if not job:
            return None

    def list_sdf_molecules(self, job_id: str):
        """Return list of dicts {index, title, vina, sa, qed} for molecules.

        SA and QED come from SDF properties SA_SCORE and QED if present.
        Vina score (if available) is read from VINA_SCORE or Vina_score property.
        """
        job = self.get_job(job_id)
        if not job or not os.path.isfile(job.outfile):
            return []
        script = (
            "import sys, json\n"
            "from rdkit import Chem\n"
            "mols = [m for m in Chem.SDMolSupplier(sys.argv[1], removeHs=False) if m]\n"
            "out = []\n"
            "for i, m in enumerate(mols):\n"
            "    title = m.GetProp('_Name') if m.HasProp('_Name') else f'mol_{i+1}'\n"
            "    sa = m.GetProp('SA_SCORE') if m.HasProp('SA_SCORE') else None\n"
            "    qed = m.GetProp('QED') if m.HasProp('QED') else None\n"
            "    vina = None\n"
            "    pname = None\n"
            "    for key in m.GetPropNames():\n"
            "        if key.upper() in ['VINA_SCORE', 'VINA', 'VINA_SCORE_KCAL_MOL']:\n"
            "            pname = key\n"
            "            break\n"
            "    if pname is not None:\n"
            "        vina = m.GetProp(pname)\n"
            "    out.append({'index': i, 'title': title, 'sa': sa, 'qed': qed, 'vina': vina})\n"
            "print(json.dumps(out))\n"
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


# Singleton helper
_runner_instance: Optional[DiffSBDDRunner] = None


def get_runner(uploads_dir: str) -> DiffSBDDRunner:
    global _runner_instance
    if _runner_instance is None:
        _runner_instance = DiffSBDDRunner(uploads_dir=uploads_dir)
    return _runner_instance
