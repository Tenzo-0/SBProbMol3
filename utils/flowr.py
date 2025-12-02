import json
import logging
import os
import subprocess
import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Dict, Optional

from urllib.request import urlretrieve

try:
    from flowr_config import FLOWR_SETTINGS  # type: ignore
except Exception:
    FLOWR_SETTINGS = {}


logger = logging.getLogger(__name__)


@dataclass
class FlowrJob:
    id: str
    state: str = "queued"
    message: str = ""
    started_at: float = field(default_factory=time.time)
    finished_at: Optional[float] = None
    target_name: str = ""
    outfile: str = ""
    log_file: str = ""
    requested_molecules: int = 0
    proc: Optional[subprocess.Popen] = None


class FlowrRunner:
    """Wrapper around FLOWR's generate_from_pdb script."""

    def __init__(self,
                 uploads_dir: str,
                 repo_path: Optional[str] = None,
                 python_path: Optional[str] = None,
                 checkpoint_path: Optional[str] = None):
        self.uploads_dir = os.path.abspath(uploads_dir)
        cfg = FLOWR_SETTINGS or {}
        self.repo_path = repo_path or cfg.get('repo_path') or os.environ.get('FLOWR_REPO', '')
        self.python_path = python_path or cfg.get('python_path') or os.environ.get('FLOWR_PYTHON', '')
        self.checkpoint_path = checkpoint_path or cfg.get('checkpoint_path') or os.environ.get('FLOWR_CHECKPOINT', '')
        self.data_path = cfg.get('data_path') or os.environ.get('FLOWR_DATA_PATH', '')
        self.dataset = cfg.get('dataset') or os.environ.get('FLOWR_DATASET', '')
        self.arch = cfg.get('arch') or os.environ.get('FLOWR_ARCH', 'pocket')
        self.pocket_type = cfg.get('pocket_type') or os.environ.get('FLOWR_POCKET_TYPE', 'holo')
        self.batch_cost = int(cfg.get('batch_cost') or os.environ.get('FLOWR_BATCH_COST', '100') or 100)
        self.num_workers = int(cfg.get('num_workers') or os.environ.get('FLOWR_NUM_WORKERS', '8') or 8)
        self.gpus = int(cfg.get('gpus') or os.environ.get('FLOWR_GPUS', '1') or 1)
        self.jobs: Dict[str, FlowrJob] = {}

    @property
    def proteins_dir(self) -> str:
        return os.path.join(self.uploads_dir, 'proteins')

    @property
    def ligands_dir(self) -> str:
        return os.path.join(self.uploads_dir, 'ligands')

    def is_configured(self) -> bool:
        return all([
            self.repo_path and os.path.isdir(self.repo_path),
            self.python_path and os.path.isfile(self.python_path),
            self.checkpoint_path and os.path.isfile(self.checkpoint_path),
        ])

    def ensure_pdb_local(self, pdb_id: Optional[str], filename: Optional[str]) -> str:
        if filename:
            local = os.path.join(self.proteins_dir, filename)
            if os.path.isfile(local):
                return local
            # legacy location (uploads root)
            legacy = os.path.join(self.uploads_dir, filename)
            if os.path.isfile(legacy):
                return legacy
        if pdb_id:
            candidate = os.path.join(self.proteins_dir, f"{pdb_id}.pdb")
            if os.path.isfile(candidate):
                return candidate
            legacy = os.path.join(self.uploads_dir, f"{pdb_id}.pdb")
            if os.path.isfile(legacy):
                return legacy
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            dest = os.path.join(self.proteins_dir, f"{pdb_id}.pdb")
            os.makedirs(self.proteins_dir, exist_ok=True)
            logger.info("Downloading PDB %s -> %s", url, dest)
            urlretrieve(url, dest)
            if os.path.isfile(dest):
                return dest
        raise FileNotFoundError('Unable to locate target protein structure locally. Upload a PDB file first.')

    def ensure_ligand_local(self, filename: str) -> str:
        if not filename:
            raise FileNotFoundError('Ligand filename is missing')
        path = os.path.join(self.ligands_dir, filename)
        if os.path.isfile(path):
            return path
        legacy = os.path.join(self.uploads_dir, filename)
        if os.path.isfile(legacy):
            return legacy
        raise FileNotFoundError(f'Ligand {filename} not found in uploads/ligands')

    def start_job(self,
                  pdb_id: Optional[str],
                  filename: Optional[str],
                  ligand_filename: str,
                  sample_n_molecules: int,
                  max_sample_iter: int,
                  coord_noise_std: float,
                  pocket_cutoff: float,
                  cut_pocket: bool,
                  sample_mol_sizes: bool,
                  filter_valid_unique: bool,
                  compute_interactions: bool,
                  compute_interaction_recovery: bool,
                  protonate_generated_ligands: bool) -> FlowrJob:
        if not self.is_configured():
            raise RuntimeError('FLOWR is not configured. Set FLOWR_REPO, FLOWR_PYTHON, FLOWR_CHECKPOINT env vars.')

        pdb_path = self.ensure_pdb_local(pdb_id, filename)
        ligand_path = self.ensure_ligand_local(ligand_filename)
        target_name = os.path.splitext(os.path.basename(pdb_path))[0]

        job_id = str(uuid.uuid4())
        job_dir = os.path.join(self.uploads_dir, 'flowr', job_id)
        os.makedirs(job_dir, exist_ok=True)
        outfile = os.path.join(job_dir, f'samples_{target_name}.sdf')
        log_file = os.path.join(job_dir, 'flowr.log')

        cmd = [
            self.python_path,
            '-m', 'flowr.gen.generate_from_pdb',
            '--pdb_file', pdb_path,
            '--ligand_file', ligand_path,
            '--ckpt_path', self.checkpoint_path,
            '--save_dir', job_dir,
            '--arch', self.arch,
            '--pocket_type', self.pocket_type,
            '--gpus', str(self.gpus),
            '--num_workers', str(self.num_workers),
            '--batch_cost', str(self.batch_cost),
            '--sample_n_molecules_per_target', str(int(sample_n_molecules)),
            '--max_sample_iter', str(int(max_sample_iter)),
            '--coord_noise_std', str(coord_noise_std),
            '--pocket_cutoff', str(pocket_cutoff),
        ]
        if self.data_path:
            cmd += ['--data_path', self.data_path]
        if self.dataset:
            cmd += ['--dataset', self.dataset]
        if cut_pocket:
            cmd.append('--cut_pocket')
        if sample_mol_sizes:
            cmd.append('--sample_mol_sizes')
        if filter_valid_unique:
            cmd.append('--filter_valid_unique')
        if compute_interactions:
            cmd.append('--compute_interactions')
        if compute_interaction_recovery:
            cmd.append('--compute_interaction_recovery')
        if protonate_generated_ligands:
            cmd.append('--protonate_generated_ligands')

        env = os.environ.copy()
        existing_py_path = env.get('PYTHONPATH', '')
        env['PYTHONPATH'] = self.repo_path + (os.pathsep + existing_py_path if existing_py_path else '')

        job = FlowrJob(
            id=job_id,
            state='running',
            message='Starting FLOWR...',
            target_name=target_name,
            outfile=outfile,
            log_file=log_file,
            requested_molecules=sample_n_molecules,
        )
        self.jobs[job_id] = job

        def _run():
            try:
                with open(log_file, 'w', encoding='utf-8', errors='ignore') as lf:
                    pretty_cmd = ' '.join([f'"{c}"' if ' ' in c else c for c in cmd])
                    lf.write('Command: ' + pretty_cmd + '\n')
                    lf.flush()
                    proc = subprocess.Popen(cmd, cwd=self.repo_path, env=env,
                                            stdout=lf, stderr=subprocess.STDOUT, shell=False)
                    job.proc = proc
                    ret = proc.wait()
                job.finished_at = time.time()
                if ret == 0 and os.path.isfile(outfile):
                    job.state = 'success'
                    job.message = 'FLOWR generation completed'
                else:
                    job.state = 'failed'
                    job.message = f'FLOWR exited with code {ret}'
            except Exception as exc:
                job.finished_at = time.time()
                job.state = 'failed'
                job.message = f'Error: {exc}'
                logger.exception('FLOWR job failed: %s', exc)

        threading.Thread(target=_run, daemon=True).start()
        return job

    def _count_sdf_molecules(self, path: str) -> int:
        if not os.path.isfile(path):
            return 0
        try:
            count = 0
            with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
                for line in handle:
                    if line.strip() == '$$$$':
                        count += 1
            return count
        except Exception:
            return 0

    def get_status(self, job_id: str) -> Dict:
        job = self.jobs.get(job_id)
        if not job:
            return {'exists': False}
        molecules = self._count_sdf_molecules(job.outfile)
        progress = 5
        if job.requested_molecules:
            progress = min(95, int(5 + 90 * (molecules / max(1, job.requested_molecules))))
        if job.state == 'success':
            progress = 100
        elif job.state == 'failed':
            progress = min(progress, 99)
        tail = ''
        try:
            if os.path.isfile(job.log_file):
                with open(job.log_file, 'r', encoding='utf-8', errors='ignore') as lf:
                    lines = lf.readlines()
                    tail = ''.join(lines[-20:])
        except Exception:
            pass
        rel_outfile = os.path.relpath(job.outfile, self.uploads_dir).replace('\\', '/')
        return {
            'exists': True,
            'job_id': job.id,
            'state': job.state,
            'message': job.message,
            'progress': progress,
            'molecules': molecules,
            'outfile_url': f"/uploads/{rel_outfile}" if os.path.isfile(job.outfile) else None,
            'log_tail': tail,
        }

    def get_job(self, job_id: str) -> Optional[FlowrJob]:
        return self.jobs.get(job_id)

    def get_job_dir(self, job_id: str) -> str:
        return os.path.join(self.uploads_dir, 'flowr', job_id)

    def list_sdf_molecules(self, job_id: str):
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
            result = subprocess.run([self.python_path, '-c', script, job.outfile], cwd=self.repo_path,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
            if result.returncode == 0 and result.stdout:
                return json.loads(result.stdout.strip())
            logger.error('FLOWR list_sdf failed rc=%s err=%s', result.returncode, result.stderr)
            return []
        except Exception as exc:
            logger.exception('FLOWR list_sdf error: %s', exc)
            return []

    def convert_sdf_to_pdb_index(self, job_id: str, index: int) -> Optional[str]:
        job = self.get_job(job_id)
        if not job or not os.path.isfile(job.outfile):
            return None
        outdir = os.path.dirname(job.outfile)
        pdb_path = os.path.join(outdir, f'sample_{index+1}.pdb')
        if os.path.isfile(pdb_path):
            return pdb_path
        script = (
            "import sys; from rdkit import Chem; "
            "idx=int(sys.argv[2]); mols=[m for m in Chem.SDMolSupplier(sys.argv[1], removeHs=False) if m]; "
            "assert 0 <= idx < len(mols), 'Index out of range'; Chem.MolToPDBFile(mols[idx], sys.argv[3])"
        )
        try:
            result = subprocess.run([self.python_path, '-c', script, job.outfile, str(index), pdb_path],
                                    cwd=self.repo_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True, check=False)
            if result.returncode == 0 and os.path.isfile(pdb_path):
                return pdb_path
            logger.error('FLOWR convert index failed rc=%s err=%s', result.returncode, result.stderr)
            return None
        except Exception as exc:
            logger.exception('FLOWR convert index error: %s', exc)
            return None


_runner_instance: Optional[FlowrRunner] = None


def get_runner(uploads_dir: str) -> FlowrRunner:
    global _runner_instance
    if _runner_instance is None:
        _runner_instance = FlowrRunner(uploads_dir=uploads_dir)
    return _runner_instance
