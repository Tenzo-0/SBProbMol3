from flask import render_template, request, jsonify, session, redirect, url_for
import os
import sys
import uuid
import shutil
import subprocess
from urllib.parse import urlparse
from app import app
from data.proteins import get_proteins, get_protein_by_id, search_proteins, add_new_protein, detect_binding_pockets
from flask import send_from_directory
# fpocket integration removed — rely on built-in detector
# NOTE: Keep this import runnable when executing `python app.py` from this folder
from utils.diffsbdd import get_runner

MD_UPLOAD_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database/MD'))
MD_TEMP_DIR = os.path.join(MD_UPLOAD_DIR, 'temp')
os.makedirs(MD_UPLOAD_DIR, exist_ok=True)
os.makedirs(MD_TEMP_DIR, exist_ok=True)

def _get_ligands_from_uploads():
    """Build a simple ligand list from files in the database folder.
    Currently treats any .pdb/.mol/.sdf file in database/ as a ligand entry.
    """
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database/ligands'))
    os.makedirs(uploads_dir, exist_ok=True)
    ligands = []
    for fname in os.listdir(uploads_dir):
        lower = fname.lower()
        if not (lower.endswith('.pdb') or lower.endswith('.mol') or lower.endswith('.sdf')):
            continue
        base, ext = os.path.splitext(fname)
        ligands.append({
            'id': base.lower().replace(' ', '_').replace('-', '_'),
            'name': base,
            'filename': fname,
            'description': f'Uploaded ligand file ({ext.lstrip(".").upper()})',
            'tags': ['Ligand', 'Uploaded']
        })
    # Sort alphabetically by name for stable ordering
    ligands.sort(key=lambda x: x['name'].lower())
    return ligands


@app.route('/')
def index():
    """Main application page"""
    proteins = get_proteins()
    ligands = _get_ligands_from_uploads()
    return render_template('index.html', proteins=proteins, ligands=ligands)


@app.route('/md')
def md_view():
    """Molecular dynamics viewer page.

    Accepts query params:
      - job_id: DiffSBDD job identifier (optional)
      - index: molecule index to render (default 0)
      - pdb_url: explicit URL to a PDB trajectory (overrides job_id)
    """
    job_id = request.args.get('job_id', '').strip()
    index = request.args.get('index', '0').strip() or '0'
    pdb_url = request.args.get('pdb_url')

    # If a job id is provided and no explicit pdb_url, use generated PDB route
    if job_id and not pdb_url:
        try:
            idx_int = int(index)
        except ValueError:
            idx_int = 0
        pdb_url = url_for('diffsbdd_result_pdb_index', job_id=job_id, index=idx_int)

    # Fallback demo path
    if not pdb_url:
        pdb_url = '/database/md/md.pdb'

    return render_template('md.html', pdb_url=pdb_url, job_id=job_id, index=index)

@app.route('/search_proteins')
def search_proteins_route():
    """Search proteins by query"""
    query = request.args.get('q', '').strip()
    if query:
        proteins = search_proteins(query)
    else:
        proteins = get_proteins()
    
    return render_template('protein_list.html', proteins=proteins)

@app.route('/protein/<protein_id>')
def get_protein_info(protein_id):
    """Get detailed protein information"""
    protein = get_protein_by_id(protein_id)
    if protein:
        return render_template('protein_info.html', protein=protein)
    else:
        return render_template('protein_info.html', protein=None)

# Removed /add_protein route; file uploads handled by /upload_pdb

@app.route('/upload_pdb', methods=['POST'])
def upload_pdb():
    """Upload a PDB file to the database folder and register a simple protein entry"""
    try:
        if 'pdb_file' not in request.files:
            return jsonify({'success': False, 'error': 'No file part'}), 400
        file = request.files['pdb_file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No selected file'}), 400

        filename = file.filename
        if not filename.lower().endswith('.pdb'):
            return jsonify({'success': False, 'error': 'Invalid file type'}), 400

        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database/proteins'))
        os.makedirs(uploads_dir, exist_ok=True)
        save_path = os.path.join(uploads_dir, filename)
        file.save(save_path)

        # Build a simple protein record for the response only
        base = os.path.splitext(filename)[0]
        protein = {
            'id': base.lower().replace(' ', '_').replace('-', '_'),
            'name': base,
            'pdb_id': base[:4].upper(),
            'description': f'Uploaded from {filename}',
            'tags': ['Uploaded'],
            'resolution': 'N/A',
            'chains': 1,
            'molecular_weight': 'Unknown',
            'organism': 'Unknown',
            'function': '',
            'filename': filename,
        }

        return jsonify({'success': True, 'protein': protein})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/upload_ligand', methods=['POST'])
def upload_ligand():
    """Upload a ligand file (PDB/MOL/SDF) into database and return a simple ligand entry."""
    try:
        if 'ligand_file' not in request.files:
            return jsonify({'success': False, 'error': 'No file part'}), 400
        file = request.files['ligand_file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No selected file'}), 400

        filename = file.filename
        lower = filename.lower()
        if not (lower.endswith('.pdb') or lower.endswith('.mol') or lower.endswith('.sdf')):
            return jsonify({'success': False, 'error': 'Invalid file type'}), 400

        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database/ligands'))
        os.makedirs(uploads_dir, exist_ok=True)
        save_path = os.path.join(uploads_dir, filename)
        file.save(save_path)

        base, ext = os.path.splitext(filename)
        ligand = {
            'id': base.lower().replace(' ', '_').replace('-', '_'),
            'name': base,
            'filename': filename,
            'description': f'Uploaded ligand file ({ext.lstrip(".").upper()})',
            'tags': ['Ligand', 'Uploaded']
        }

        return jsonify({'success': True, 'ligand': ligand})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/md/upload', methods=['POST'])
def md_upload():
    """Upload an MD ligand/trajectory file into database/MD and return its URL."""
    try:
        if 'md_file' not in request.files:
            return jsonify({'success': False, 'error': 'No file part'}), 400
        file = request.files['md_file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No selected file'}), 400

        filename = file.filename
        lower = filename.lower()
        if not (lower.endswith('.pdb') or lower.endswith('.pdbqt') or lower.endswith('.mol') or lower.endswith('.sdf')):
            return jsonify({'success': False, 'error': 'Invalid file type'}), 400

        os.makedirs(MD_UPLOAD_DIR, exist_ok=True)
        save_path = os.path.join(MD_UPLOAD_DIR, filename)
        file.save(save_path)

        url = url_for('serve_md_file', filename=filename)
        base, _ = os.path.splitext(filename)
        return jsonify({
            'success': True,
            'name': base,
            'filename': filename,
            'url': url,
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/database/<path:filename>')
def serve_upload(filename):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database/proteins'))
    return send_from_directory(uploads_dir, filename)


@app.route('/database/md/<path:filename>')
def serve_md_file(filename):
    return send_from_directory(MD_UPLOAD_DIR, filename)


@app.route('/database/md/temp/<path:filename>')
def serve_md_temp_file(filename):
    # Basic traversal guard
    if '..' in filename.split('/'):
        return jsonify({'success': False, 'error': 'Invalid path'}), 400
    return send_from_directory(MD_TEMP_DIR, filename)


@app.route('/md/list')
def md_list():
    """List available MD files in the MD upload directory."""
    try:
        files = []
        for fname in os.listdir(MD_UPLOAD_DIR):
            lower = fname.lower()
            if not (lower.endswith('.pdb') or lower.endswith('.pdbqt') or lower.endswith('.mol') or lower.endswith('.sdf')):
                continue
            url = url_for('serve_md_file', filename=fname)
            base, ext = os.path.splitext(fname)
            files.append({
                'name': base,
                'filename': fname,
                'ext': ext.lstrip('.').lower(),
                'url': url,
            })
        files.sort(key=lambda x: x['name'].lower())
        return jsonify({'success': True, 'files': files})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


def _resolve_md_path_from_url(md_url: str) -> str:
    """Convert a /database/md/... URL to an absolute file path inside MD_UPLOAD_DIR."""
    if not md_url:
        raise ValueError('Missing MD file url')
    parsed = urlparse(md_url)
    path = parsed.path or ''
    prefix = '/database/md/'
    if not path.startswith(prefix):
        raise ValueError('MD file must reside under /database/md/')
    rel = os.path.normpath(path[len(prefix):])
    if rel.startswith('..'):
        raise ValueError('Invalid MD path')
    abs_path = os.path.join(MD_UPLOAD_DIR, rel)
    if not os.path.isfile(abs_path):
        raise FileNotFoundError('Input file not found')
    return abs_path


def _split_sdf_to_files(sdf_path: str, out_dir: str, prefix: str):
    os.makedirs(out_dir, exist_ok=True)
    with open(sdf_path, 'r', encoding='utf-8', errors='ignore') as fh:
        content = fh.read()
    # Split by $$$$ delimiters, keep non-empty blocks
    chunks = [c.strip() for c in content.split('$$$$') if c.strip()]
    written = []
    for i, chunk in enumerate(chunks, 1):
        fname = f"{prefix}_{i}.sdf"
        fpath = os.path.join(out_dir, fname)
        with open(fpath, 'w', encoding='utf-8') as out:
            out.write(chunk)
            out.write('\n$$$$\n')
        written.append(fname)
    return written


@app.route('/md/run', methods=['POST'])
def md_run():
    """Run FairChem MD on a selected file, store outputs in database/MD/temp/<job_id>, convert to PDB, and return its URL."""
    try:
        payload = request.get_json(force=True, silent=True) or {}
        md_url = payload.get('url') or payload.get('file')
        steps = int(payload.get('steps', 2000))
        dt_fs = float(payload.get('dt_fs', 0.5))
        temp_k = float(payload.get('temp_k', 300.0))
        friction_fs = float(payload.get('friction_fs', 0.01))
        device = payload.get('device', 'cuda')
        task = payload.get('task', 'omol')
        model = payload.get('model', 'uma-s-1p1')

        input_path = _resolve_md_path_from_url(md_url)
        ext = os.path.splitext(input_path)[1].lower()

        # If the user already provided a PDB, just return it directly for NGL.
        if ext == '.pdb':
            return jsonify({
                'success': True,
                'job_id': None,
                'pdb_url': md_url,
                'passthrough': True,
                'cached': True,
            })

        if ext != '.sdf':
            return jsonify({'success': False, 'error': 'Only SDF is supported for MD simulation'}), 400

        base_name = os.path.splitext(os.path.basename(input_path))[0]
        cache_pdb = os.path.join(MD_TEMP_DIR, f"{base_name}.pdb")
        if os.path.isfile(cache_pdb):
            pdb_url = url_for('serve_md_temp_file', filename=os.path.basename(cache_pdb))
            return jsonify({
                'success': True,
                'job_id': base_name,
                'pdb_url': pdb_url,
                'passthrough': False,
                'cached': True,
            })

        job_id = uuid.uuid4().hex[:12]
        job_dir = os.path.join(MD_TEMP_DIR, f"{base_name}-{job_id}")
        os.makedirs(job_dir, exist_ok=True)

        # Copy input locally to keep outputs grouped per job
        local_input = os.path.join(job_dir, os.path.basename(input_path))
        shutil.copy2(input_path, local_input)

        md_script = os.path.join(os.path.dirname(__file__), 'utils', 'MD.py')
        if not os.path.isfile(md_script):
            return jsonify({'success': False, 'error': 'MD script missing'}), 500

        md_cmd = [
            sys.executable, md_script,
            '--input', local_input,
            '--traj', 'md.traj',
            '--log', 'md.log',
            '--steps', str(steps),
            '--dt_fs', str(dt_fs),
            '--temp', str(temp_k),
            '--friction_fs', str(friction_fs),
            '--device', device,
            '--task', task,
            '--model', model,
        ]

        run_result = subprocess.run(md_cmd, cwd=job_dir, capture_output=True, text=True)
        if run_result.returncode != 0:
            return jsonify({'success': False, 'error': run_result.stderr or 'MD run failed'}), 500

        convert_cmd = ['ase', 'convert', 'md.traj', cache_pdb]
        convert_result = subprocess.run(convert_cmd, cwd=job_dir, capture_output=True, text=True)
        if convert_result.returncode != 0:
            return jsonify({'success': False, 'error': convert_result.stderr or 'Failed to convert traj to pdb'}), 500

        if not os.path.isfile(cache_pdb):
            return jsonify({'success': False, 'error': 'PDB output missing'}), 500

        pdb_url = url_for('serve_md_temp_file', filename=os.path.basename(cache_pdb))
        return jsonify({
            'success': True,
            'job_id': job_id,
            'pdb_url': pdb_url,
            'passthrough': False,
            'cached': False,
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/md/import_generated', methods=['POST'])
def md_import_generated():
    """Split a generated job SDF into per-ligand SDFs and place them in MD database."""
    try:
        payload = request.get_json(force=True, silent=True) or {}
        job_id = (payload.get('job_id') or '').strip()
        if not job_id:
            return jsonify({'success': False, 'error': 'job_id required'}), 400

        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database'))
        runner = get_runner(uploads_dir)
        job = runner.get_job(job_id)
        if not job or not os.path.isfile(job.outfile):
            return jsonify({'success': False, 'error': 'Job or SDF not found'}), 404

        written = _split_sdf_to_files(job.outfile, MD_UPLOAD_DIR, job_id)
        files = [{
            'name': fname,
            'url': url_for('serve_md_file', filename=fname),
        } for fname in written]

        return jsonify({'success': True, 'count': len(written), 'files': files})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/sampling_parameters')
def sampling_parameters():
    """Show sampling parameters interface"""
    # Single-page app: redirect to main with step=3
    return redirect(url_for('index', step=3))

@app.route('/save_sampling_config', methods=['POST'])
def save_sampling_config():
    """Save sampling configuration"""
    # Clamp n_samples to a maximum of 100 to match UI limits
    n_samples = int(request.form.get('n_samples', 10))
    if n_samples > 100:
        n_samples = 100

    config = {
        'n_samples': n_samples,
        'batch_cost': int(request.form.get('batch_cost', 20)),
        'num_workers': int(request.form.get('num_workers', 12)),
        'coord_noise_scale': float(request.form.get('coord_noise_scale', 0.1)),
        'model': request.form.get('model', 'Conditional model (Binding MOAD)'),
        'timesteps': int(request.form.get('timesteps', 100)),
        'resamplings': int(request.form.get('resamplings', 1)),
        'jump_length': int(request.form.get('jump_length', 1)),
        'keep_all_fragments': request.form.get('keep_all_fragments') == 'on',
        'sanitize': request.form.get('sanitize') == 'on',
        'relax': request.form.get('relax') == 'on'
    }
    
    # Store in session
    session['sampling_config'] = config
    return jsonify({'success': True, 'config': config})

@app.route('/pocket_definition')
def pocket_definition():
    """Show pocket definition interface"""
    # Single-page app: redirect to main with step=2
    return redirect(url_for('index', step=2))


@app.route('/diffsbdd/start', methods=['POST'])
def diffsbdd_start():
    """Start a DiffSBDD generation job using the saved sampling config and current selection.
    Expects form fields: ref_ligand (e.g., 'A:300'), pdb_id, filename (optional)
    """
    try:
        # Gather inputs
        ref_ligand = request.form.get('ref_ligand', '').strip()
        pdb_id = request.form.get('pdb_id')
        filename = request.form.get('filename')
        if not ref_ligand:
            return jsonify({'success': False, 'error': 'Missing ref_ligand'}), 400

        cfg = session.get('sampling_config', {})
        n_samples = int(cfg.get('n_samples', 10))
        batch_cost = int(cfg.get('batch_cost', 20))
        num_workers = int(cfg.get('num_workers', 12))
        coord_noise_scale = float(cfg.get('coord_noise_scale', 0.1))
        timesteps = int(cfg.get('timesteps', 100))
        resamplings = int(cfg.get('resamplings', 1))
        jump_length = int(cfg.get('jump_length', 1))
        keep_all_fragments = bool(cfg.get('keep_all_fragments', False))
        sanitize = bool(cfg.get('sanitize', True))
        relax = bool(cfg.get('relax', True))

        # Initialize runner
        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database'))
        runner = get_runner(uploads_dir)
        if not runner.is_configured():
            return jsonify({
                'success': False,
                'error': 'DiffSBDD not configured on server. Set DIFFSBDD_REPO, DIFFSBDD_PYTHON, DIFFSBDD_CHECKPOINT environment variables.'
            }), 500

        job = runner.start_job(
            pdb_id=pdb_id,
            filename=filename,
            ref_ligand=ref_ligand,
            n_samples=n_samples,
            batch_cost=batch_cost,
            num_workers=num_workers,
            coord_noise_scale=coord_noise_scale,
            timesteps=timesteps,
            resamplings=resamplings,
            jump_length=jump_length,
            keep_all_fragments=keep_all_fragments,
            sanitize=sanitize,
            relax=relax,
        )
        return jsonify({'success': True, 'job_id': job.id})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/diffsbdd/status/<job_id>')
def diffsbdd_status(job_id):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database'))
    runner = get_runner(uploads_dir)
    status = runner.get_status(job_id)
    if not status.get('exists'):
        return jsonify({'success': False, 'error': 'Job not found'}), 404
    return jsonify({'success': True, 'status': status})


@app.route('/diffsbdd/result/<job_id>/sdf')
def diffsbdd_result_sdf(job_id):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database'))
    runner = get_runner(uploads_dir)
    job = runner.get_job(job_id)
    if not job or not os.path.isfile(job.outfile):
        return jsonify({'success': False, 'error': 'Result not found'}), 404
    rel = os.path.relpath(job.outfile, uploads_dir)
    d, f = os.path.split(rel)
    return send_from_directory(os.path.join(uploads_dir, d), f, mimetype='chemical/x-mdl-sdfile')


@app.route('/diffsbdd/result/<job_id>/pdb')
def diffsbdd_result_pdb(job_id):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database'))
    runner = get_runner(uploads_dir)
    pdb_path = runner.convert_sdf_to_pdb(job_id)
    if not pdb_path or not os.path.isfile(pdb_path):
        return jsonify({'success': False, 'error': 'PDB not available'}), 404
    rel = os.path.relpath(pdb_path, uploads_dir)
    d, f = os.path.split(rel)
    return send_from_directory(os.path.join(uploads_dir, d), f, mimetype='chemical/x-pdb')


@app.route('/view_result/<job_id>')
def view_result(job_id):
    """Minimal viewer page that loads the converted PDB into NGL."""
    # Inline simple HTML to avoid new template files
    pdb_url = url_for('diffsbdd_result_pdb', job_id=job_id)
    html = """
    <!doctype html>
    <html>
    <head>
        <meta charset='utf-8'>
        <title>Result Viewer</title>
        <script src='https://unpkg.com/ngl@latest/dist/ngl.js'></script>
        <style>html,body,#stage{width:100%;height:100%;margin:0;background:#000;}</style>
    </head>
    <body>
        <div id='stage'></div>
        <script>
            var stage = new NGL.Stage('stage', { backgroundColor: 'black' });
            stage.loadFile('{pdb_url}', { ext: 'pdb' }).then(function(comp){
                comp.addRepresentation('cartoon');
                comp.addRepresentation('ball+stick', { sele: 'ligand' });
                comp.autoView();
            }).catch(function(e){ console.error(e); alert('Failed to load result'); });
        </script>
    </body>
    </html>
    """
    return html.replace('{pdb_url}', pdb_url)


@app.route('/diffsbdd/result/<job_id>/list')
def diffsbdd_result_list(job_id):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database'))
    runner = get_runner(uploads_dir)
    lst = runner.list_sdf_molecules(job_id)
    return jsonify({'success': True, 'molecules': lst})


@app.route('/diffsbdd/result/<job_id>/pdb/<int:index>')
def diffsbdd_result_pdb_index(job_id, index):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'database'))
    runner = get_runner(uploads_dir)
    pdb_path = runner.convert_sdf_to_pdb_index(job_id, index)
    if not pdb_path or not os.path.isfile(pdb_path):
        return jsonify({'success': False, 'error': 'PDB index not available'}), 404
    rel = os.path.relpath(pdb_path, uploads_dir)
    d, f = os.path.split(rel)
    return send_from_directory(os.path.join(uploads_dir, d), f, mimetype='chemical/x-pdb')


@app.route('/view_result/<job_id>/<int:index>')
def view_result_index(job_id, index):
    pdb_url = url_for('diffsbdd_result_pdb_index', job_id=job_id, index=index)
    html = """
    <!doctype html>
    <html>
    <head>
        <meta charset='utf-8'>
        <title>Result Viewer</title>
        <script src='https://unpkg.com/ngl@latest/dist/ngl.js'></script>
        <style>html,body,#stage{width:100%;height:100%;margin:0;background:#000;}</style>
    </head>
    <body>
        <div id='stage'></div>
        <script>
            var stage = new NGL.Stage('stage', { backgroundColor: 'black' });
            stage.loadFile('{pdb_url}', { ext: 'pdb' }).then(function(comp){
                comp.addRepresentation('cartoon');
                comp.addRepresentation('ball+stick', { sele: 'ligand' });
                comp.autoView();
            }).catch(function(e){ console.error(e); alert('Failed to load result'); });
        </script>
    </body>
    </html>
    """
    return html.replace('{pdb_url}', pdb_url)
