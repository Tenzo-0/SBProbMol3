from flask import render_template, request, jsonify, session, redirect, url_for
import os
import sys
import subprocess
from app import app
from data.proteins import get_proteins, get_protein_by_id, search_proteins, add_new_protein, detect_binding_pockets
from flask import send_from_directory
# fpocket integration removed — rely on built-in detector
from utils.flowr import get_runner

def _get_ligands_from_uploads():
    """Build a simple ligand list from files in the uploads folder.
    Currently treats any .pdb/.mol/.sdf file in uploads/ as a ligand entry.
    """
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads/ligands'))
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
    """Upload a PDB file to the uploads folder and register a simple protein entry"""
    try:
        if 'pdb_file' not in request.files:
            return jsonify({'success': False, 'error': 'No file part'}), 400
        file = request.files['pdb_file']
        if file.filename == '':
            return jsonify({'success': False, 'error': 'No selected file'}), 400

        filename = file.filename
        if not filename.lower().endswith('.pdb'):
            return jsonify({'success': False, 'error': 'Invalid file type'}), 400

        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads/proteins'))
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
    """Upload a ligand file (PDB/MOL/SDF) into uploads and return a simple ligand entry."""
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

        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads/ligands'))
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

@app.route('/uploads/<path:filename>')
def serve_upload(filename):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
    return send_from_directory(uploads_dir, filename)

@app.route('/sampling_parameters')
def sampling_parameters():
    """Show sampling parameters interface"""
    # Single-page app: redirect to main with step=3
    return redirect(url_for('index', step=3))

@app.route('/save_sampling_config', methods=['POST'])
def save_sampling_config():
    """Save sampling configuration"""
    def _int(field, default, max_value=None, min_value=0):
        try:
            value = int(request.form.get(field, default))
        except (TypeError, ValueError):
            value = default
        if max_value is not None:
            value = min(value, max_value)
        if min_value is not None:
            value = max(value, min_value)
        return value

    def _float(field, default, min_value=None):
        try:
            value = float(request.form.get(field, default))
        except (TypeError, ValueError):
            value = default
        if min_value is not None:
            value = max(value, min_value)
        return value

    config = {
        'sample_n_molecules': _int('sample_n_molecules', 32, max_value=500, min_value=1),
        'max_sample_iter': _int('max_sample_iter', 20, max_value=200, min_value=1),
        'coord_noise_std': _float('coord_noise_std', 0.0, min_value=0.0),
        'pocket_cutoff': _float('pocket_cutoff', 6.0, min_value=1.0),
        'cut_pocket': request.form.get('cut_pocket') == 'on',
        'sample_mol_sizes': request.form.get('sample_mol_sizes') == 'on',
        'filter_valid_unique': request.form.get('filter_valid_unique') == 'on',
        'compute_interactions': request.form.get('compute_interactions') == 'on',
        'compute_interaction_recovery': request.form.get('compute_interaction_recovery') == 'on',
        'protonate_generated_ligands': request.form.get('protonate_generated_ligands') == 'on'
    }
    
    # Store in session
    session['sampling_config'] = config
    return jsonify({'success': True, 'config': config})

@app.route('/pocket_definition')
def pocket_definition():
    """Show pocket definition interface"""
    # Single-page app: redirect to main with step=2
    return redirect(url_for('index', step=2))


@app.route('/flowr/start', methods=['POST'])
def flowr_start():
    """Kick off a FLOWR job using the saved sampling configuration."""
    try:
        ligand_filename = request.form.get('ligand_filename', '').strip()
        if not ligand_filename:
            return jsonify({'success': False, 'error': 'Missing ligand file selection'}), 400
        pdb_id = request.form.get('pdb_id')
        filename = request.form.get('filename')

        cfg = session.get('sampling_config', {})
        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
        runner = get_runner(uploads_dir)
        if not runner.is_configured():
            return jsonify({
                'success': False,
                'error': 'FLOWR not configured. Set FLOWR_REPO, FLOWR_PYTHON, FLOWR_CHECKPOINT (and optional FLOWR_* settings).'
            }), 500

        job = runner.start_job(
            pdb_id=pdb_id,
            filename=filename,
            ligand_filename=ligand_filename,
            sample_n_molecules=int(cfg.get('sample_n_molecules', 32)),
            max_sample_iter=int(cfg.get('max_sample_iter', 20)),
            coord_noise_std=float(cfg.get('coord_noise_std', 0.0)),
            pocket_cutoff=float(cfg.get('pocket_cutoff', 6.0)),
            cut_pocket=bool(cfg.get('cut_pocket', True)),
            sample_mol_sizes=bool(cfg.get('sample_mol_sizes', True)),
            filter_valid_unique=bool(cfg.get('filter_valid_unique', True)),
            compute_interactions=bool(cfg.get('compute_interactions', False)),
            compute_interaction_recovery=bool(cfg.get('compute_interaction_recovery', False)),
            protonate_generated_ligands=bool(cfg.get('protonate_generated_ligands', False)),
        )
        return jsonify({'success': True, 'job_id': job.id})
    except Exception as exc:
        return jsonify({'success': False, 'error': str(exc)}), 500


@app.route('/flowr/status/<job_id>')
def flowr_status(job_id):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
    runner = get_runner(uploads_dir)
    status = runner.get_status(job_id)
    if not status.get('exists'):
        return jsonify({'success': False, 'error': 'Job not found'}), 404
    return jsonify({'success': True, 'status': status})


@app.route('/flowr/result/<job_id>/sdf')
def flowr_result_sdf(job_id):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
    runner = get_runner(uploads_dir)
    job = runner.get_job(job_id)
    if not job or not os.path.isfile(job.outfile):
        return jsonify({'success': False, 'error': 'Result not found'}), 404
    rel = os.path.relpath(job.outfile, uploads_dir)
    directory, filename = os.path.split(rel)
    return send_from_directory(os.path.join(uploads_dir, directory), filename, mimetype='chemical/x-mdl-sdfile')


@app.route('/flowr/result/<job_id>/list')
def flowr_result_list(job_id):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
    runner = get_runner(uploads_dir)
    molecules = runner.list_sdf_molecules(job_id)
    return jsonify({'success': True, 'molecules': molecules})


@app.route('/flowr/result/<job_id>/pdb/<int:index>')
def flowr_result_pdb(job_id, index):
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
    runner = get_runner(uploads_dir)
    pdb_path = runner.convert_sdf_to_pdb_index(job_id, index)
    if not pdb_path or not os.path.isfile(pdb_path):
        return jsonify({'success': False, 'error': 'PDB not available at that index'}), 404
    rel = os.path.relpath(pdb_path, uploads_dir)
    directory, filename = os.path.split(rel)
    return send_from_directory(os.path.join(uploads_dir, directory), filename, mimetype='chemical/x-pdb')


@app.route('/view_result/<job_id>/<int:index>')
def view_result(job_id, index):
    pdb_url = url_for('flowr_result_pdb', job_id=job_id, index=index)
    html = f"""
    <!doctype html>
    <html>
    <head>
        <meta charset='utf-8'>
        <title>Result Viewer</title>
        <script src='https://unpkg.com/ngl@latest/dist/ngl.js'></script>
        <style>html,body,#stage{{width:100%;height:100%;margin:0;background:#000;}}</style>
    </head>
    <body>
        <div id='stage'></div>
        <script>
            var stage = new NGL.Stage('stage', {{ backgroundColor: 'black' }});
            stage.loadFile('{pdb_url}', {{ ext: 'pdb' }}).then(function(comp){{
                comp.addRepresentation('cartoon');
                comp.addRepresentation('ball+stick', {{ sele: 'ligand' }});
                comp.autoView();
            }}).catch(function(e){{ console.error(e); alert('Failed to load result'); }});
        </script>
    </body>
    </html>
    """
    return html
