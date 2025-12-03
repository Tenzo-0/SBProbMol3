from flask import render_template, request, jsonify, session, redirect, url_for
import os
import sys
import subprocess
from app import app
from data.proteins import get_proteins, get_protein_by_id, search_proteins, add_new_protein, detect_binding_pockets
from flask import send_from_directory
# fpocket integration removed — rely on built-in detector
from utils.diffsbdd import get_runner

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
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads/proteins'))
    return send_from_directory(uploads_dir, filename)

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
        'legend_modes': int(request.form.get('legend_modes', 20)),
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
        legend_modes = int(cfg.get('legend_modes', 20))
        timesteps = int(cfg.get('timesteps', 100))
        resamplings = int(cfg.get('resamplings', 1))
        jump_length = int(cfg.get('jump_length', 1))
        keep_all_fragments = bool(cfg.get('keep_all_fragments', False))
        sanitize = bool(cfg.get('sanitize', True))
        relax = bool(cfg.get('relax', True))

        # Initialize runner
        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
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
            num_nodes_lig=legend_modes,
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
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
    runner = get_runner(uploads_dir)
    status = runner.get_status(job_id)
    if not status.get('exists'):
        return jsonify({'success': False, 'error': 'Job not found'}), 404
    return jsonify({'success': True, 'status': status})


@app.route('/diffsbdd/result/<job_id>/sdf')
def diffsbdd_result_sdf(job_id):
        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
        runner = get_runner(uploads_dir)
        job = runner.get_job(job_id)
        if not job or not os.path.isfile(job.outfile):
                return jsonify({'success': False, 'error': 'Result not found'}), 404
        rel = os.path.relpath(job.outfile, uploads_dir)
        d, f = os.path.split(rel)
        return send_from_directory(os.path.join(uploads_dir, d), f, mimetype='chemical/x-mdl-sdfile')


@app.route('/diffsbdd/result/<job_id>/pdb')
def diffsbdd_result_pdb(job_id):
        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
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


@app.route('/diffsbdd/result/<job_id>/list')
def diffsbdd_result_list(job_id):
        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
        runner = get_runner(uploads_dir)
        lst = runner.list_sdf_molecules(job_id)
        return jsonify({'success': True, 'molecules': lst})


@app.route('/diffsbdd/result/<job_id>/scores')
def diffsbdd_result_scores(job_id):
    """Return VINAScore, QED and SA for each generated ligand.

    The VINAScore field is currently a heuristic placeholder based on heavy
    atom count unless a proper docking workflow is wired in on the
    DiffSBDD side. This endpoint is meant purely for ranking and display
    in the web UI.
    """
    uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
    runner = get_runner(uploads_dir)
    try:
        scores = runner.score_sdf_molecules(job_id)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

    out = []
    for s in scores:
        out.append({
            'index': s.index,
            'title': s.title,
            'vina_score': s.vina_score,
            'qed': s.qed,
            'sa': s.sa,
        })
    return jsonify({'success': True, 'scores': out})


@app.route('/diffsbdd/result/<job_id>/pdb/<int:index>')
def diffsbdd_result_pdb_index(job_id, index):
        uploads_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'uploads'))
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
