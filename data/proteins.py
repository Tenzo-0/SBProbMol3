"""
Protein data utilities.

This version derives proteins dynamically from files in the uploads/ folder,
so there's no in-memory PROTEINS_DATA to maintain.
"""

import os


def _uploads_dir():
    # data/proteins.py -> project_root/uploads
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    return os.path.normpath(os.path.join(root, 'uploads'))


def _protein_from_filename(filename: str):
    base = os.path.splitext(filename)[0]
    return {
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


def get_proteins():
    """Return all available protein structures (derived from uploads/*.pdb)."""
    uploads = _uploads_dir()
    try:
        files = [f for f in os.listdir(uploads) if f.lower().endswith('.pdb')]
    except FileNotFoundError:
        files = []
    return [_protein_from_filename(f) for f in sorted(files)]


def get_protein_by_id(protein_id):
    """Get a specific protein by its ID (derived)."""
    for protein in get_proteins():
        if protein['id'] == protein_id:
            return protein
    return None


def search_proteins(query):
    """Search proteins by name, PDB ID, or description."""
    proteins = get_proteins()
    if not query:
        return proteins

    q = str(query).lower().strip()
    results = []
    for p in proteins:
        searchable_text = f"{p['name']} {p['pdb_id']} {p['description']} {' '.join(p.get('tags', []))}".lower()
        if q in searchable_text:
            results.append(p)
    return results
def get_protein_details(protein_id):
    """Get detailed information about a protein including structural data"""
    protein = get_protein_by_id(protein_id)
    if not protein:
        return None
    
    # Add additional computed details
    details = protein.copy()
    details.update({
        "binding_sites": get_binding_sites(protein_id),
        "known_ligands": get_known_ligands(protein_id),
        "druggability_score": calculate_druggability_score(protein_id),
        "similar_structures": find_similar_structures(protein_id)
    })
    
    return details

def get_binding_sites(protein_id):
    """Get known binding sites for a protein"""
    binding_sites = {
        "hiv1_protease": [
            {"name": "Active Site", "residues": ["ASP25", "ASP125"], "volume": "650 Ų"},
            {"name": "Flap Region", "residues": ["ILE50", "ILE150"], "volume": "200 Ų"}
        ],
        "sars_cov2_mpro": [
            {"name": "Active Site", "residues": ["CYS145", "HIS41"], "volume": "800 Ų"},
            {"name": "Allosteric Site", "residues": ["ASP289", "GLU290"], "volume": "300 Ų"}
        ],
        "human_cdk2": [
            {"name": "ATP Binding Site", "residues": ["LYS33", "ASP145"], "volume": "550 Ų"},
            {"name": "Cyclin Binding", "residues": ["HIS84", "LEU85"], "volume": "400 Ų"}
        ]
    }
    return binding_sites.get(protein_id, [])

def get_known_ligands(protein_id):
    """Get known ligands that bind to this protein"""
    ligands = {
        "hiv1_protease": ["Saquinavir", "Ritonavir", "Indinavir", "Nelfinavir"],
        "sars_cov2_mpro": ["Nirmatrelvir", "Ensitrelvir", "Simeprevir"],
        "human_cdk2": ["Roscovitine", "Purvalanol", "ATP"],
        "egfr_kinase": ["Gefitinib", "Erlotinib", "Lapatinib"],
        "bcl2_protein": ["Venetoclax", "Navitoclax", "Gossypol"]
    }
    return ligands.get(protein_id, [])

def calculate_druggability_score(protein_id):
    """Calculate a mock druggability score"""
    scores = {
        "hiv1_protease": 0.85,
        "sars_cov2_mpro": 0.78,
        "human_cdk2": 0.72,
        "egfr_kinase": 0.88,
        "bcl2_protein": 0.65,
        "p53_tumor": 0.42,
        "ace2_receptor": 0.58,
        "insulin_receptor": 0.51
    }
    return scores.get(protein_id, 0.50)

def find_similar_structures(protein_id):
    """Find structurally similar proteins"""
    similarities = {
        "hiv1_protease": ["1HTM", "1HVH", "1HOS"],
        "sars_cov2_mpro": ["2HOB", "3V3M", "1Q2W"],
        "human_cdk2": ["1FIN", "1GII", "1H1R"],
        "egfr_kinase": ["1M17", "2ITY", "2J6M"]
    }
    return similarities.get(protein_id, [])

def add_new_protein(protein_data):
    """No-op placeholder retained for compatibility.

    Proteins are discovered from files in uploads/, so there's nothing to add here.
    Returns the input unchanged.
    """
    # Ensure an id field exists for downstream code that may rely on it
    if protein_data and 'id' not in protein_data and 'name' in protein_data:
        protein_data['id'] = protein_data['name'].lower().replace(' ', '_').replace('-', '_')
    return protein_data

def detect_binding_pockets(pdb_id, sequence=None):
    """Auto-detect binding pockets from PDB structure or sequence"""
    # Mock implementation - in reality this would use computational tools
    mock_pockets = [
        {
            "id": "pocket_1",
            "name": "Primary Binding Site",
            "center": {"x": 12.5, "y": 8.3, "z": -4.2},
            "radius": 8.5,
            "volume": "720 Ų",
            "druggability_score": 0.85,
            "residues": ["ASP25", "ILE50", "GLY48", "VAL82"],
            "confidence": 0.92
        },
        {
            "id": "pocket_2", 
            "name": "Allosteric Site",
            "center": {"x": -5.1, "y": 15.7, "z": 2.8},
            "radius": 6.2,
            "volume": "480 Ų",
            "druggability_score": 0.67,
            "residues": ["GLU166", "ARG188", "TYR181"],
            "confidence": 0.78
        },
        {
            "id": "pocket_3",
            "name": "Secondary Site", 
            "center": {"x": 8.9, "y": -3.4, "z": 11.6},
            "radius": 5.8,
            "volume": "350 Ų",
            "druggability_score": 0.54,
            "residues": ["LYS103", "ASN142", "SER144"],
            "confidence": 0.65
        }
    ]
    return mock_pockets

def validate_protein_sequence(sequence):
    """Validate protein sequence format"""
    if not sequence:
        return False
    
    # Check if sequence contains only valid amino acid codes
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    sequence_clean = sequence.upper().replace(' ', '').replace('\n', '')
    
    return all(aa in valid_aa for aa in sequence_clean)

def predict_protein_properties(sequence):
    """Predict basic protein properties from sequence"""
    if not validate_protein_sequence(sequence):
        return None
    
    sequence_clean = sequence.upper().replace(' ', '').replace('\n', '')
    
    # Mock property predictions
    return {
        "length": len(sequence_clean),
        "molecular_weight": f"{len(sequence_clean) * 0.11:.1f} kDa",  # Rough estimate
        "isoelectric_point": 7.2,
        "hydrophobicity": 0.45,
        "predicted_structure": "Mixed α/β",
        "signal_peptide": False,
        "transmembrane_regions": 0
    }
