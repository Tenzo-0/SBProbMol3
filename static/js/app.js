/**
 *
 * MoleculeAI - Drug Design Platform
 * Main JavaScript functionality
 */

class MoleculeAI {
    constructor() {
        this.selectedProtein = null;
        this.currentWorkflowStep = 1;
        // NGL handles
        this.nglComp = null;
        this.currentLigandComp = null;
        this.reprs = { cartoon: null, ligand: null, surface: null, proteinBallStick: null, vertices: null };
        // Pocket state
        this.detectedPockets = [];
        this.selectedPocket = null;
        this.init();
    }

    init() {
        this.setupEventListeners();
        this.initializeSliders();
        this.setupProteinSearch();
        this.setupLigandSearch();
        this.initNGLStage();
        this.setWorkflowStep(1);
        try {
            const params = new URLSearchParams(window.location.search);
            const stepParam = parseInt(params.get('step'));
            // Load previous selection if present
            const savedPid = sessionStorage.getItem('selectedProteinId');
            const savedPdb = sessionStorage.getItem('selectedPdbId');
            const savedFile = sessionStorage.getItem('selectedFilename');
            if (savedPid && (savedPdb || savedFile)) {
                this.loadProteinViewerNGL({ 
                    proteinId: savedPid, 
                    pdbId: savedPdb || '', 
                    filename: savedFile || '' 
                });
                // Enable continue step1
                const c1 = document.getElementById('continue-step1');
                if (c1) c1.disabled = false;
                this.selectedProtein = savedPid;
            }
            if ([1,2,3,4].includes(stepParam)) {
                if (stepParam > 1 && !this.selectedProtein) {
                    this.setWorkflowStep(1);
                } else {
                    this.setWorkflowStep(stepParam);
                }
            }
        } catch (_) {}
    }

    setupEventListeners() {
        // Protein selection
        document.querySelectorAll('.protein-card').forEach(card => {
            card.addEventListener('click', (e) => {
                this.selectProtein(e.currentTarget);
            });
        });

        // Workflow step navigation
        document.querySelectorAll('.workflow-step').forEach(step => {
            step.addEventListener('click', (e) => {
                const stepNumber = parseInt(e.currentTarget.dataset.step);
                if (stepNumber <= this.currentWorkflowStep || this.selectedProtein) {
                    this.setWorkflowStep(stepNumber);
                }
            });
        });

        // Step 1 continue
        const continue1 = document.getElementById('continue-step1');
        if (continue1) {
            continue1.addEventListener('click', () => {
                if (this.selectedProtein) this.setWorkflowStep(2);
            });
        }

        // Display toggles
        document.querySelectorAll('.form-check-input').forEach(toggle => {
            toggle.addEventListener('change', (e) => {
                this.handleDisplayToggle(e.target);
                // Reflect immediately in NGL
                if (e.target.id === 'show-vertices') {
                    this.ensureVerticesRepresentation(e.target.checked);
                } else {
                    this.syncUIToReprs();
                }
            });
        });

        // Representation select -> map to NGL
        const reprSelect = document.getElementById('representation-select');
        if (reprSelect) {
            reprSelect.addEventListener('change', (e) => {
                const value = e.target.value;
                this.applyRepresentation(value);
                this.syncUIToReprs();
            });
        }

        // Viewer controls
        document.querySelectorAll('.viewer-controls .btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                this.handleViewerControl(e.currentTarget);
            });
        });

        // Export buttons
        document.querySelectorAll('[class*="btn"][class*="export"], [class*="btn"][class*="download"], [class*="btn"][class*="share"]').forEach(btn => {
            btn.addEventListener('click', (e) => {
                this.handleExportAction(e.currentTarget);
            });
        });

        // Add New protein -> open file picker
        const addBtn = document.getElementById('add-new-btn');
        const fileInput = document.getElementById('pdb-file-input');
        if (addBtn && fileInput) {
            addBtn.addEventListener('click', () => fileInput.click());
            fileInput.addEventListener('change', async () => {
                if (!fileInput.files || !fileInput.files[0]) return;
                const file = fileInput.files[0];
                if (!file.name.toLowerCase().endsWith('.pdb')) {
                    this.showNotification('Please select a .pdb file', 'info');
                    return;
                }
                try {
                    const fd = new FormData();
                    fd.append('pdb_file', file);
                    const res = await fetch('/upload_pdb', { method: 'POST', body: fd });
                    const data = await res.json();
                    if (data.success) {
                        this.showNotification('PDB uploaded successfully', 'success');
                        // Simple approach: reload to refresh protein list
                        setTimeout(() => window.location.reload(), 1000);
                    } else {
                        this.showNotification(data.error || 'Upload failed', 'info');
                    }
                } catch (err) {
                    console.error(err);
                    this.showNotification('Upload error', 'info');
                } finally {
                    fileInput.value = '';
                }
            });
        }

        // Add New ligand -> open file picker
        const addLigandBtn = document.getElementById('add-new-ligand-btn');
        const ligandInput = document.getElementById('ligand-file-input');
        if (addLigandBtn && ligandInput) {
            addLigandBtn.addEventListener('click', () => ligandInput.click());
            ligandInput.addEventListener('change', async () => {
                if (!ligandInput.files || !ligandInput.files[0]) return;
                const file = ligandInput.files[0];
                const name = file.name.toLowerCase();
                if (!(name.endsWith('.pdb') || name.endsWith('.mol') || name.endsWith('.sdf'))) {
                    this.showNotification('Please select a .pdb, .mol or .sdf file', 'info');
                    return;
                }
                try {
                    const fd = new FormData();
                    fd.append('ligand_file', file);
                    const res = await fetch('/upload_ligand', { method: 'POST', body: fd });
                    const data = await res.json();
                    if (data.success) {
                        this.showNotification('Ligand uploaded successfully', 'success');
                        // For now simply reload to refresh ligand list
                        setTimeout(() => window.location.reload(), 1000);
                    } else {
                        this.showNotification(data.error || 'Ligand upload failed', 'info');
                    }
                } catch (err) {
                    console.error(err);
                    this.showNotification('Ligand upload error', 'info');
                } finally {
                    ligandInput.value = '';
                }
            });
        }

        // Pocket method cards: only manual method remains
        document.querySelectorAll('.method-card').forEach(card => {
            card.addEventListener('click', (e) => {
                document.querySelectorAll('.method-card').forEach(c => c.classList.remove('active'));
                const el = e.currentTarget;
                el.classList.add('active');
                this.currentPocketMethod = 'manual';
                document.getElementById('manual-input').style.display = 'block';
            });
        });

        // Manual reference input: create a simple pocket entry when user types a reference
        const manualRef = document.getElementById('manual-ref');
        if (manualRef) {
            manualRef.addEventListener('input', (e) => {
                const val = e.target.value.trim();
                const contBtn = document.getElementById('continue-step2');
                if (!val) {
                    this.detectedPockets = [];
                    this.selectedPocket = null;
                    this.renderPocketsList();
                    if (contBtn) contBtn.disabled = true;
                    return;
                }
                // Create a simple placeholder pocket from the reference string
                const pocket = {
                    id: `ref_${Date.now()}`,
                    name: `Ref: ${val}`,
                    center: { x: 0, y: 0, z: 0 },
                    radius: 6.0,
                    volume: `${Math.round((4/3)*Math.PI*6*6*6)}`,
                    druggability_score: 0.5,
                    residues: [val],
                    confidence: 0.6
                };
                this.detectedPockets = [pocket];
                this.selectedPocket = pocket;
                this.renderPocketsList();
                this.renderPocketProperties();
                if (contBtn) contBtn.disabled = false;
            });
        }

        // Continue step 2
        const continue2 = document.getElementById('continue-step2');
        if (continue2) {
            continue2.addEventListener('click', () => {
                if (this.selectedPocket || sessionStorage.getItem('selectedLigandFilename')) this.setWorkflowStep(3);
            });
        }

        // Sampling form (step 3)
        const samplingForm = document.getElementById('sampling-form');
        if (samplingForm) {
            samplingForm.addEventListener('submit', (e) => this.submitSamplingForm(e));
            const resetBtn = document.getElementById('reset-defaults');
            if (resetBtn) resetBtn.addEventListener('click', () => this.resetSamplingDefaults());
            // Range display
            const ns = document.getElementById('n-samples');
            const nsv = document.getElementById('n-samples-value');
            if (ns && nsv) ns.addEventListener('input', () => nsv.textContent = ns.value);
            const ts = document.getElementById('timesteps');
            const tsv = document.getElementById('timesteps-value');
            if (ts && tsv) ts.addEventListener('input', () => tsv.textContent = ts.value);
        }
    }

    initializeSliders() {
        const qualitySlider = document.getElementById('quality-slider');
        const opacitySlider = document.getElementById('opacity-slider');

        qualitySlider.addEventListener('input', (e) => {
            const value = parseInt(e.target.value);
            const labels = ['Low', 'Medium', 'High'];
            document.getElementById('quality-value').textContent = labels[value];
            this.updateRenderingQuality(value);
        });

        opacitySlider.addEventListener('input', (e) => {
            const value = e.target.value;
            document.getElementById('opacity-value').textContent = `${value}%`;
            this.updateSurfaceOpacity(value);
        });
    }

    setupProteinSearch() {
        const searchInput = document.getElementById('protein-search');
        if (searchInput) {
            let searchTimeout;
            searchInput.addEventListener('input', (e) => {
                clearTimeout(searchTimeout);
                searchTimeout = setTimeout(() => {
                    this.searchProteins(e.target.value);
                }, 300);
            });
        }
    }

    searchProteins(query) {
        const proteinCards = document.querySelectorAll('.protein-card');
        const searchQuery = query.toLowerCase().trim();

        proteinCards.forEach(card => {
            const name = card.querySelector('.protein-name').textContent.toLowerCase();
            const pdbId = card.querySelector('.pdb-id').textContent.toLowerCase();
            const description = card.querySelector('.protein-description').textContent.toLowerCase();
            
            const matches = name.includes(searchQuery) || 
                          pdbId.includes(searchQuery) || 
                          description.includes(searchQuery);

            card.style.display = matches || !searchQuery ? 'block' : 'none';
        });
    }

    setupLigandSearch() {
        const searchInput = document.getElementById('ligand-search');
        if (searchInput) {
            let searchTimeout;
            searchInput.addEventListener('input', (e) => {
                clearTimeout(searchTimeout);
                const query = e.target.value.toLowerCase().trim();
                searchTimeout = setTimeout(() => {
                    const cards = document.querySelectorAll('.ligand-card');
                    cards.forEach(card => {
                        const name = (card.querySelector('.protein-name')?.textContent || '').toLowerCase();
                        const fileSpan = card.querySelector('.pdb-id');
                        const fileTxt = (fileSpan ? fileSpan.textContent : '').toLowerCase();
                        const descEl = card.querySelector('.protein-description');
                        const desc = (descEl ? descEl.textContent : '').toLowerCase();
                        const matches = !query || name.includes(query) || fileTxt.includes(query) || desc.includes(query);
                        card.style.display = matches ? 'block' : 'none';
                    });
                }, 200);
            });
        }
    }

    initNGLStage() {
        const container = document.getElementById('protein-viewer');
        if (!container || typeof NGL === 'undefined') return;
        // Clear placeholder
        container.innerHTML = '';
        // Create NGL Stage
        this.nglStage = new NGL.Stage('protein-viewer', { backgroundColor: "#000010" });
        window.addEventListener('resize', () => {
            this.nglStage.handleResize();
        });
    }

    // --- Ligand selection helpers (Step 2) ---
    setupLigandSelection() {
        const cards = document.querySelectorAll('.ligand-card');
        const contBtn = document.getElementById('continue-step2');
        cards.forEach(card => {
            card.addEventListener('click', () => {
                // Clear previous selection
                document.querySelectorAll('.ligand-card').forEach(c => c.classList.remove('selected'));
                card.classList.add('selected');
                // Store selected ligand filename for later (e.g., DiffSBDD ref ligand)
                const filename = card.dataset.filename || '';
                try {
                    sessionStorage.setItem('selectedLigandFilename', filename);
                } catch (_) {}
                if (contBtn) contBtn.disabled = false;
            });
        });
    }

    selectProtein(proteinCard) {
        // Remove previous selection
        document.querySelectorAll('.protein-card').forEach(card => {
            card.classList.remove('selected');
        });

        // Add selection to clicked card
        proteinCard.classList.add('selected');
        
        const proteinId = proteinCard.dataset.proteinId;
        this.selectedProtein = proteinId;

        // Persist selection across workflow pages
        try {
            sessionStorage.setItem('selectedProteinId', proteinId);
            if (proteinCard.dataset.pdbId) {
                sessionStorage.setItem('selectedPdbId', proteinCard.dataset.pdbId);
            }
            sessionStorage.setItem('selectedFilename', proteinCard.dataset.filename || '');
        } catch (_) {}

        // Enable continue button
    const continue1 = document.getElementById('continue-step1');
    if (continue1) continue1.disabled = false;

    // Update viewer (real 3D with NGL)
    const pdbId = proteinCard.dataset.pdbId;
    const filename = proteinCard.dataset.filename || '';
    this.loadProteinViewerNGL({ proteinId, pdbId, filename });
        
        // Update molecule information
        this.updateMoleculeInfo(proteinCard);

        console.log(`Selected protein: ${proteinId}`);
    }

    async loadProteinViewerNGL({ proteinId, pdbId, filename }) {
        if (!this.nglStage) {
            this.initNGLStage();
            if (!this.nglStage) return;
        }
        try {
            this.nglStage.removeAllComponents();
            this.currentLigandComp = null;
            // Try load from uploads using provided filename first, then PDBID.pdb
            const localUrl = filename ? `/uploads/${filename}` : `/uploads/${pdbId}.pdb`;
            let comp;
            try {
                comp = await this.nglStage.loadFile(localUrl, { ext: 'pdb' });
            } catch (_) {
                // Fallback to RCSB
                comp = await this.nglStage.loadFile(`rcsb://${pdbId}`);
            }
            this.nglComp = comp;
            this.reprs.cartoon = comp.addRepresentation('cartoon', { colorScheme: 'element' });
            this.reprs.ligand = comp.addRepresentation('ball+stick', { sele: 'ligand', colorScheme: 'element' });
            this.reprs.surface = comp.addRepresentation('surface', { colorScheme: 'uniform', opacity: 0.7 });
            // Default visibility: protein+ligand on, surface off
            if (this.reprs.surface && typeof this.reprs.surface.setVisibility === 'function') {
                this.reprs.surface.setVisibility(false);
            }
            this.syncUIToReprs();
            comp.autoView();
        } catch (e) {
            console.error('NGL load error:', e);
            this.showNotification('Failed to load 3D structure', 'info');
        }
    }

    applyRepresentation(mode) {
        // Toggle between main protein reps: cartoon or ball+stick for entire structure
        const showCartoon = mode === 'cartoon';
        if (this.reprs.cartoon && typeof this.reprs.cartoon.setVisibility === 'function') {
            this.reprs.cartoon.setVisibility(showCartoon);
        }
        // Add or reuse a full-structure ball+stick for protein backbone when selected
        if (this.nglStage) {
            // Build/replace a full-structure ball+stick if needed
            const comp = this.nglComp;
            if (comp) {
                // Ensure a protein ball+stick repr exists
                if (!this.reprs.proteinBallStick) {
                    this.reprs.proteinBallStick = comp.addRepresentation('ball+stick', { sele: 'protein', colorScheme: 'element' });
                }
                const showProteinBS = !showCartoon;
                if (this.reprs.proteinBallStick.setVisibility) {
                    this.reprs.proteinBallStick.setVisibility(showProteinBS);
                }
            }
        }
    }
    syncUIToReprs() {
        const reprMode = document.getElementById('representation-select')?.value;
        const toggle = id => document.getElementById(id)?.checked;

        this.reprs.cartoon?.setVisibility(reprMode === 'cartoon' && toggle('show-protein'));
        this.reprs.proteinBallStick?.setVisibility(reprMode === 'ball+stick' && toggle('show-protein'));
        this.reprs.ligand?.setVisibility(toggle('show-ligands'));
        this.reprs.surface?.setVisibility(toggle('show-surface'));
    }

    ensureVerticesRepresentation(show) {
        if (!this.nglComp) return;
        if (show) {
            if (!this.reprs.vertices) {
                this.reprs.vertices = this.nglComp.addRepresentation('point', { sele: 'protein', pointSize: 0.2, colorScheme: 'element' });
            }
            if (this.reprs.vertices.setVisibility) this.reprs.vertices.setVisibility(true);
        } else if (this.reprs.vertices) {
            if (this.reprs.vertices.setVisibility) this.reprs.vertices.setVisibility(false);
        }
    }

    loadProteinViewer(proteinId) {
        const viewer = document.getElementById('protein-viewer');
        const noProteinState = viewer.querySelector('.no-protein-selected');
        
        if (noProteinState) {
            // Show loading state
            viewer.classList.add('loading');
            
            // Simulate loading time
            setTimeout(() => {
                noProteinState.style.display = 'none';
                viewer.classList.remove('loading');
                viewer.classList.add('protein-loaded');
                
                // Create mock 3D viewer content
                viewer.innerHTML = `
                    <div class="protein-viewer-content">
                        <div class="protein-structure">
                            <div class="protein-representation">
                                <div class="molecular-structure">
                                    <!-- Mock molecular visualization -->
                                    <div class="protein-backbone"></div>
                                    <div class="active-site"></div>
                                    <div class="binding-pocket"></div>
                                </div>
                            </div>
                        </div>
                        <div class="viewer-overlay">
                            <div class="structure-info">
                                <span class="structure-label">Protein Structure Loaded</span>
                            </div>
                        </div>
                    </div>
                `;
            }, 1500);
        }
    }

    updateMoleculeInfo(proteinCard) {
        const infoContainer = document.getElementById('molecule-info');
        const proteinName = proteinCard.querySelector('.protein-name').textContent;
        const pdbId = proteinCard.querySelector('.pdb-id').textContent;
        const description = proteinCard.querySelector('.protein-description').textContent;
        
        infoContainer.innerHTML = `
            <div class="molecule-details">
                <div class="detail-row mb-2">
                    <span class="detail-label">PDB ID:</span>
                    <span class="detail-value">${pdbId}</span>
                </div>
                <div class="detail-row mb-2">
                    <span class="detail-label">Description:</span>
                    <span class="detail-value">${description}</span>
                </div>
                <div class="detail-row mb-2">
                    <span class="detail-label">Resolution:</span>
                    <span class="detail-value">1.8 Å</span>
                </div>
                <div class="detail-row">
                    <span class="detail-label">Chain Count:</span>
                    <span class="detail-value">2</span>
                </div>
            </div>
        `;
    }

    searchProteins(query) {
        const proteinCards = document.querySelectorAll('.protein-card');
        const searchQuery = query.toLowerCase().trim();

        proteinCards.forEach(card => {
            const name = card.querySelector('.protein-name').textContent.toLowerCase();
            const pdbId = card.querySelector('.pdb-id').textContent.toLowerCase();
            const description = card.querySelector('.protein-description').textContent.toLowerCase();
            
            const matches = name.includes(searchQuery) || 
                          pdbId.includes(searchQuery) || 
                          description.includes(searchQuery);

            card.style.display = matches || !searchQuery ? 'block' : 'none';
        });
    }

    setWorkflowStep(stepNumber) {
        // Update workflow step visual state
        document.querySelectorAll('.workflow-step').forEach(step => {
            step.classList.remove('active');
        });
        
        document.querySelector(`.workflow-step[data-step="${stepNumber}"]`).classList.add('active');
        this.currentWorkflowStep = stepNumber;

        // Update UI based on step
        this.updateWorkflowUI(stepNumber);
        // Mirror left sidebar content to right sidebar for context
        this.mirrorRightSidebar(stepNumber);
    }

    mirrorRightSidebar(stepNumber) {
        const right = document.getElementById('sidebar-right');
        if (!right) return;
        // Decide which left panel to mirror
        const leftId = stepNumber === 1 ? 'sidebar-step-1' : (stepNumber === 2 ? 'sidebar-step-2' : (stepNumber === 3 ? 'sidebar-step-3' : 'sidebar-step-4'));
        const left = document.getElementById(leftId);
        if (!left) { right.innerHTML = ''; return; }
        // Clone node and strip ids from children to avoid duplicates
        const clone = left.cloneNode(true);
        // Remove any ids from clone subtree
        clone.querySelectorAll && clone.querySelectorAll('*').forEach(n => { if (n.id) n.removeAttribute('id'); });
        // Also remove id on root clone
        clone.removeAttribute && clone.removeAttribute('id');
        // Add a small header to indicate mirrored content
        right.innerHTML = '';
        const header = document.createElement('div');
        header.className = 'mirrored-header mb-2';
        header.innerHTML = `<h6 class="section-title">Preview</h6>`;
        right.appendChild(header);
        right.appendChild(clone);
    }

    updateWorkflowUI(stepNumber) {
        // Toggle sidebar sections
        const s1 = document.getElementById('sidebar-step-1');
        const s2 = document.getElementById('sidebar-step-2');
        const s3 = document.getElementById('sidebar-step-3');
        const s4 = document.getElementById('sidebar-step-4');
        if (s1 && s2 && s3 && s4) {
            s1.style.display = stepNumber === 1 ? 'block' : 'none';
            s2.style.display = stepNumber === 2 ? 'block' : 'none';
            s3.style.display = stepNumber === 3 ? 'block' : 'none';
            s4.style.display = stepNumber === 4 ? 'block' : 'none';
        }
        // Toggle right panels
        const pInfo = document.getElementById('panel-info');
        const pPockets = document.getElementById('panel-pockets');
        if (pInfo && pPockets) {
            pInfo.style.display = stepNumber === 2 ? 'none' : 'block';
            pPockets.style.display = stepNumber === 2 ? 'block' : 'none';
        }
        // Update workflow step active style
        document.querySelectorAll('.workflow-step').forEach(step => step.classList.remove('active'));
        const cur = document.querySelector(`.workflow-step[data-step="${stepNumber}"]`);
        if (cur) cur.classList.add('active');
    }

    handleDisplayToggle(toggle) {
        const isChecked = toggle.checked;
        const toggleId = toggle.id;
        
        console.log(`Display toggle ${toggleId}: ${isChecked ? 'ON' : 'OFF'}`);
        
        // Add visual feedback
        const label = toggle.closest('.form-check').querySelector('.form-check-label');
        if (label) {
            label.style.opacity = isChecked ? '1' : '0.6';
        }

        // Simulate updating 3D viewer based on toggles
        this.updateViewer3D();
    }

    updateViewer3D() {
        // Mock 3D viewer update based on current display settings
        const viewer = document.getElementById('protein-viewer');
        if (viewer.classList.contains('protein-loaded')) {
            viewer.style.filter = 'brightness(1.1)';
            setTimeout(() => {
                viewer.style.filter = '';
            }, 200);
        }
    }

    handleViewerControl(button) {
        const icon = button.querySelector('svg');
        const action = icon ? icon.getAttribute('data-feather') : '';
        
        // Add button press animation
        button.style.transform = 'scale(0.95)';
        setTimeout(() => {
            button.style.transform = '';
        }, 150);

        switch(action) {
            case 'rotate-ccw':
                console.log('Rotating left');ư
                break;
            case 'rotate-cw':
                console.log('Rotating right');
                break;
            case 'zoom-in':
                console.log('Zooming');
                break;
        }
    }

    updateRenderingQuality(quality) {
        // Map slider [0,1,2] -> NGL quality ['low','medium','high'] and stage sampleLevel
        const qMap = ['low', 'medium', 'high'];
        const q = qMap[Number(quality)] ?? 'medium';
        try {
            // Update stage sampling (antialias/supersampling)
            if (this.nglStage && typeof this.nglStage.setParameters === 'function') {
                this.nglStage.setParameters({
                    sampleLevel: Number(quality) 
                });
            }
            // Update all known representations
            Object.values(this.reprs).forEach(r => {
                if (r && typeof r.setParameters === 'function') {
                    r.setParameters({
                        quality: q 
                    });
                }
            });
        } catch (e) {
            console.warn('Failed to apply NGL quality:', e);
        }
    }

    updateSurfaceOpacity(opacity) {
        const alpha = Number(opacity) / 100;
        try {
            if (this.reprs.surface && typeof this.reprs.surface.setParameters === 'function') {
                this.reprs.surface.setParameters({
                    opacity: alpha 
                });
            }
        } catch (e) {
            console.warn('Failed to set surface opacity:', e);
        }
    }

    handleExportAction(button) {
        const buttonText = button.textContent.trim().toLowerCase();
        
        // Add loading state
        const originalText = button.innerHTML;
        button.innerHTML = '<i data-feather="loader" class="me-2"></i>Processing...';
        button.disabled = true;
        feather.replace();

        setTimeout(() => {
            button.innerHTML = originalText;
            button.disabled = false;
            feather.replace();

            if (buttonText.includes('png')) {
                this.showNotification('PNG export completed', 'success');
            } else if (buttonText.includes('download')) {
                this.showNotification('PDB file downloaded', 'success');
            } else if (buttonText.includes('share')) {
                this.showNotification('View link copied to clipboard', 'success');
            }
        }, 2000);
    }

    showNotification(message, type = 'info') {
        // Create notification element
        const notification = document.createElement('div');
        notification.className = `alert alert-${type === 'success' ? 'success' : 'info'} position-fixed`;
        notification.style.cssText = `
            top: 20px;
            right: 20px;
            z-index: 9999;
            min-width: 300px;
            border-radius: 8px;
            animation: slideIn 0.3s ease;
        `;
        notification.textContent = message;

        document.body.appendChild(notification);

        // Auto remove after 3 seconds
        setTimeout(() => {
            notification.style.animation = 'slideOut 0.3s ease';
            setTimeout(() => {
                if (notification.parentNode) {
                    notification.parentNode.removeChild(notification);
                }
            }, 300);
        }, 3000);
    }
}

// Initialize the application when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    window.moleculeAI = new MoleculeAI();
});

// Add CSS animations
const style = document.createElement('style');
style.textContent = `
    @keyframes slideIn {
        from { transform: translateX(100%); opacity: 0; }
        to { transform: translateX(0); opacity: 1; }
    }
    @keyframes slideOut {
        from { transform: translateX(0); opacity: 1; }
        to { transform: translateX(100%); opacity: 0; }
    }
    
    .detail-row {
        display: flex;
        justify-content: space-between;
        align-items: center;
    }
    
    .detail-label {
        font-size: 0.75rem;
        color: var(--text-muted);
        font-weight: 500;
    }
    
    .detail-value {
        font-size: 0.75rem;
        color: var(--text-primary);
        font-weight: 600;
    }
    
    .protein-viewer-content {
        position: relative;
        width: 100%;
        height: 100%;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    .molecular-structure {
        position: relative;
        width: 200px;
        height: 200px;
        border-radius: 50%;
        background: radial-gradient(circle, rgba(74, 144, 226, 0.3) 0%, rgba(80, 200, 163, 0.2) 50%, transparent 70%);
        animation: moleculeRotate 10s linear infinite;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    
    .protein-backbone {
        position: absolute;
        width: 80%;
        height: 80%;
        border: 2px solid rgba(74, 144, 226, 0.6);
        border-radius: 50%;
        border-style: dashed;
    }
    
    .active-site {
        position: absolute;
        width: 40%;
        height: 40%;
        background: radial-gradient(circle, rgba(80, 200, 163, 0.8) 0%, transparent 70%);
        border-radius: 50%;
        animation: pulse 2s ease-in-out infinite;
    }
    
    .binding-pocket {
        position: absolute;
        width: 20px;
        height: 20px;
        background-color: rgba(245, 158, 11, 0.8);
        border-radius: 50%;
        top: 30%;
        left: 60%;
        animation: glow 1.5s ease-in-out infinite alternate;
    }
    
    .viewer-overlay {
        position: absolute;
        bottom: 20px;
        left: 20px;
    }
    
    .structure-label {
        background-color: rgba(0, 0, 0, 0.7);
        color: var(--text-primary);
        padding: 0.5rem 1rem;
        border-radius: 6px;
        font-size: 0.875rem;
        font-weight: 500;
    }
    
    @keyframes moleculeRotate {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
    }
    
    @keyframes pulse {
        0%, 100% { opacity: 0.6; transform: scale(1); }
        50% { opacity: 1; transform: scale(1.1); }
    }
    
    @keyframes glow {
        0% { box-shadow: 0 0 5px rgba(245, 158, 11, 0.5); }
        100% { box-shadow: 0 0 20px rgba(245, 158, 11, 0.8), 0 0 30px rgba(245, 158, 11, 0.4); }
    }
`;
document.head.appendChild(style);

// --- Step 2: Pocket detection helpers ---
// auto-detect removed; manual pocket creation used instead

MoleculeAI.prototype.renderPocketsList = function() {
    const list = document.getElementById('pockets-list');
    const count = document.getElementById('pocket-count');
    const contBtn = document.getElementById('continue-step2');
    if (!list || !count) return;
    if (this.detectedPockets.length === 0) {
        list.innerHTML = '<p class="text-muted text-center py-4">No pockets detected</p>';
        count.textContent = '0 found';
        if (contBtn) contBtn.disabled = true;
        return;
    }
    count.textContent = `${this.detectedPockets.length} found`;
    list.innerHTML = this.detectedPockets.map((p, i) => `
        <div class="pocket-item ${i===0?'selected':''}" data-pocket-id="${p.id}">
            <div class="pocket-header"><div class="pocket-name">${p.name}</div><div class="pocket-score">${(p.druggability_score*100).toFixed(0)}%</div></div>
            <div class="pocket-info"><span class="pocket-volume">${p.volume}</span><span class="pocket-confidence">Conf: ${(p.confidence*100).toFixed(0)}%</span></div>
        </div>
    `).join('');
    list.querySelectorAll('.pocket-item').forEach(item => item.addEventListener('click', (e) => {
        list.querySelectorAll('.pocket-item').forEach(i => i.classList.remove('selected'));
        const el = e.currentTarget; el.classList.add('selected');
        const pid = el.dataset.pocketId; this.selectedPocket = this.detectedPockets.find(x => x.id === pid) || null;
        this.renderPocketProperties();
        if (contBtn) contBtn.disabled = !this.selectedPocket;
    }));
    // Auto-select first
    this.selectedPocket = this.detectedPockets[0];
    this.renderPocketProperties();
    if (contBtn) contBtn.disabled = false;
};

MoleculeAI.prototype.renderPocketProperties = function() {
    const props = document.getElementById('pocket-properties');
    if (!props) return;
    if (!this.selectedPocket) {
        props.innerHTML = '<p class="text-muted text-center py-4">Select a pocket to view properties</p>';
        return;
    }
    const p = this.selectedPocket;
    props.innerHTML = `
        <div class="pocket-details">
            <div class="property-row"><span class="property-label">Volume:</span><span class="property-value">${p.volume}</span></div>
            <div class="property-row"><span class="property-label">Radius:</span><span class="property-value">${p.radius} Å</span></div>
            <div class="property-row"><span class="property-label">Center:</span><span class="property-value">${p.center.x.toFixed(1)}, ${p.center.y.toFixed(1)}, ${p.center.z.toFixed(1)}</span></div>
            <div class="property-row"><span class="property-label">Druggability:</span><span class="property-value">${(p.druggability_score*100).toFixed(0)}%</span></div>
            <div class="property-row"><span class="property-label">Confidence:</span><span class="property-value">${(p.confidence*100).toFixed(0)}%</span></div>
            <div class="residues-section mt-2"><div class="property-label mb-1">Key Residues:</div><div class="residue-tags">${p.residues.map(r=>`<span class="residue-tag">${r}</span>`).join('')}</div></div>
        </div>
    `;
};

MoleculeAI.prototype.renderPocketOverlay = function() {
    // Optional: could add NGL shape objects, but we keep UI simple here
};

// --- Step 3: Sampling form helpers ---
MoleculeAI.prototype.resetSamplingDefaults = function() {
    const v = (id,val)=>{const el=document.getElementById(id); if(el){ if(el.type==='checkbox'){el.checked=!!val;} else {el.value=val;}}};
    v('n-samples',10); const nsv=document.getElementById('n-samples-value'); if(nsv) nsv.textContent='10';
    // Ensure slider max reflects new UI limit (100)
    const nsEl = document.getElementById('n-samples'); if (nsEl) nsEl.max = '100';
    v('legend-modes',20); v('model','Conditional model (Binding MOAD)');
    v('timesteps',100); const tsv=document.getElementById('timesteps-value'); if(tsv) tsv.textContent='100';
    v('resamplings',1); v('jump-length',1);
    const k=document.getElementById('keep-all-fragments'); if(k) k.checked=false;
    const s=document.getElementById('sanitize'); if(s) s.checked=true;
    const r=document.getElementById('relax'); if(r) r.checked=true;
};

MoleculeAI.prototype.submitSamplingForm = async function(e) {
    e.preventDefault();
    const form = e.currentTarget;
    const submitBtn = document.getElementById('run-generation');
    const original = submitBtn ? submitBtn.innerHTML : '';
    if (submitBtn) { submitBtn.innerHTML = '<i data-feather="loader" class="me-2"></i>Running...'; submitBtn.disabled = true; feather.replace(); }
    try {
        const fd = new FormData(form);
        const res = await fetch('/save_sampling_config', { method: 'POST', body: fd });
        const data = await res.json();
        if (data.success) {
            this.showNotification('Generation started successfully!', 'success');
            this.startProgressTracking();
        } else {
            throw new Error('Failed to start generation');
        }
    } catch (err) {
        this.showNotification('Error starting generation. Please try again.', 'info');
        if (submitBtn) { submitBtn.innerHTML = original; submitBtn.disabled = false; feather.replace(); }
    }
};

MoleculeAI.prototype.startProgressTracking = function() {
    const section = document.getElementById('progress-section');
    const bar = document.getElementById('progress-bar');
    const label = document.getElementById('progress-label');
    const pct = document.getElementById('progress-percentage');
    const mol = document.getElementById('molecules-count');
    const tEl = document.getElementById('time-elapsed');
    const tRem = document.getElementById('time-remaining');
    if (!section || !bar) return;
    section.style.display = 'block';
    let start = Date.now();

    // Kick off backend job using current selection
    // Prefer explicit ligand selection if available, otherwise fall back to manual text
    const selectedLigFile = sessionStorage.getItem('selectedLigandFilename') || '';
    const manualRefInput = document.getElementById('manual-ref');
    const manualRef = manualRefInput ? manualRefInput.value.trim() : '';
    const ref = selectedLigFile || manualRef;
    const pdbId = sessionStorage.getItem('selectedPdbId') || '';
    const filename = sessionStorage.getItem('selectedFilename') || '';
    const form = document.getElementById('sampling-form');
    const fd = new FormData(form);
    fd.append('ref_ligand', ref || '');
    fd.append('pdb_id', pdbId);
    fd.append('filename', filename);

    const submitBtn = document.getElementById('run-generation');
    const original = submitBtn ? submitBtn.innerHTML : '';

    fetch('/diffsbdd/start', { method: 'POST', body: fd }).then(r=>r.json()).then(data=>{
        if (!data.success) throw new Error(data.error || 'Failed to start DiffSBDD');
        const jobId = data.job_id;
        if (label) label.textContent = 'Initializing...';
        const poll = setInterval(async () => {
            try {
                const res = await fetch(`/diffsbdd/status/${jobId}`);
                const js = await res.json();
                if (!js.success) throw new Error('Status error');
                const st = js.status;
                const progress = st.progress ?? 0;
                bar.style.width = `${progress}%`;
                if (pct) pct.textContent = `${progress}%`;
                if (mol) mol.textContent = st.molecules ?? 0;
                const elapsed = Math.floor((Date.now()-start)/1000);
                if (tEl) tEl.textContent = `${elapsed}s`;
                if (label) label.textContent = st.message || st.state;
                if (progress >= 100 || st.state === 'success') {
                    clearInterval(poll);
                    bar.classList.remove('progress-bar-animated');
                    bar.classList.add('bg-success');
                    if (submitBtn) { submitBtn.innerHTML = original; submitBtn.disabled = false; feather.replace(); }
                    this.showNotification('Molecule generation completed!', 'success');
                    if (st.outfile_url) {
                        this.showNotification('Download: ' + st.outfile_url, 'success');
                    }
                    // Let user pick which molecule to view
                    this.showResultsPicker(jobId);
                }
                if (st.state === 'failed') {
                    clearInterval(poll);
                    if (submitBtn) { submitBtn.innerHTML = original; submitBtn.disabled = false; feather.replace(); }
                    this.showNotification('Generation failed. Check logs.', 'info');
                }
            } catch (e) {
                clearInterval(poll);
                if (submitBtn) { submitBtn.innerHTML = original; submitBtn.disabled = false; feather.replace(); }
                this.showNotification('Error polling status', 'info');
            }
        }, 2000);
    }).catch(err => {
        if (submitBtn) { submitBtn.innerHTML = original; submitBtn.disabled = false; feather.replace(); }
        this.showNotification(err.message || 'Error starting DiffSBDD', 'info');
    });
};

MoleculeAI.prototype.showResultsPicker = async function(jobId){
    try {
        const res = await fetch(`/diffsbdd/result/${jobId}/list`);
        const data = await res.json();
        if (!data.success) throw new Error('Failed to list results');
        const mols = data.molecules || [];
        if (mols.length === 0) return;
        // Persist jobId for future interactions
        this.lastJobId = jobId;
        // Populate persistent results list in step 4 panel
        const list = document.getElementById('results-list');
        if (list) {
            list.innerHTML = mols.map(m=>
                `<button data-idx="${m.index}" class="btn btn-outline-primary btn-sm w-100 mb-1 text-start">#${m.index+1} - ${m.title || 'molecule'}</button>`
            ).join('');
            list.querySelectorAll('button[data-idx]').forEach(btn=>{
                btn.addEventListener('click', async ()=>{
                    const idx = parseInt(btn.getAttribute('data-idx'));
                    await this.overlayGeneratedLigand(jobId, idx);
                    this.setWorkflowStep(4);
                });
            });
        }
        const dl = document.getElementById('download-sdf-link');
        if (dl) {
            dl.href = `/diffsbdd/result/${jobId}/sdf`;
            dl.classList.remove('d-none');
        }
        // Navigate to results panel
        this.setWorkflowStep(4);
    } catch (e) {
        this.showNotification('Failed to load results list', 'info');
    }
};

MoleculeAI.prototype.loadGeneratedPdbIntoViewer = async function(jobId, index){
    try {
        if (!this.nglStage) this.initNGLStage();
        if (!this.nglStage) return;
        const url = `/diffsbdd/result/${jobId}/pdb/${index}`;
        this.nglStage.removeAllComponents();
        const comp = await this.nglStage.loadFile(url, { ext: 'pdb' });
        this.nglComp = comp;
        this.reprs.cartoon = comp.addRepresentation('cartoon', { colorScheme: 'element' });
        this.reprs.ligand = comp.addRepresentation('ball+stick', { sele: 'ligand', colorScheme: 'element' });
        this.reprs.surface = comp.addRepresentation('surface', { colorScheme: 'uniform', opacity: 0.7 });
        if (this.reprs.surface && typeof this.reprs.surface.setVisibility === 'function') {
            this.reprs.surface.setVisibility(false);
        }
        this.syncUIToReprs();
        comp.autoView();
    } catch (e) {
        console.error(e);
        this.showNotification('Failed to load generated molecule', 'info');
    }
};

// Overlay the selected generated ligand on top of the currently selected protein
MoleculeAI.prototype.overlayGeneratedLigand = async function(jobId, index){
    try {
        if (!this.nglStage) this.initNGLStage();
        if (!this.nglStage) return;
        const pdbId = sessionStorage.getItem('selectedPdbId') || '';
        const filename = sessionStorage.getItem('selectedFilename') || '';
        // If no protein loaded yet, load it first
        if (!this.nglComp) {
            await this.loadProteinViewerNGL({ proteinId: this.selectedProtein || '', pdbId, filename });
        }
        if (!this.nglComp) return;
    // Remove previously loaded ligand if any
    try { if (this.currentLigandComp) { this.nglStage.removeComponent(this.currentLigandComp); } } catch(_) {}
    this.currentLigandComp = null;
    // Load ligand as a separate component and add representation
        const ligUrl = `/diffsbdd/result/${jobId}/pdb/${index}`;
    const ligComp = await this.nglStage.loadFile(ligUrl, { ext: 'pdb' });
    this.currentLigandComp = ligComp;
    ligComp.addRepresentation('ball+stick', { sele: 'ligand', colorScheme: 'element' });
    // Focus view to include both protein and ligand
    this.nglStage.autoView();
    } catch (e) {
        console.error(e);
        this.showNotification('Failed to overlay ligand with protein', 'info');
    }
};
