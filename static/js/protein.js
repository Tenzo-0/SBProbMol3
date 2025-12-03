// JS helpers for integrating DiffSBDD results with the ligand property table.

let currentJobId = null;

// Called from existing DiffSBDD workflow when a job is started
export function setCurrentDiffsbddJob(jobId) {
	currentJobId = jobId;
}

// Fetch scores for current job and render into the table
export async function refreshLigandScores() {
	if (!currentJobId) return;
	try {
		const resp = await fetch(`/diffsbdd/result/${currentJobId}/scores`);
		if (!resp.ok) return;
		const data = await resp.json();
		if (!data.success || !Array.isArray(data.scores)) return;
		const tbody = document.querySelector('#ligand-scores-table tbody');
		if (!tbody) return;
		tbody.innerHTML = '';
		data.scores.forEach((s, i) => {
			const tr = document.createElement('tr');
			const fmt = (v, digits = 3) => (v === null || v === undefined ? '-' : Number(v).toFixed(digits));
			tr.innerHTML = `
				<th scope="row">${i + 1}</th>
				<td>${s.title || `mol_${i + 1}`}</td>
				<td>${fmt(s.vina_score, 2)}</td>
				<td>${fmt(s.qed, 3)}</td>
				<td>${fmt(s.sa, 2)}</td>
			`;
			tbody.appendChild(tr);
		});
	} catch (e) {
		console.error('Failed to refresh ligand scores', e);
	}
}

// Optional: poll for updates while generation is running
export function startLigandScorePolling(intervalMs = 5000) {
	if (!currentJobId) return;
	refreshLigandScores();
	return setInterval(refreshLigandScores, intervalMs);
}

