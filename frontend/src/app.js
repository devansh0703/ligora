/**
 * Ligora App: binds the static index.html layout to real behavior.
 *
 * - 3D viewer: 3Dmol.js rendering the actual structure file delivered by
 *   the backend (mmCIF/PDB content; no re-fetching from third parties here).
 * - Panels: ligand card, contacts, evidence, jobs, notes - populated only
 *   from real backend responses.
 * - V2 surfaces: batch analysis, water network, 2D ligand editor,
 *   scripting console and pose comparison - all wired to real backend
 *   implementations (PLIP, CCD classification, RDKit, live APIs).
 */

import $3Dmol from '3dmol'
import { open, save } from '@tauri-apps/plugin-dialog'
import { invoke } from '@tauri-apps/api/core'
import { listen } from '@tauri-apps/api/event'

const CHAIN_COLORS = [
  '#4e79ba', '#f28e2b', '#59a14f', '#e15759', '#76b7b2',
  '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac',
]

const ELEMENT_COLORS = {
  H: '#e6e6e6', C: '#4ecdc4', N: '#58a6ff', O: '#ff6b6b',
  S: '#edc948', P: '#f28e2b', F: '#9ce0e0', CL: '#59a14f', BR: '#a371f7',
}

function escapeHtml(s) {
  return String(s ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;')
}

export function createApp({ sendCommand }) {
  const state = {
    structure: null,
    selectedLigandId: null,
    selectedLigand: null,
    contacts: [],
    representation: 'cartoon',
    colorScheme: 'chain',
    measureMode: false,
    picks: [],
  }

  let viewer = null
  let viewerElement = null

  const el = (id) => document.getElementById(id)

  // ------------------------------------------------------------------
  // Status toast
  // ------------------------------------------------------------------

  function showStatus(message, isError = false) {
    let toast = el('status-toast')
    if (!toast) {
      toast = document.createElement('div')
      toast.id = 'status-toast'
      document.body.appendChild(toast)
    }
    toast.textContent = message
    toast.className = isError ? 'status-toast error' : 'status-toast'
    toast.style.display = 'block'
    clearTimeout(toast._timer)
    toast._timer = setTimeout(() => { toast.style.display = 'none' }, 4000)
  }

  function setBusy(busy, text = 'Loading...') {
    const overlay = el('loading-overlay')
    if (overlay) {
      overlay.style.display = busy ? 'flex' : 'none'
      const label = overlay.querySelector('.loading-text')
      if (label) label.textContent = text
    }
  }

  // ------------------------------------------------------------------
  // 3D viewer (3Dmol.js)
  // ------------------------------------------------------------------

  function initViewer() {
    viewerElement = el('scene-container')
    if (!viewerElement || viewer) return
    viewer = $3Dmol.createViewer(viewerElement, {
      backgroundColor: '#1a1a2e',
      callback: onAtomPicked,
    })
    const ro = new ResizeObserver(() => {
      if (viewer) viewer.resize()
    })
    ro.observe(viewerElement)
  }

  async function loadStructureIntoViewer() {
    if (!viewer) return
    const file = await sendCommand('get_structure_file', {})
    viewer.clear()
    viewer.addModel(file.content, file.format === 'pdb' ? 'pdb' : 'cif')
    applyView()
    viewer.zoomTo()
    viewer.render()
  }

  function proteinStyle() {
    const clickable = state.measureMode ? { clickable: true } : {}
    const sel = { not: { elem: 'H' } }
    switch (state.representation) {
      case 'stick':
        return [{ sel, style: { stick: { radius: 0.25, ...clickable } } }]
      case 'sphere':
        return [{ sel, style: { sphere: { scale: 0.4, ...clickable } } }]
      case 'line':
        return [{ sel, style: { line: { ...clickable } } }]
      case 'licorice':
        return [{ sel, style: { stick: { radius: 0.12, ...clickable } } }]
      case 'surface':
        return [{ sel, style: { cartoon: { color: 'spectrum', ...clickable } } }]
      case 'cartoon':
      default:
        return [{ sel, style: { cartoon: { color: 'spectrum', ...clickable } } }]
    }
  }

  function chainStyle() {
    // Real per-chain colors from the backend structure's chain list.
    const clickable = state.measureMode ? { clickable: true } : {}
    const chains = (state.structure?.chains || []).filter(c => c.is_polymer)
    const styles = chains.map((chain, i) => ({
      sel: { chain: chain.id, not: { elem: 'H' } },
      style: {
        cartoon: { color: CHAIN_COLORS[i % CHAIN_COLORS.length], ...clickable },
        stick: { color: CHAIN_COLORS[i % CHAIN_COLORS.length], ...clickable },
      },
    }))
    if (styles.length === 0) {
      return [{ sel: { not: { elem: 'H' } }, style: { cartoon: { color: '#4e79ba', ...clickable } } }]
    }
    return styles
  }

  function applyView() {
    if (!viewer) return
    viewer.removeAllSurfaces()
    const styles = state.colorScheme === 'chain' ? chainStyle() : proteinStyle()
    for (const { sel, style } of styles) {
      viewer.setStyle(sel, style)
    }
    if (state.representation === 'surface') {
      const type = $3Dmol.SurfaceType ? $3Dmol.SurfaceType.VDW : 'VDW'
      viewer.addSurface(type, { opacity: 0.75 }, { not: { elem: 'H' } })
    }
    // Ligand always visible as sticks + spheres, highlighted when selected.
    const ligSel = state.selectedLigand
      ? { resn: state.selectedLigand.residue_name }
      : null
    if (ligSel) {
      viewer.addStyle(ligSel, {
        stick: { radius: 0.3, color: '#4ecdc4' },
        sphere: { scale: 0.35, color: '#4ecdc4' },
      })
    }
    // Contacted protein residues highlighted from real PLIP results.
    for (const c of state.contacts) {
      viewer.addStyle(
        { chain: c.protein_chain_id, resi: c.protein_residue_id },
        { stick: { radius: 0.28, color: '#ff6b6b' } })
    }
    viewer.render()
  }

  // ------------------------------------------------------------------
  // Measurement (real coordinates via backend math)
  // ------------------------------------------------------------------

  function onAtomPicked(atom) {
    if (!state.measureMode || !atom) return
    state.picks.push({
      chain: atom.chain,
      residue_id: atom.resi,
      residue_name: atom.resn,
      name: atom.name,
      x: atom.x, y: atom.y, z: atom.z,
    })
    if (state.picks.length > 3) state.picks = state.picks.slice(-3)
    updateMeasureOverlay()
    if (state.picks.length === 2) {
      measure('measure_distance', 'distance')
    } else if (state.picks.length === 3) {
      measure('measure_angle', 'angle')
    }
  }

  async function measure(command, field) {
    try {
      const refs = state.picks.map(p => ({
        chain: p.chain, residue_id: p.residue_id, name: p.name,
      }))
      const payload = command === 'measure_distance'
        ? { atom1: refs[0], atom2: refs[1] }
        : { atom1: refs[0], atom2: refs[1], atom3: refs[2] }
      const result = await sendCommand(command, payload)
      updateMeasureOverlay(result[field])
    } catch (e) {
      showStatus(`Measurement failed: ${e.message}`, true)
    }
  }

  function updateMeasureOverlay(value) {
    let overlay = el('measure-overlay')
    if (!overlay) {
      overlay = document.createElement('div')
      overlay.id = 'measure-overlay'
      viewerElement?.appendChild(overlay)
    }
    const labels = state.picks.map(p =>
      `${escapeHtml(p.chain || '?')}:${p.residue_id ?? '?'}:${escapeHtml(p.name || '?')}`)
    overlay.innerHTML = `
      <div class="measure-points">${labels.map(l => `<span>${l}</span>`).join(' → ')}</div>
      ${value != null
        ? `<div class="measure-result">${Number(value).toFixed(3)}
           ${state.picks.length === 3 ? '°' : ' Å'}</div>`
        : '<div class="measure-result">pick 1 more atom</div>'}
    `
  }

  // ------------------------------------------------------------------
  // Panels
  // ------------------------------------------------------------------

  function updateLigandSelect() {
    const select = el('ligand-select')
    if (!select) return
    const ligands = state.structure?.ligands || []
    const seen = new Set()
    select.innerHTML = ''
    for (const lig of ligands) {
      const key = `${lig.residue_name}:${lig.id}`
      if (seen.has(lig.residue_name)) continue
      seen.add(lig.residue_name)
      const opt = document.createElement('option')
      opt.value = lig.id
      opt.textContent = `${lig.residue_name} (${lig.name || ''})`
      select.appendChild(opt)
    }
    if (state.selectedLigandId) select.value = state.selectedLigandId
  }

  function updateLigandCard(ligand) {
    const card = el('ligand-card')
    if (!card || !ligand) return
    card.innerHTML = `
      <div class="ligand-header">
        <span class="ligand-name">${escapeHtml(ligand.name || ligand.residue_name)}</span>
        <span class="ligand-badge">${escapeHtml(ligand.classification_hint || 'unclassified')}</span>
      </div>
      <div class="ligand-properties">
        ${ligandRow('Formula', ligand.formula)}
        ${ligandRow('MW', ligand.molecular_weight != null ? Number(ligand.molecular_weight).toFixed(2) : null)}
        ${ligandRow('SMILES', ligand.smiles)}
        ${ligandRow('InChI Key', ligand.inchi_key)}
        ${ligandRow('IUPAC', ligand.iupac_name)}
        ${ligandRow('PubChem CID', ligand.pubchem_cid ? String(ligand.pubchem_cid) : null)}
        ${ligandRow('ChEMBL', ligand.chembl_id)}
        ${ligandRow('Status', ligand.resolution_status)}
      </div>`
  }

  function ligandRow(label, value) {
    return `<div class="ligand-property">
      <div class="ligand-property-label">${escapeHtml(label)}</div>
      <div class="ligand-property-value">${value ? escapeHtml(value) : '—'}</div>
    </div>`
  }

  function updateContactsPanel(contacts, plipAvailable) {
    const pane = el('contacts-table')
    if (!pane) return
    if (!contacts || contacts.length === 0) {
      pane.innerHTML = `<div class="contacts-placeholder">${
        plipAvailable === false
          ? 'PLIP is not installed - contact analysis unavailable'
          : 'No contacts found'}</div>`
      return
    }
    let html = '<table><thead><tr><th>Type</th><th>Protein</th><th>Lig atom</th><th>Å</th></tr></thead><tbody>'
    for (const c of contacts.slice(0, 100)) {
      html += `<tr class="contact-row" data-chain="${escapeHtml(c.protein_chain_id)}" data-resi="${escapeHtml(String(c.protein_residue_id))}">
        <td><span class="contact-badge">${escapeHtml(c.contact_type || '')}</span></td>
        <td>${escapeHtml(c.protein_residue_name)}${escapeHtml(String(c.protein_residue_id))} (${escapeHtml(c.protein_chain_id)})</td>
        <td>${escapeHtml(c.ligand_atom || '')}</td>
        <td>${c.distance != null ? Number(c.distance).toFixed(2) : '—'}</td>
      </tr>`
    }
    html += '</tbody></table>'
    if (contacts.length > 100) {
      html += `<div class="contacts-more">... and ${contacts.length - 100} more</div>`
    }
    pane.innerHTML = html
    pane.querySelectorAll('.contact-row').forEach(row => {
      row.addEventListener('click', () => {
        if (!viewer) return
        viewer.zoomTo({
          chain: row.dataset.chain,
          resi: parseInt(row.dataset.resi, 10),
        })
        viewer.render()
      })
    })
  }

  function updateEvidencePanel(evidence) {
    const pane = el('evidence-pane')
    if (!pane) return
    if (!evidence || evidence.length === 0) {
      pane.innerHTML = '<div class="evidence-placeholder">No enrichment evidence available</div>'
      return
    }
    pane.innerHTML = evidence.map(e => {
      const val = typeof e.value === 'object' && e.value !== null
        ? JSON.stringify(e.value) : String(e.value ?? '')
      return `<div class="evidence-item">
        <div class="evidence-source">${escapeHtml(e.source)}</div>
        <div class="evidence-value">${escapeHtml(val)}</div>
        ${e.url ? `<a class="evidence-url" href="${escapeHtml(e.url)}" target="_blank" rel="noopener">${escapeHtml(e.url)}</a>` : ''}
      </div>`
    }).join('')
  }

  function updateJobsPanel(status) {
    const pane = el('jobs-pane')
    if (!pane) return
    if (!status) {
      pane.innerHTML = '<div class="jobs-placeholder">No running jobs</div>'
      return
    }
    let html = ''
    const engines = status.engines || {}
    html += '<div class="jobs-section"><div class="jobs-section-title">Engines</div><div class="jobs-list">'
    for (const [name, info] of Object.entries(engines)) {
      const ok = info && info.available === true
      html += `<div class="job-item">
        <span class="job-status ${ok ? 'completed' : 'failed'}"></span>
        <span>${escapeHtml(name)}</span>
        <span class="job-detail">${ok ? 'available' : 'not installed'}</span>
      </div>`
    }
    html += '</div></div>'
    const sources = status.data_sources || {}
    const sourceNames = Object.keys(sources)
    if (sourceNames.length > 0) {
      html += '<div class="jobs-section"><div class="jobs-section-title">Data sources</div><div class="jobs-list">'
      for (const [name, st] of Object.entries(sources)) {
        html += `<div class="job-item">
          <span class="job-status ${st === 'ok' ? 'completed' : 'failed'}"></span>
          <span>${escapeHtml(name)}</span>
          <span class="job-detail">${escapeHtml(st)}</span>
        </div>`
      }
      html += '</div></div>'
    }
    pane.innerHTML = html
  }

  function showJobResult(job) {
    const pane = el('jobs-pane')
    if (!pane || !job) return
    const poses = job.result?.poses || []
    if (poses.length > 0) {
      const rows = poses.map(p =>
        `<div class="job-item"><span class="job-status completed"></span>
         <span>pose ${p.pose_id}</span>
         <span class="job-detail">${Number(p.affinity).toFixed(2)} kcal/mol</span></div>`).join('')
      pane.insertAdjacentHTML('beforeend',
        `<div class="jobs-section"><div class="jobs-section-title">Docking poses</div><div class="jobs-list">${rows}</div></div>`)
    } else if (job.error) {
      showStatus(`Job failed: ${job.error}`, true)
    }
  }

  // ------------------------------------------------------------------
  // Actions
  // ------------------------------------------------------------------

  async function openLocalFile() {
    const selected = await open({
      multiple: false,
      filters: [{ name: 'Structure files', extensions: ['cif', 'mmcif', 'cif.gz', 'pdb', 'ent'] }],
    })
    if (!selected) return
    setBusy(true, 'Opening file...')
    try {
      const data = await sendCommand('open_local_file', { file_path: selected })
      await onStructureLoaded(data)
    } catch (e) {
      showStatus(`Failed to open file: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  /** Real in-app modal prompt (no browser prompt dialogs). */
  function appPrompt(title, placeholder) {
    return new Promise((resolve) => {
      let modal = el('app-modal')
      if (!modal) {
        modal = document.createElement('div')
        modal.id = 'app-modal'
        modal.className = 'modal-backdrop'
        document.body.appendChild(modal)
      }
      modal.innerHTML = `
        <div class="modal">
          <div class="modal-title">${escapeHtml(title)}</div>
          <input class="input modal-input" id="modal-input"
              placeholder="${escapeHtml(placeholder || '')}" autocomplete="off">
          <div class="modal-actions">
            <button class="btn btn-secondary" id="modal-cancel">Cancel</button>
            <button class="btn btn-primary" id="modal-ok">Open</button>
          </div>
        </div>`
      modal.style.display = 'flex'
      const input = el('modal-input')
      input.focus()
      const close = (value) => {
        modal.style.display = 'none'
        modal.innerHTML = ''
        document.removeEventListener('keydown', onKey)
        resolve(value)
      }
      const onKey = (e) => {
        if (e.key === 'Enter') close(input.value.trim())
        if (e.key === 'Escape') close(null)
      }
      document.addEventListener('keydown', onKey)
      el('modal-cancel').addEventListener('click', () => close(null))
      el('modal-ok').addEventListener('click', () => close(input.value.trim()))
      modal.addEventListener('mousedown', (e) => {
        if (e.target === modal) close(null)
      })
    })
  }

  async function openPdbId() {
    const pdbId = await appPrompt('Open PDB ID from RCSB', 'e.g. 3W85')
    if (!pdbId) return
    if (!/^[0-9][A-Za-z0-9]{3}$/.test(pdbId)) {
      showStatus('Invalid PDB ID format (expect 4 characters, e.g. 3W85)', true)
      return
    }
    setBusy(true, `Fetching ${pdbId.toUpperCase()} from RCSB...`)
    try {
      const data = await sendCommand('open_pdb_id', { pdb_id: pdbId.trim().toUpperCase() })
      await onStructureLoaded(data)
    } catch (e) {
      showStatus(`Failed to open ${pdbId.toUpperCase()}: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  async function onStructureLoaded(data) {
    state.structure = data.structure
    state.contacts = []
    state.selectedLigand = null
    state.selectedLigandId = null
    state.picks = []
    updateLigandSelect()
    updateContactsPanel([], true)
    updateEvidencePanel([])
    await loadStructureIntoViewer()
    // Auto-select the first ligand that is not solvent/ion per its CCD class.
    const first = (state.structure.ligands || []).find(
      l => l.classification_hint && !['HETAS', 'HETAI'].includes(l.classification_hint))
      || (state.structure.ligands || [])[0]
    if (first) await selectLigand(first.id)
    showStatus(`Loaded ${state.structure.id}`)
  }

  async function selectLigand(ligandId) {
    state.selectedLigandId = ligandId
    try {
      const data = await sendCommand('select_ligand', { ligand_id: ligandId })
      state.selectedLigand = data.ligand
      updateLigandCard(data.ligand)
      updateLigandSelect()
      applyView()
    } catch (e) {
      showStatus(`Ligand selection failed: ${e.message}`, true)
    }
  }

  async function runAnalysis() {
    setBusy(true, 'Running PLIP contact analysis...')
    try {
      const data = await sendCommand('run_contact_analysis', {})
      state.contacts = data.contacts || []
      if (data.ligand) state.selectedLigand = data.ligand
      updateContactsPanel(data.contacts, data.plip_available)
      if (data.ligand) updateLigandCard(data.ligand)
      updateEvidencePanel(data.evidence)
      applyView()
      showStatus(`${data.contact_count} contacts found (PLIP)`)
    } catch (e) {
      showStatus(`Analysis failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  async function refreshStatus() {
    try {
      const status = await sendCommand('get_status', {})
      updateJobsPanel(status)
    } catch (e) {
      showStatus(`Status failed: ${e.message}`, true)
    }
  }

  async function exportArtifacts() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    try {
      const saved = await sendCommand('save_artifact', {})
      const dest = await save({
        defaultPath: `analysis_${state.structure.id}.json`,
        filters: [{ name: 'JSON', extensions: ['json'] }],
      })
      if (dest) {
        await invoke('copy_file', { src: saved.path, dest })
      }
      try {
        const sdf = await sendCommand('export_ligand_sdf', {})
        const sdfDest = await save({
          defaultPath: `ligand_${state.selectedLigand?.residue_name || 'ligand'}.sdf`,
          filters: [{ name: 'SDF', extensions: ['sdf'] }],
        })
        if (sdfDest) {
          await invoke('copy_file', { src: sdf.path, dest: sdfDest })
        }
      } catch (e) {
        showStatus(`SDF export skipped: ${e.message}`, true)
      }
      showStatus('Export complete')
    } catch (e) {
      showStatus(`Export failed: ${e.message}`, true)
    }
  }

  function toggleMeasure() {
    state.measureMode = !state.measureMode
    if (!state.measureMode) {
      state.picks = []
      updateMeasureOverlay(null)
      const overlay = el('measure-overlay')
      if (overlay) overlay.style.display = 'none'
    } else if (overlay) {
      overlay.style.display = 'block'
    }
    const btn = el('btn-measure')
    if (btn) btn.classList.toggle('active', state.measureMode)
    applyView()
  }

  // ==================================================================
  // V2: batch analysis
  // ==================================================================

  let lastBatchId = null

  async function runBatch() {
    const raw = el('batch-sources')?.value || ''
    const sources = raw.split(',').map(s => s.trim()).filter(Boolean)
    if (sources.length === 0) {
      showStatus('Enter at least one PDB ID or file path', true)
      return
    }
    const sourceType = sources.every(s => /^[0-9][a-z0-9]{3}$/i.test(s))
      ? 'pdb_id' : 'local'
    setBusy(true, 'Running batch analysis (PLIP + live enrichment)...')
    try {
      const added = await sendCommand('batch_add', { sources, source_type: sourceType })
      lastBatchId = added.batch_id
      const summary = await sendCommand('batch_run', { batch_id: added.batch_id })
      renderBatchResults(summary)
      showStatus(`Batch done: ${summary.succeeded}/${summary.total} succeeded`)
    } catch (e) {
      showStatus(`Batch failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  function renderBatchResults(summary) {
    const pane = el('batch-results')
    if (!pane) return
    let html = `<div class="v2-summary">${summary.succeeded} succeeded / ${summary.failed} failed of ${summary.total}</div>`
    html += '<table><thead><tr><th>Structure</th><th>Ligand</th><th>Contacts</th><th></th></tr></thead><tbody>'
    for (const r of summary.results || []) {
      html += `<tr>
        <td>${escapeHtml(String(r.structure_id || ''))}</td>
        <td>${escapeHtml(r.ligand_name || '—')}</td>
        <td>${r.contact_count != null ? String(r.contact_count) : '—'}</td>
        <td>${r.success
          ? '<span class="job-status completed"></span>'
          : `<span class="job-status failed"></span><span class="v2-error" title="${escapeHtml(r.error || '')}">failed</span>`}</td>
      </tr>`
    }
    html += '</tbody></table>'
    pane.innerHTML = html
  }

  async function exportBatch() {
    if (!lastBatchId) {
      showStatus('Run a batch first', true)
      return
    }
    try {
      const result = await sendCommand('batch_export', { batch_id: lastBatchId })
      const dest = await save({
        defaultPath: 'batch_results',
        filters: [{ name: 'Directory', extensions: ['*'] }],
      })
      if (dest) {
        await invoke('copy_file', { src: result.path, dest })
      }
      showStatus(`Batch exported: ${result.path}`)
    } catch (e) {
      showStatus(`Batch export failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2: water network
  // ==================================================================

  async function analyzeWaterNetwork() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    setBusy(true, 'Analyzing water network...')
    try {
      const data = await sendCommand('analyze_water_network', {})
      renderWaterNetwork(data)
      showStatus(`${data.water_count} waters, ${data.network.length} network edges`)
    } catch (e) {
      showStatus(`Water analysis failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  function renderWaterNetwork(data) {
    const pane = el('water-results')
    if (!pane) return
    let html = `<div class="v2-summary">${data.water_count} waters · ${data.network.length} water-water edges · ${data.water_contacts.length} water-protein contacts</div>`
    // Clusters
    if (data.clusters.length > 0) {
      html += `<div class="v2-section-title">Clusters</div><div class="v2-chips">`
      for (const c of data.clusters.slice(0, 20)) {
        html += `<span class="v2-chip">cluster of ${c.size}</span>`
      }
      html += '</div>'
    }
    // Nearest protein contacts (real distances from the structure)
    html += '<div class="v2-section-title">Water-protein contacts (nearest)</div>'
    html += '<table><thead><tr><th>Water</th><th>Protein atom</th><th>Å</th></tr></thead><tbody>'
    const contacts = [...data.water_contacts].sort((a, b) => a.distance - b.distance).slice(0, 40)
    for (const c of contacts) {
      html += `<tr class="contact-row" data-chain="${escapeHtml(c.protein.chain)}" data-resi="${escapeHtml(String(c.protein.residue_id))}">
        <td>${escapeHtml(c.water.chain)}:${escapeHtml(String(c.water.residue))}</td>
        <td>${escapeHtml(c.protein.residue)}${escapeHtml(String(c.protein.residue_id))}/${escapeHtml(c.protein.atom)} (${escapeHtml(c.protein.chain)})</td>
        <td>${Number(c.distance).toFixed(2)}</td>
      </tr>`
    }
    html += '</tbody></table>'
    pane.innerHTML = html
    pane.querySelectorAll('.contact-row').forEach(row => {
      row.addEventListener('click', () => {
        if (!viewer) return
        viewer.zoomTo({
          chain: row.dataset.chain,
          resi: parseInt(row.dataset.resi, 10),
        })
        viewer.render()
      })
    })
  }

  // ==================================================================
  // V2: 2D ligand editor
  // ==================================================================

  const editorState = { atoms: [], bonds: [], selection: [], dirty: false }

  async function editorLoad() {
    if (!state.structure) {
      showStatus('Open a structure and select a ligand first', true)
      return
    }
    try {
      const data = await sendCommand('get_ligand_2d', {})
      editorState.atoms = data.atoms || []
      editorState.bonds = data.bonds || []
      editorState.selection = []
      editorState.dirty = false
      drawEditorCanvas()
      el('editor-info').innerHTML =
        `<div class="v2-summary">${escapeHtml(data.residue_name || '')}: ${editorState.atoms.length} atoms, ${editorState.bonds.length} bonds (CCD bond orders)</div>`
      showStatus('Ligand loaded into 2D editor')
    } catch (e) {
      showStatus(`Editor load failed: ${e.message}`, true)
    }
  }

  function drawEditorCanvas() {
    const canvas = el('editor-canvas')
    if (!canvas) return
    const ctx = canvas.getContext('2d')
    ctx.clearRect(0, 0, canvas.width, canvas.height)
    const atoms = editorState.atoms
    if (atoms.length === 0) return
    // Real 2D coordinates from the editor state (x/y), auto-scaled to fit.
    const xs = atoms.map(a => a.x), ys = atoms.map(a => a.y)
    const minX = Math.min(...xs), maxX = Math.max(...xs)
    const minY = Math.min(...ys), maxY = Math.max(...ys)
    const pad = 28
    const scale = Math.min(
      (canvas.width - 2 * pad) / Math.max(maxX - minX, 0.001),
      (canvas.height - 2 * pad) / Math.max(maxY - minY, 0.001))
    const toPx = a => ({
      px: pad + (a.x - minX) * scale + (canvas.width - 2 * pad - (maxX - minX) * scale) / 2,
      py: canvas.height - (pad + (a.y - minY) * scale + (canvas.height - 2 * pad - (maxY - minY) * scale) / 2),
    })
    // Bonds (double/triple drawn as parallel lines)
    for (const b of editorState.bonds) {
      const a1 = atoms.find(a => a.id === b.from)
      const a2 = atoms.find(a => a.id === b.to)
      if (!a1 || !a2) continue
      const p1 = toPx(a1), p2 = toPx(a2)
      const order = b.order || 1
      ctx.strokeStyle = '#8b949e'
      for (let k = 0; k < order; k++) {
        const off = (k - (order - 1) / 2) * 3.2
        const dx = p2.px - p1.px, dy = p2.py - p1.py
        const len = Math.hypot(dx, dy) || 1
        ctx.beginPath()
        ctx.moveTo(p1.px - dy / len * off, p1.py + dx / len * off)
        ctx.lineTo(p2.px - dy / len * off, p2.py + dx / len * off)
        ctx.stroke()
      }
    }
    // Atoms
    for (const a of atoms) {
      const { px, py } = toPx(a)
      const selected = editorState.selection.includes(a.id)
      const color = ELEMENT_COLORS[(a.element || '').toUpperCase()] || '#c9d1d9'
      ctx.beginPath()
      ctx.arc(px, py, selected ? 8 : 6, 0, 2 * Math.PI)
      ctx.fillStyle = selected ? '#58a6ff' : color
      ctx.fill()
      if (a.element && a.element !== 'C') {
        ctx.fillStyle = '#0d1117'
        ctx.font = 'bold 8px sans-serif'
        ctx.textAlign = 'center'
        ctx.textBaseline = 'middle'
        ctx.fillText(a.element, px, py)
      }
    }
  }

  function editorPick(event) {
    const canvas = el('editor-canvas')
    if (!canvas || editorState.atoms.length === 0) return
    const rect = canvas.getBoundingClientRect()
    const mx = (event.clientX - rect.left) * (canvas.width / rect.width)
    const my = (event.clientY - rect.top) * (canvas.height / rect.height)
    // Find nearest atom in canvas space using the same transform as drawing.
    const atoms = editorState.atoms
    const xs = atoms.map(a => a.x), ys = atoms.map(a => a.y)
    const minX = Math.min(...xs), maxX = Math.max(...xs)
    const minY = Math.min(...ys), maxY = Math.max(...ys)
    const pad = 28
    const scale = Math.min(
      (canvas.width - 2 * pad) / Math.max(maxX - minX, 0.001),
      (canvas.height - 2 * pad) / Math.max(maxY - minY, 0.001))
    let best = null, bestD = 14 * 14
    for (const a of atoms) {
      const px = pad + (a.x - minX) * scale + (canvas.width - 2 * pad - (maxX - minX) * scale) / 2
      const py = canvas.height - (pad + (a.y - minY) * scale + (canvas.height - 2 * pad - (maxY - minY) * scale) / 2)
      const d = (px - mx) ** 2 + (py - my) ** 2
      if (d < bestD) { best = a; bestD = d }
    }
    if (!best) return
    // Two-step selection: first pick selects, second pick applies the
    // pending operation (bond/atom removal or bond creation).
    if (editorState.pendingOp) {
      const op = editorState.pendingOp
      editorState.pendingOp = null
      applyEditorOp(op, best.id)
      return
    }
    editorState.selection = [best.id]
    drawEditorCanvas()
    el('editor-info').innerHTML =
      `<div class="v2-summary">selected atom ${best.id} (${escapeHtml(best.element || '?')}) — pick an operation, then a second atom when prompted</div>`
  }

  async function applyEditorOp(op, secondAtomId) {
    try {
      let result
      if (op.kind === 'add_bond') {
        result = await sendCommand('editor_add_bond', {
          from_atom: op.atomId, to_atom: secondAtomId,
          order: op.order,
        })
      } else if (op.kind === 'remove_bond') {
        result = await sendCommand('editor_remove_bond', {
          from_atom: op.atomId, to_atom: secondAtomId,
        })
      }
      if (result?.state) {
        editorState.atoms = result.state.atoms
        editorState.bonds = result.state.bonds
        editorState.selection = []
        editorState.dirty = true
        drawEditorCanvas()
      }
    } catch (e) {
      showStatus(`Edit failed: ${e.message}`, true)
    }
  }

  async function editorAddAtom() {
    const element = el('editor-element')?.value || 'C'
    // Place the new atom at the centroid of the current layout (real 2D
    // coordinates; the user drags positions via update_position later).
    const atoms = editorState.atoms
    if (!atoms.length) {
      showStatus('Load a ligand first', true)
      return
    }
    const cx = atoms.reduce((s, a) => s + a.x, 0) / atoms.length
    const cy = atoms.reduce((s, a) => s + a.y, 0) / atoms.length
    try {
      const result = await sendCommand('editor_add_atom', {
        element, x: cx + 1.5, y: cy + 1.5,
      })
      editorState.atoms = result.state.atoms
      editorState.bonds = result.state.bonds
      editorState.dirty = true
      drawEditorCanvas()
      el('editor-info').innerHTML =
        `<div class="v2-summary">added ${escapeHtml(element)} atom ${result.atom_id}</div>`
    } catch (e) {
      showStatus(`Add atom failed: ${e.message}`, true)
    }
  }

  function editorRemoveAtom() {
    if (editorState.selection.length !== 1) {
      showStatus('Select one atom first (click it on the canvas)', true)
      return
    }
    sendCommand('editor_remove_atom', { atom_id: editorState.selection[0] })
      .then(result => {
        editorState.atoms = result.state.atoms
        editorState.bonds = result.state.bonds
        editorState.selection = []
        editorState.dirty = true
        drawEditorCanvas()
      })
      .catch(e => showStatus(`Remove atom failed: ${e.message}`, true))
  }

  function editorBeginOp(kind) {
    if (editorState.selection.length !== 1) {
      showStatus('Select one atom first (click it on the canvas)', true)
      return
    }
    const order = kind === 'add_bond'
      ? parseInt(el('editor-bond-order')?.value || '1', 10) : 1
    editorState.pendingOp = { kind, atomId: editorState.selection[0], order }
    showStatus('Now click the second atom')
  }

  async function editorExportSdf() {
    try {
      const result = await sendCommand('editor_export_sdf', {})
      const dest = await save({
        defaultPath: 'edited_ligand.sdf',
        filters: [{ name: 'SDF', extensions: ['sdf'] }],
      })
      if (dest) {
        await invoke('copy_file', { src: result.path, dest })
        showStatus('Edited ligand SDF saved')
      }
    } catch (e) {
      showStatus(`Editor export failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2: scripting console
  // ==================================================================

  async function runScript() {
    const code = el('script-input')?.value || ''
    if (!code.trim()) return
    try {
      const result = await sendCommand('run_script', { code })
      const out = el('script-output')
      if (out) {
        const printed = result.output || ''
        const vars = (result.variables || []).join(', ') || '(none)'
        out.textContent = `${printed || '(no output)'}\n# variables now defined: ${vars}`
      }
      showStatus('Script executed')
    } catch (e) {
      const out = el('script-output')
      if (out) out.textContent = `Error: ${e.message}`
      showStatus(`Script failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2: pose comparison
  // ==================================================================

  async function refreshComparisonJobs() {
    try {
      const data = await sendCommand('list_jobs', {})
      const completed = (data.jobs || []).filter(j => j.status === 'completed' && j.pose_count > 0)
      for (const id of ['compare-job-a', 'compare-job-b']) {
        const sel = el(id)
        if (!sel) continue
        const current = sel.value
        sel.innerHTML = '<option value="">Job...</option>' +
          completed.map(j =>
            `<option value="${escapeHtml(j.job_id)}">${escapeHtml(j.job_id.slice(0, 8))} (${j.pose_count} poses${j.best_affinity != null ? `, ${Number(j.best_affinity).toFixed(1)} kcal/mol` : ''})</option>`
          ).join('')
        if (current) sel.value = current
      }
    } catch {
      // No jobs yet - selection stays empty.
    }
  }

  async function runComparison() {
    const a = el('compare-job-a')?.value
    const b = el('compare-job-b')?.value
    if (!a || !b) {
      showStatus('Select two completed docking jobs', true)
      return
    }
    try {
      const data = await sendCommand('compare_results', { job_id_a: a, job_id_b: b })
      renderComparison(data)
      showStatus(`Compared: avg RMSD ${data.avg_rmsd ?? '—'} Å`)
    } catch (e) {
      showStatus(`Comparison failed: ${e.message}`, true)
    }
  }

  function renderComparison(data) {
    const pane = el('compare-results')
    if (!pane) return
    const best = data.best_matches || []
    if (best.length === 0) {
      pane.innerHTML = '<div class="v2-placeholder">No comparable poses</div>'
      return
    }
    let html = `<div class="v2-summary">average RMSD ${Number(data.avg_rmsd).toFixed(2)} Å (Kabsch-superposed, atom-name matched)</div>`
    html += '<table><thead><tr><th>A pose</th><th>B pose</th><th>RMSD Å</th><th>A aff.</th><th>B aff.</th><th>Δ kcal/mol</th></tr></thead><tbody>'
    for (const m of best) {
      html += `<tr>
        <td>${m.pose1}</td><td>${m.pose2}</td>
        <td>${Number(m.rmsd).toFixed(2)}</td>
        <td>${Number(m.pose1_affinity).toFixed(2)}</td>
        <td>${Number(m.pose2_affinity).toFixed(2)}</td>
        <td>${Number(m.energy_diff).toFixed(2)}</td>
      </tr>`
    }
    html += '</tbody></table>'
    pane.innerHTML = html
  }

  // ==================================================================
  // V2: tab switching
  // ==================================================================

  function initV2Tabs() {
    const dock = el('v2-dock')
    if (!dock) return
    dock.querySelectorAll('.v2-tab').forEach(tab => {
      tab.addEventListener('click', () => {
        dock.querySelectorAll('.v2-tab').forEach(t => t.classList.remove('active'))
        dock.querySelectorAll('.v2-panel').forEach(p => p.classList.add('v2-hidden'))
        tab.classList.add('active')
        el(`v2-panel-${tab.dataset.tab}`)?.classList.remove('v2-hidden')
        if (tab.dataset.tab === 'compare') refreshComparisonJobs()
        if (tab.dataset.tab === 'editor') drawEditorCanvas()
      })
    })
    dock.querySelector('.v2-tab')?.classList.add('active')
    el('v2-panel-batch')?.classList.remove('v2-hidden')
  }

  // ------------------------------------------------------------------
  // Wiring
  // ------------------------------------------------------------------

  function upgrade(container) {
    // Inject the ligand selector above the ligand card.
    const card = el('ligand-card')
    if (card && !el('ligand-select')) {
      const select = document.createElement('select')
      select.id = 'ligand-select'
      select.className = 'select'
      card.parentElement.insertBefore(select, card)
      select.addEventListener('change', (e) => selectLigand(e.target.value))
    }
    // Hidden until measure mode is on.
    container.addEventListener('keydown', (e) => {
      if (e.key === 'Escape' && state.measureMode) toggleMeasure()
    })
    // 2D editor canvas: atom picking.
    el('editor-canvas')?.addEventListener('click', editorPick)
  }

  async function start() {
    initViewer()

    el('btn-open-file')?.addEventListener('click', openLocalFile)
    el('btn-open-pdb')?.addEventListener('click', openPdbId)
    el('btn-run-analysis')?.addEventListener('click', runAnalysis)
    el('btn-export')?.addEventListener('click', exportArtifacts)
    el('btn-measure')?.addEventListener('click', toggleMeasure)
    el('btn-get-status')?.addEventListener('click', refreshStatus)

    el('select-representation')?.addEventListener('change', (e) => {
      state.representation = e.target.value
      applyView()
    })
    el('select-colorscheme')?.addEventListener('change', (e) => {
      state.colorScheme = e.target.value
      applyView()
    })

    let notesTimer = null
    el('notes-area')?.addEventListener('input', (e) => {
      clearTimeout(notesTimer)
      notesTimer = setTimeout(() => {
        sendCommand('set_notes', { notes: e.target.value })
          .catch(err => showStatus(`Notes not saved: ${err.message}`, true))
      }, 600)
    })

    // Job progress events from the backend (forwarded by the Rust bridge).
    // Tauri-only API: browser (dev/automation) mode skips event listening.
    const IN_TAURI = typeof window !== 'undefined' && '__TAURI_INTERNALS__' in window
    if (IN_TAURI) {
      listen('backend-event', (event) => {
        const payload = event.payload?.data || event.payload || {}
        if (payload.job_id) {
          sendCommand('get_job', { job_id: payload.job_id })
            .then(showJobResult)
            .catch(() => {})
        }
      })
    }

    // V2 surfaces.
    initV2Tabs()
    el('btn-batch-run')?.addEventListener('click', runBatch)
    el('btn-batch-export')?.addEventListener('click', exportBatch)
    el('btn-water-analyze')?.addEventListener('click', analyzeWaterNetwork)
    el('btn-editor-load')?.addEventListener('click', editorLoad)
    el('btn-editor-export')?.addEventListener('click', editorExportSdf)
    el('btn-editor-add-atom')?.addEventListener('click', editorAddAtom)
    el('btn-editor-del-atom')?.addEventListener('click', editorRemoveAtom)
    el('btn-editor-add-bond')?.addEventListener('click', () => editorBeginOp('add_bond'))
    el('btn-editor-del-bond')?.addEventListener('click', () => editorBeginOp('remove_bond'))
    el('btn-script-run')?.addEventListener('click', runScript)
    el('btn-compare-run')?.addEventListener('click', runComparison)

    await refreshStatus()
  }

  return { upgrade, start }
}
