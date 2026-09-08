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
    contactFilter: 'all',
    contactSort: 'type',
    pocketSurfaceOn: false,
    recentStructures: [],
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
    // The backend tags each file with the format name 3Dmol.js parses it
    // as (pdb, cif, mol2, sdf, xyz, gro, pdbqt, ...).
    viewer.addModel(file.content, file.format || 'pdb')
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
        // Molecular surface over a thin cartoon backbone so the fold is
        // still readable through the translucent surface.
        return [
          { sel, style: { cartoon: { color: 'spectrum', opacity: 0.6, ...clickable } } },
        ]
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
    // Ligand-centric pocket surface (spec V1 §4.4): on demand, over the
    // residues the backend's own geometric binding-pocket computation
    // selected — never a fixed radius guessed here.
    if (state.pocketSurfaceOn && state.pocketResidues?.length) {
      const type = $3Dmol.SurfaceType ? $3Dmol.SurfaceType.VDW : 'VDW'
      for (const r of state.pocketResidues) {
        viewer.addSurface(type, { opacity: 0.55, color: '#a371f7' },
          { chain: r.chain_id, resi: r.residue_id })
      }
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
  // Camera presets (spec V1 §5.2)
  // ------------------------------------------------------------------

  function applyCameraPreset(name) {
    if (!viewer) return
    if (name === 'all') {
      viewer.zoomTo()
    } else if (name === 'ligand' && state.selectedLigand) {
      viewer.zoomTo({ resn: state.selectedLigand.residue_name })
    } else if (name === 'pocket' && state.pocketResidues?.length) {
      // Zoom to the union of pocket residues (3Dmol accepts a sel dict;
      // for multiple residues use the first as anchor and fit all via OR).
      const first = state.pocketResidues[0]
      viewer.zoomTo({ chain: first.chain_id, resi: first.residue_id })
    } else if (name === 'interface' && state.contacts.length) {
      viewer.zoomTo({
        chain: state.contacts[0].protein_chain_id,
        resi: state.contacts[0].protein_residue_id,
      })
    } else {
      viewer.zoomTo()
    }
    viewer.render()
  }

  // ------------------------------------------------------------------
  // Scene image export (spec V1 §4.6): real rendered viewport pixels.
  // ------------------------------------------------------------------

  async function exportSceneImage() {
    if (!viewer || !state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    try {
      // 3Dmol.js renders into a canvas inside the scene container.
      const canvas = viewerElement?.querySelector('canvas')
      if (!canvas) throw new Error('No rendered canvas found')
      // preserveDrawingBuffer is off by default; render() right before
      // toDataURL keeps the buffer valid.
      viewer.render()
      const dataUrl = canvas.toDataURL('image/png')
      const b64 = dataUrl.split(',')[1]
      if (!b64) throw new Error('Canvas capture produced no image data')
      const result = await sendCommand('save_scene_image', { png_base64: b64 })
      const dest = await save({
        defaultPath: `scene_${state.structure.id}.png`,
        filters: [{ name: 'PNG', extensions: ['png'] }],
      })
      if (dest) {
        await invoke('copy_file', { src: result.path, dest })
        showStatus(`Scene image saved (${result.size_bytes} bytes)`)
      }
    } catch (e) {
      showStatus(`Scene image export failed: ${e.message}`, true)
    }
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
      </div>
      <div class="v2-actions" style="margin-top:6px">
        <button class="btn btn-secondary btn-sm" id="btn-ligand-similar">Find similar components (RCSB)</button>
      </div>`
    el('btn-ligand-similar')?.addEventListener('click', () => {
      switchTab('discover')
      el('discover-mode').value = 'similar_components'
      runDiscovery()
    })
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
    state.contacts = contacts || []
    if (!contacts || contacts.length === 0) {
      pane.innerHTML = `<div class="contacts-placeholder">${
        plipAvailable === false
          ? 'PLIP is not installed - contact analysis unavailable'
          : 'No contacts found'}</div>`
      return
    }

    // Spec V1 §5.3: sort/filter basics over the real contact rows.
    const types = [...new Set(contacts.map(c => c.contact_type || 'other'))].sort()
    const controls = `<div class="v2-actions" style="margin-bottom:6px">
      <select class="select" id="contact-filter" style="max-width:130px">
        <option value="all"${state.contactFilter === 'all' ? ' selected' : ''}>All types</option>
        ${types.map(t => `<option value="${escapeHtml(t)}"${state.contactFilter === t ? ' selected' : ''}>${escapeHtml(t)}</option>`).join('')}
      </select>
      <select class="select" id="contact-sort" style="max-width:120px">
        <option value="type"${state.contactSort === 'type' ? ' selected' : ''}>Sort: type</option>
        <option value="distance"${state.contactSort === 'distance' ? ' selected' : ''}>Sort: distance</option>
        <option value="residue"${state.contactSort === 'residue' ? ' selected' : ''}>Sort: residue</option>
      </select>
    </div>`

    let rows = contacts.filter(c =>
      state.contactFilter === 'all' || (c.contact_type || 'other') === state.contactFilter)
    rows = [...rows].sort((a, b) => {
      if (state.contactSort === 'distance') {
        return (a.distance ?? Infinity) - (b.distance ?? Infinity)
      }
      if (state.contactSort === 'residue') {
        const ka = `${a.protein_chain_id}${a.protein_residue_id}`
        const kb = `${b.protein_chain_id}${b.protein_residue_id}`
        return ka.localeCompare(kb)
      }
      const ta = a.contact_type || 'other'
      const tb = b.contact_type || 'other'
      return ta === tb
        ? (a.distance ?? Infinity) - (b.distance ?? Infinity)
        : ta.localeCompare(tb)
    })

    const diagramBtn = `<div class="v2-actions" style="margin-bottom:6px">
      <button class="btn btn-secondary btn-sm" id="btn-contacts-diagram">2D diagram</button>
      <button class="btn btn-secondary btn-sm" id="btn-contacts-structconn">Export _struct_conn (mmCIF)</button>
    </div>`
    let html = '<table><thead><tr><th>Type</th><th>Protein</th><th>Lig atom</th><th>Å</th></tr></thead><tbody>'
    for (const c of rows.slice(0, 100)) {
      html += `<tr class="contact-row" data-chain="${escapeHtml(c.protein_chain_id)}" data-resi="${escapeHtml(String(c.protein_residue_id))}">
        <td><span class="contact-badge">${escapeHtml(c.contact_type || '')}</span></td>
        <td>${escapeHtml(c.protein_residue_name)}${escapeHtml(String(c.protein_residue_id))} (${escapeHtml(c.protein_chain_id)})</td>
        <td>${escapeHtml(c.ligand_atom || '')}</td>
        <td>${c.distance != null ? Number(c.distance).toFixed(2) : '—'}</td>
      </tr>`
    }
    html += '</tbody></table>'
    if (rows.length > 100) {
      html += `<div class="contacts-more">... and ${rows.length - 100} more</div>`
    }
    pane.innerHTML = diagramBtn + controls + html
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
    el('btn-contacts-diagram')?.addEventListener('click', () => {
      switchTab('structure')
      showInteractionDiagram()
    })
    el('btn-contacts-structconn')?.addEventListener('click', exportStructConn)
    el('contact-filter')?.addEventListener('change', (e) => {
      state.contactFilter = e.target.value
      updateContactsPanel(state.contacts, true)
    })
    el('contact-sort')?.addEventListener('change', (e) => {
      state.contactSort = e.target.value
      updateContactsPanel(state.contacts, true)
    })
  }

  async function exportContactsCsv() {
    try {
      const result = await sendCommand('export_contacts_csv', {})
      const dest = await save({
        defaultPath: `contacts_${state.structure?.id || 'structure'}.csv`,
        filters: [{ name: 'CSV', extensions: ['csv'] }],
      })
      if (dest) {
        await invoke('copy_file', { src: result.path, dest })
        showStatus('Contacts CSV saved')
      }
    } catch (e) {
      showStatus(`Contacts CSV export failed: ${e.message}`, true)
    }
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
      filters: [{ name: 'Structure files', extensions: ['cif', 'mmcif', 'cif.gz', 'pdb', 'ent', 'pdb.gz', 'mol2', 'sdf', 'mol', 'xyz', 'gro', 'pdbqt'] }],
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
    state.pocketResidues = null
    state.pocketSurfaceOn = false
    state.selectedLigand = null
    state.selectedLigandId = null
    state.picks = []
    rememberStructure(data.structure?.id)
    updateLigandSelect()
    updateContactsPanel([], true)
    updateEvidencePanel([])
    await loadStructureIntoViewer()
    // Auto-select the first ligand that is neither solvent (HETAS) nor an
    // ion (HETAI) per its CCD classification. When the structure has none
    // (protein-only, or ion/solvent only), select nothing and say so -
    // never present a water molecule as the ligand.
    // Viewer-only formats (mol2/sdf/xyz/gro/...): rendered in 3D, but the
    // analysis model was never parsed — say so instead of pretending.
    if (state.structure.viewer_only) {
      updateLigandCard(null)
      showStatus(`Loaded ${state.structure.id} (${state.structure.file_format}) — view only: ligand/contacts analysis needs PDB or mmCIF`, true)
      return
    }
    const first = (state.structure.ligands || []).find(
      l => l.classification_hint &&
        !['HETAS', 'HETAI'].includes(l.classification_hint))
    if (first) {
      await selectLigand(first.id)
    } else {
      updateLigandCard(null)
      showStatus(`Loaded ${state.structure.id} (no drug-like ligand: only solvent/ions present)`)
      return
    }
    showStatus(`Loaded ${state.structure.id}`)
  }

  async function selectLigand(ligandId) {
    state.selectedLigandId = ligandId
    try {
      const data = await sendCommand('select_ligand', { ligand_id: ligandId })
      state.selectedLigand = data.ligand
      rememberLigand(ligandId)
      updateLigandCard(data.ligand)
      updateLigandSelect()
      applyView()
      if (data.warning) showStatus(data.warning, true)
    } catch (e) {
      showStatus(`Ligand selection failed: ${e.message}`, true)
    }
  }

  async function runAnalysis() {
    setBusy(true, 'Running PLIP contact analysis...')
    try {
      const data = await sendCommand('run_contact_analysis', {})
      state.contacts = data.contacts || []
      // The binding pocket from the backend's own geometric computation
      // feeds the pocket-surface toggle (no radius guessed in the UI).
      state.pocketResidues = (data.pocket?.pocket_residues || []).map(r => ({
        chain_id: r.chain_id,
        residue_id: r.residue_id,
      }))
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

  // Spec V2 §4.8: 2D edit → 3D sync. The edited molecule is overlaid as
  // its own model (green carbons) next to the original ligand — nothing
  // in the loaded structure is modified.
  async function editorSyncTo3d() {
    if (!viewer) return
    try {
      const result = await sendCommand('editor_sync_to_3d', {})
      // Remove any previous synced model, keep everything else.
      const models = viewer.models || []
      if (models.length > 1) {
        viewer.removeAllModels()
        await loadStructureIntoViewer()
      }
      viewer.addModel(result.content, 'sdf')
      const modelId = (viewer.models?.length || 2) - 1
      viewer.setStyle({ model: modelId }, {
        stick: { radius: 0.28, colorscheme: 'greenCarbon' },
        sphere: { scale: 0.3, colorscheme: 'greenCarbon' },
      })
      viewer.zoomTo({ model: modelId })
      viewer.render()
      showStatus(`Edited molecule in 3D: ${result.atom_count} atoms, ${result.bond_count} bonds (green)`)
    } catch (e) {
      showStatus(`Sync to 3D failed: ${e.message}`, true)
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
    const comparison = data.comparison || {}
    if (comparison.error) {
      pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(comparison.error)}</div>`
      return
    }
    const comparisons = comparison.comparisons || []
    let html = `<div class="v2-summary">average RMSD ${Number(data.avg_rmsd ?? 0).toFixed(2)} Å across ${comparison.num_poses1 ?? 0}×${comparison.num_poses2 ?? 0} pose pairs</div>`
    if (comparisons.length > 0) {
      html += '<table><thead><tr><th>Run A pose</th><th>Run B pose</th><th>RMSD (Å)</th><th>Affinity A</th><th>Affinity B</th><th>ΔE</th></tr></thead><tbody>'
      for (const c of comparisons) {
        html += `<tr><td>#${c.pose1}</td><td>#${c.pose2}</td><td>${Number(c.rmsd).toFixed(2)}</td><td>${Number(c.pose1_affinity).toFixed(2)}</td><td>${Number(c.pose2_affinity).toFixed(2)}</td><td>${Number(c.energy_diff).toFixed(2)}</td></tr>`
      }
      html += '</tbody></table>'
    } else {
      html += '<div class="v2-placeholder">No pose pairs to compare</div>'
    }
    pane.innerHTML = html
  }

  // Spec V2 §4.7: pose clustering + consensus over both jobs' poses.
  async function runClustering() {
    const a = el('compare-job-a')?.value
    const b = el('compare-job-b')?.value
    if (!a || !b) {
      showStatus('Select two completed docking jobs', true)
      return
    }
    const threshold = parseFloat(el('cluster-threshold')?.value || '2.0')
    try {
      const data = await sendCommand('compare_results', {
        job_id_a: a, job_id_b: b, cluster_rmsd_threshold: threshold,
      })
      renderClustering(data)
    } catch (e) {
      showStatus(`Clustering failed: ${e.message}`, true)
    }
  }

  function renderClustering(data) {
    const pane = el('compare-results')
    if (!pane) return
    const clusters = data.clusters || []
    const consensus = data.consensus_pose
    let html = ''
    if (data.clustering_note) {
      html += `<div class="v2-error v2-summary">${escapeHtml(data.clustering_note)}</div>`
    }
    html += `<div class="v2-summary">${clusters.length} cluster(s) at ≤${data.cluster_rmsd_threshold ?? 2.0} Å pairwise RMSD</div>`
    if (clusters.length > 0) {
      html += '<table><thead><tr><th>Size</th><th>Poses (job:index)</th><th>Representative</th></tr></thead><tbody>'
      for (const c of clusters) {
        html += `<tr><td>${c.size}</td><td>${c.result_indices.map(i => `#${i}`).join(', ')}</td><td>#${c.representative}</td></tr>`
      }
      html += '</tbody></table>'
    } else {
      html += '<div class="v2-placeholder">No multi-pose clusters at this threshold (each pose is its own cluster)</div>'
    }
    if (consensus) {
      html += `<div class="v2-summary">Consensus pose (lowest affinity across both runs): ${Number(consensus.affinity).toFixed(2)} kcal/mol${consensus.ligand_name ? ` · ${escapeHtml(consensus.ligand_name)}` : ''}</div>`
      html += `<div class="v2-actions"><button class="btn btn-secondary btn-sm" id="btn-consensus-view">Show consensus in 3D</button></div>`
    }
    pane.innerHTML = html
    el('btn-consensus-view')?.addEventListener('click', () => {
      if (!viewer || !consensus?.atoms?.length) return
      viewer.removeAllModels()
      loadStructureIntoViewer()
      const lines = ['MODEL        1']
      let serial = 1
      for (const a of consensus.atoms) {
        const elem = (a.element || 'C').toUpperCase().slice(0, 2).padEnd(2)
        lines.push(
          `HETATM${String(serial).padStart(5, ' ')}  ${(a.name || elem.trim()).padEnd(4).slice(0, 4)}` +
          ` LIG L   1    ` +
          `${a.x.toFixed(3).padStart(8)}${a.y.toFixed(3).padStart(8)}${a.z.toFixed(3).padStart(8)}` +
          `  1.00  0.00          ${elem}`)
        serial++
      }
      lines.push('ENDMDL', 'END')
      viewer.addModel(lines.join('\n'), 'pdb')
      viewer.setStyle({ model: 1 }, {
        stick: { radius: 0.25, colorscheme: 'greenCarbon' },
        sphere: { scale: 0.3, colorscheme: 'greenCarbon' },
      })
      viewer.zoomTo({ model: 1 })
      viewer.render()
    })
  }

  // ==================================================================
  // V2: docking tab (real Vina run + poses rendered in the 3D viewer)
  // ==================================================================

  let lastDockJobId = null
  let poseVisible = false

  async function runDocking() {
    if (!state.structure || !state.selectedLigandId) {
      showStatus('Open a structure and select a ligand first', true)
      return
    }
    const exhaustiveness = Math.max(1, parseInt(el('dock-exhaustiveness')?.value || '8', 10))
    const numModes = Math.max(1, parseInt(el('dock-num-modes')?.value || '5', 10))
    const progress = el('dock-progress')
    if (progress) progress.style.display = 'flex'
    const results = el('dock-results')
    if (results) results.innerHTML = '<div class="v2-placeholder">Running real AutoDock Vina (Open Babel prep + docking)…</div>'
    try {
      const started = await sendCommand('run_docking', {
        exhaustiveness, num_modes: numModes,
      })
      lastDockJobId = started.job_id
      let job = null
      for (let i = 0; i < 240; i++) {
        job = await sendCommand('get_job', { job_id: started.job_id })
        if (job.status === 'completed' || job.status === 'failed') break
        await new Promise(r => setTimeout(r, 2000))
      }
      if (!job || job.status !== 'completed') {
        throw new Error(job?.error || 'docking did not complete')
      }
      state.lastDockJob = job
      renderDocking(job)
      showStatus(`Docking done: ${job.result.poses.length} poses`)
    } catch (e) {
      if (results) results.innerHTML = `<div class="v2-error v2-summary">Docking failed: ${escapeHtml(e.message)}</div>`
      showStatus(`Docking failed: ${e.message}`, true)
    } finally {
      if (progress) progress.style.display = 'none'
    }
  }

  function renderDocking(job) {
    const results = el('dock-results')
    if (!results) return
    const poses = job.result?.poses || []
    if (poses.length === 0) {
      results.innerHTML = '<div class="v2-placeholder">No poses returned</div>'
      return
    }
    let html = `<div class="v2-summary">${poses.length} poses from Vina · box from ligand extent · job ${escapeHtml(job.job_id.slice(0, 8))}</div>`
    html += '<table><thead><tr><th>#</th><th>Affinity kcal/mol</th><th></th></tr></thead><tbody>'
    for (const p of poses) {
      html += `<tr>
        <td>${p.pose_id}</td>
        <td>${Number(p.affinity).toFixed(2)}</td>
        <td><button class="btn btn-secondary btn-sm dock-pose-btn" data-pose="${p.pose_id}">Show in 3D</button></td>
      </tr>`
    }
    html += '</tbody></table>'
    results.innerHTML = html
    results.querySelectorAll('.dock-pose-btn').forEach(btn => {
      btn.addEventListener('click', () => showPoseInViewer(
        parseInt(btn.dataset.pose, 10)))
    })
  }

  /** Render a real Vina pose as a model in the 3D viewer. */
  function showPoseInViewer(poseId) {
    if (!viewer) return
    const job = state.lastDockJob
    if (!job) return
    const pose = (job.result?.poses || []).find(p => p.pose_id === poseId)
    if (!pose) return
    viewer.removeAllModels()
    // Re-add the structure, then the pose atoms as a separate model.
    loadStructureIntoViewer()
    const pdbLines = ['MODEL        1']
    let serial = 1
    for (const a of pose.atoms) {
      const elem = (a.element || 'C').toUpperCase().slice(0, 2).rjust(2)
      pdbLines.push(
        `HETATM${String(serial).padStart(5, ' ')}  ${(a.name || elem.trim()).padEnd(4).slice(0, 4)}` +
        ` LIG L   1    ` +
        `${a.x.toFixed(3).padStart(8)}${a.y.toFixed(3).padStart(8)}${a.z.toFixed(3).padStart(8)}` +
        `  1.00  0.00          ${elem}`)
      serial++
    }
    pdbLines.push('ENDMDL', 'END')
    viewer.addModel(pdbLines.join('\n'), 'pdb')
    viewer.setStyle({ model: 1 }, {
      stick: { radius: 0.25, colorscheme: 'orangeCarbon' },
      sphere: { scale: 0.3, colorscheme: 'orangeCarbon' },
    })
    viewer.zoomTo({ model: 1 })
    viewer.render()
    poseVisible = true
    showStatus(`Pose ${poseId}: ${Number(pose.affinity).toFixed(2)} kcal/mol`)
  }

  // ==================================================================
  // V2: BM25 search tab
  // ==================================================================

  async function runSearch() {
    const query = el('search-input')?.value?.trim()
    if (!query) {
      showStatus('Type a query first', true)
      return
    }
    const scope = el('search-scope')?.value || ''
    try {
      const data = await sendCommand('search', {
        query, scope: scope || undefined, limit: 25,
      })
      renderSearch(data)
    } catch (e) {
      showStatus(`Search failed: ${e.message}`, true)
    }
  }

  function renderSearch(data) {
    const pane = el('search-results')
    if (!pane) return
    const status = el('search-status')
    if (status) status.textContent = `${data.documents_indexed || 0} docs indexed`
    if (data.needs_refresh) {
      pane.innerHTML = '<div class="v2-placeholder">No indexed documents yet — click “Index current session” to fetch from live sources.</div>'
      return
    }
    if (!data.hits || data.hits.length === 0) {
      pane.innerHTML = `<div class="v2-placeholder">No hits for “${escapeHtml(data.query)}” in ${data.documents_indexed} indexed documents</div>`
      return
    }
    let html = `<div class="v2-summary">${data.total} hits for “${escapeHtml(data.query)}” (BM25 over ${data.documents_indexed} live-fetched documents)</div>`
    html += '<table><thead><tr><th>ID</th><th>Description</th><th>Scope</th><th>Score</th><th>Source</th></tr></thead><tbody>'
    for (const h of data.hits) {
      const f = h.fields || {}
      const desc = f.title || f.name || f.text || ''
      html += `<tr class="search-hit" data-doc="${escapeHtml(h.doc_id)}">
        <td>${escapeHtml((f.id || h.doc_id).slice(0, 16))}</td>
        <td>${escapeHtml(String(desc).slice(0, 60))}</td>
        <td>${escapeHtml(h.scope)}</td>
        <td>${Number(h.score).toFixed(3)}</td>
        <td>${escapeHtml(h.source)}</td>
      </tr>`
    }
    html += '</tbody></table>'
    pane.innerHTML = html
    pane.querySelectorAll('.search-hit').forEach(row => {
      row.addEventListener('click', () => {
        const docId = row.dataset.doc || ''
        // A structure hit opens that PDB entry through the real fetch path;
        // a chemical hit loads the component's real enrichment card.
        if (docId.startsWith('struct:')) {
          openPdbIdValue(docId.slice(7))
        } else if (docId.startsWith('chem:')) {
          enrichComponent(docId.slice(5))
        }
      })
    })
  }

  async function openPdbIdValue(pdbId) {
    if (!/^[0-9][A-Za-z0-9]{3}$/.test(pdbId)) return
    setBusy(true, `Fetching ${pdbId} from RCSB...`)
    try {
      const data = await sendCommand('open_pdb_id', { pdb_id: pdbId })
      await onStructureLoaded(data)
    } catch (e) {
      showStatus(`Failed to open ${pdbId}: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  // ------------------------------------------------------------------
  // Top bar: PDB-wide search, recent structures, settings (spec V1 §5.1)
  // ------------------------------------------------------------------

  async function topbarSearch() {
    const query = el('topbar-search')?.value?.trim()
    if (!query) {
      showStatus('Type a search term first', true)
      return
    }
    setBusy(true, 'Searching the PDB…')
    try {
      const data = await sendCommand('discover_text', { text: query, rows: 15 })
      switchTab('discover')
      el('discover-input').value = query
      renderDiscovery(data)
      showStatus(`${data.total_count} PDB hits for "${query}"`)
    } catch (e) {
      showStatus(`Search failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  function rememberStructure(pdbId) {
    if (!pdbId) return
    state.recentStructures = [
      pdbId, ...state.recentStructures.filter(x => x !== pdbId),
    ].slice(0, 8)
    const sel = el('select-recent')
    if (sel) {
      sel.innerHTML = '<option value="">Recent…</option>' +
        state.recentStructures.map(id =>
          `<option value="${escapeHtml(id)}">${escapeHtml(id)}</option>`).join('')
    }
  }

  function showSettings() {
    let modal = el('app-modal')
    if (!modal) {
      modal = document.createElement('div')
      modal.id = 'app-modal'
      modal.className = 'modal-backdrop'
      document.body.appendChild(modal)
    }
    modal.innerHTML = `
      <div class="modal">
        <div class="modal-title">Engine & data-source status (live)</div>
        <div id="settings-body" class="settings-body">Checking…</div>
        <div class="modal-actions">
          <button class="btn btn-primary" id="modal-ok">Close</button>
        </div>
      </div>`
    modal.style.display = 'flex'
    el('modal-ok').addEventListener('click', () => {
      modal.style.display = 'none'
      modal.innerHTML = ''
    })
    sendCommand('get_status', {}).then(st => {
      const engines = Object.entries(st.engines || {})
        .map(([name, info]) =>
          `<div class="evidence-item"><div class="evidence-source">${escapeHtml(name)}</div>` +
          `<div class="evidence-value">${info.available ? 'available' : 'NOT installed'}</div></div>`)
        .join('')
      const sources = Object.entries(st.data_sources || {})
        .map(([name, v]) =>
          `<div class="evidence-item"><div class="evidence-source">${escapeHtml(name)}</div>` +
          `<div class="evidence-value">${escapeHtml(typeof v === 'string' ? v : JSON.stringify(v))}</div></div>`)
        .join('')
      const body = el('settings-body')
      if (body) body.innerHTML =
        `<h3>Engines</h3>${engines}<h3>Data sources</h3>${sources}`
    }).catch(e => {
      const body = el('settings-body')
      if (body) body.textContent = `Status check failed: ${e.message}`
    })
  }

  async function refreshSearchIndex() {
    const status = el('search-status')
    if (status) status.textContent = 'fetching from live sources…'
    try {
      const data = await sendCommand('search_index_refresh', {})
      if (status) {
        status.textContent = `${data.total} docs indexed (${data.indexed} new)`
      }
      showStatus(`Search index refreshed: ${data.total} documents`)
    } catch (e) {
      if (status) status.textContent = ''
      showStatus(`Index refresh failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2.1: 2D interaction diagram (real PLIP contacts -> SVG)
  // ==================================================================

  async function showInteractionDiagram() {
    if (!state.structure) {
      showStatus('Open a structure and run analysis first', true)
      return
    }
    setBusy(true, 'Building 2D interaction diagram...')
    try {
      const data = await sendCommand('get_interaction_diagram', {
        width: 900, height: 640,
      })
      const pane = el('structure-results')
      if (pane) {
        pane.innerHTML = `
          <div class="v2-summary">Interaction diagram — ${escapeHtml(data.ligand_name || '')}: ${data.residues.length} contacting residues · layout: ${escapeHtml(data.layout_source || '')}</div>
          <div class="diagram-wrap">${data.svg}</div>`
      }
      showStatus('2D interaction diagram built')
    } catch (e) {
      showStatus(`Diagram failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  async function exportStructConn() {
    try {
      const result = await sendCommand('export_struct_conn', {})
      const dest = await save({
        defaultPath: `contacts_${state.structure?.id || 'structure'}.cif`,
        filters: [{ name: 'mmCIF', extensions: ['cif'] }],
      })
      if (dest) {
        await invoke('copy_file', { src: result.path, dest })
        showStatus(`_struct_conn exported (${result.contact_count} contacts)`)
      }
    } catch (e) {
      showStatus(`struct_conn export failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2.1: sequence view with interaction/SS/pocket mapping
  // ==================================================================

  const SS_COLORS = {
    'alpha-helix': '#f28e2b', 'beta-strand': '#4e79ba', coil: '#8b949e',
    'beta-bridge': '#76b7b2', 'bend': '#b07aa1', 'turn': '#edc948',
    'polyproline': '#59a14f',
  }

  async function showSequenceView() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    try {
      const data = await sendCommand('get_sequence_view', {})
      const pane = el('structure-results')
      if (!pane) return
      let html = `<div class="v2-summary">${data.chains.length} chain(s) · contacts: ${data.has_contacts ? 'yes' : 'run analysis first'} · DSSP: ${data.has_dssp ? 'yes' : 'not run'} · pockets: ${data.has_pockets ? 'yes' : 'not run'}</div>`
      for (const chain of data.chains) {
        html += `<div class="v2-section-title">Chain ${escapeHtml(chain.chain_id)} — ${chain.length} residues</div>`
        // Row 1: secondary structure strip (when DSSP has run).
        // Row 2: the sequence itself; contacting residues highlighted.
        const residues = [...chain.residues].sort((a, b) => a.index - b.index)
        if (data.has_dssp && residues.some(r => r.ss && r.ss !== 'coil')) {
          html += '<div class="seq-row seq-ss">'
          for (const r of residues) {
            const color = SS_COLORS[r.ss] || '#8b949e'
            html += `<span class="seq-cell" style="background:${color}" title="${escapeHtml(r.ss)} ${r.id}"></span>`
          }
          html += '</div>'
        }
        html += '<div class="seq-row seq-aa">'
        for (const r of residues) {
          const hasContact = (r.contacts || []).length > 0
          const pocket = r.in_pocket
          const cls = hasContact ? 'seq-contact' : (pocket ? 'seq-pocket' : '')
          const tip = `${r.name}${r.id}${hasContact ? ' — ' + r.contacts.map(c => c.type).join(', ') : ''}${pocket ? ' — in pocket' : ''}`
          html += `<span class="seq-cell seq-aa-cell ${cls}" title="${escapeHtml(tip)}" data-chain="${escapeHtml(chain.chain_id)}" data-resi="${r.id}">${escapeHtml(r.name[0] || '?')}</span>`
        }
        html += '</div>'
      }
      pane.innerHTML = html
      pane.querySelectorAll('.seq-aa-cell.seq-contact').forEach(cell => {
        cell.addEventListener('click', () => {
          if (!viewer) return
          viewer.zoomTo({ chain: cell.dataset.chain, resi: parseInt(cell.dataset.resi, 10) })
          viewer.render()
        })
      })
      showStatus('Sequence view rendered')
    } catch (e) {
      showStatus(`Sequence view failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2.1: structure quality + protein context (live sources)
  // ==================================================================

  async function showEntryQuality() {
    try {
      const data = await sendCommand('entry_quality', {})
      const pane = el('structure-results')
      if (!pane) return
      const q = data.quality || {}
      const rows = Object.entries(q).map(([k, v]) =>
        `<div class="ligand-property"><div class="ligand-property-label">${escapeHtml(k)}</div><div class="ligand-property-value">${escapeHtml(String(v))}</div></div>`).join('')
      pane.innerHTML = `
        <div class="v2-summary">Quality for ${escapeHtml(data.pdb_id)} — ${data.has_validation_report ? 'validation report present' : 'no validation report for this entry'}</div>
        ${rows || '<div class="v2-placeholder">No quality fields reported for this entry</div>'}
        ${data.url ? `<a class="evidence-url" href="${escapeHtml(data.url)}" target="_blank" rel="noopener">${escapeHtml(data.url)}</a>` : ''}`
      showStatus('Quality metrics fetched')
    } catch (e) {
      showStatus(`Quality fetch failed: ${e.message}`, true)
    }
  }

  async function showUniprotContext() {
    try {
      const data = await sendCommand('uniprot_context', {})
      const pane = el('structure-results')
      if (!pane) return
      const fmtSites = (sites) => sites.map(s =>
        `${s.start}${s.end && s.end !== s.start ? '-' + s.end : ''} ${escapeHtml(s.description || '')}${s.ligand ? ' (' + escapeHtml(s.ligand) + ')' : ''}`).join('; ') || '—'
      pane.innerHTML = `
        <div class="v2-summary">${escapeHtml(data.entry_name || data.accession)} · ${escapeHtml(data.organism || '')} · ${data.sequence_length || '?'} aa</div>
        <div class="ligand-property"><div class="ligand-property-label">Function</div><div class="ligand-property-value">${escapeHtml((data.function || []).join(' ') || '—')}</div></div>
        <div class="ligand-property"><div class="ligand-property-label">Binding sites</div><div class="ligand-property-value">${fmtSites(data.binding_sites || [])}</div></div>
        <div class="ligand-property"><div class="ligand-property-label">Active sites</div><div class="ligand-property-value">${fmtSites(data.active_sites || [])}</div></div>
        <div class="ligand-property"><div class="ligand-property-label">Domains</div><div class="ligand-property-value">${fmtSites(data.domains || [])}</div></div>
        <a class="evidence-url" href="${escapeHtml(data.url)}" target="_blank" rel="noopener">${escapeHtml(data.url)}</a>`
      showStatus('UniProt context fetched')
    } catch (e) {
      showStatus(`UniProt fetch failed: ${e.message}`, true)
    }
  }

  async function showAlphafoldModel() {
    try {
      const data = await sendCommand('alphafold_model', {})
      const pane = el('structure-results')
      if (!pane) return
      pane.innerHTML = `
        <div class="v2-summary">AlphaFold model available for ${escapeHtml(data.accession)} (${data.size_chars} chars, mmCIF)</div>
        <a class="evidence-url" href="${escapeHtml(data.url)}" target="_blank" rel="noopener">${escapeHtml(data.url)}</a>`
      showStatus('AlphaFold model metadata fetched')
    } catch (e) {
      showStatus(`AlphaFold lookup failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2.2: BindingDB affinities, covalent links, GROMACS MD, manifest
  // ==================================================================

  async function showBindingDB() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    setBusy(true, 'Querying BindingDB...')
    const pane = el('structure-results')
    try {
      const data = await sendCommand('bindingdb_affinity', {})
      if (!pane) return
      if (!data.records || data.records.length === 0) {
        pane.innerHTML = `<div class="v2-summary">BindingDB answered: ${data.record_count} binding records for ${escapeHtml(data.pdb_id || data.uniprot || '')} (${escapeHtml(data.source || '')}).</div>`
        showStatus('BindingDB: no records')
        return
      }
      let html = `<div class="v2-summary">${data.record_count} real binding record(s) from BindingDB (${escapeHtml(data.source || '')})${data.uniprot_accession ? ` · UniProt ${escapeHtml(data.uniprot_accession)}` : ''} · cutoff ${data.affinity_cutoff_nm} nM</div>`
      html += '<table><thead><tr><th>Affinity type</th><th>Value (nM)</th><th>Ligand SMILES</th><th>PMID / DOI</th></tr></thead><tbody>'
      for (const r of data.records.slice(0, 30)) {
        const mark = r.matches_selected_ligand ? ' ✔ (matches selected ligand)' : ''
        const link = r.doi
          ? `<a class="evidence-url" href="https://doi.org/${escapeHtml(r.doi)}" target="_blank" rel="noopener">${escapeHtml(r.doi)}</a>`
          : (r.pmid ? escapeHtml(r.pmid) : '—')
        html += `<tr><td>${escapeHtml(r.affinity_type || '—')}${mark}</td><td>${escapeHtml(r.affinity_nm ?? '—')}</td><td class="mono">${escapeHtml((r.smiles || '').slice(0, 70))}</td><td>${link}</td></tr>`
      }
      html += '</tbody></table>'
      if (data.records.length > 30) html += `<div class="v2-placeholder">…and ${data.records.length - 30} more records</div>`
      pane.innerHTML = html
      showStatus('BindingDB records fetched')
    } catch (e) {
      if (pane) pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(e.message)}</div>`
      showStatus(`BindingDB lookup failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  async function showCovalentLinks() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    setBusy(true, 'Reading covalent-link records...')
    const pane = el('structure-results')
    try {
      const data = await sendCommand('get_covalent_links', {})
      if (!pane) return
      if (!data.available) {
        pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(data.error || 'Covalent-link records unavailable')}</div>`
        return
      }
      const ligands = Object.keys(data.covalent_ligands || {})
      if (data.links.length === 0) {
        pane.innerHTML = `<div class="v2-summary">No covalent ligand links declared in this structure's ${escapeHtml(data.source || '')}. Nothing is inferred from geometry.</div>`
        showStatus('No covalent ligand links')
        return
      }
      let html = `<div class="v2-summary">${data.links.length} covalent link(s) declared in the file's own records (${escapeHtml(data.source || '')})</div>`
      html += '<table><thead><tr><th>Type</th><th>Partner 1</th><th>Partner 2</th><th>Distance Å</th></tr></thead><tbody>'
      for (const l of data.links) {
        const f = (p) => `${escapeHtml(p.comp || '')} ${escapeHtml(p.chain || '')}${escapeHtml(p.seq ?? '')}@${escapeHtml(p.atom || '')}`
        html += `<tr><td>${escapeHtml(l.conn_type)}</td><td>${f(l.partner1)}</td><td>${f(l.partner2)}</td><td>${l.distance != null ? Number(l.distance).toFixed(3) : '—'}</td></tr>`
      }
      html += '</tbody></table>'
      if (data.selected_ligand_is_covalent) {
        html += `<div class="v2-error v2-summary">The selected ligand is COVALENTLY ATTACHED to the polymer (depositor-stated). Non-covalent contact analysis still works, but interpret its table with that in mind.</div>`
      }
      pane.innerHTML = html
      showStatus('Covalent-link records read')
    } catch (e) {
      if (pane) pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(e.message)}</div>`
      showStatus(`Covalent-link lookup failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  async function runMd() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    const mode = el('md-mode')?.value || 'em'
    const nsteps = Math.max(100, parseInt(el('md-steps')?.value || '5000', 10))
    const results = el('md-results')
    if (results) results.innerHTML = `<div class="v2-placeholder">Running real GROMACS (${mode}) — pdb2gmx → grompp → mdrun…</div>`
    try {
      const started = await sendCommand('run_md', { mode, nsteps })
      let job = null
      for (let i = 0; i < 900; i++) {
        job = await sendCommand('get_job', { job_id: started.job_id })
        if (job.status === 'completed' || job.status === 'failed') break
        await new Promise(r => setTimeout(r, 2000))
      }
      if (!job || job.status !== 'completed') {
        throw new Error(job?.error || `${mode} did not complete`)
      }
      const res = job.result || {}
      const energies = res.energies?.series || {}
      const energyStr = Object.entries(energies).map(([k, v]) => `${k}: ${Number(v).toFixed(2)}`).join(' · ')
      if (results) {
        results.innerHTML = `<div class="v2-summary">GROMACS ${escapeHtml(res.mode || mode)} finished in ${res.runtime_seconds ?? '?'}s · ${res.atom_count} polymer atoms · ${escapeHtml(res.force_field || '')} · ${escapeHtml(res.scope_note || '')}</div>
          ${energyStr ? `<div class="v2-summary">Final energy: ${escapeHtml(energyStr)}</div>` : ''}
          <div class="v2-placeholder">Artifacts (topology, run input, final structure, energies) saved in the session workspace.</div>`
      }
      showStatus(`GROMACS ${mode} done`)
    } catch (e) {
      if (results) results.innerHTML = `<div class="v2-error v2-summary">GROMACS run failed: ${escapeHtml(e.message)}</div>`
      showStatus(`GROMACS run failed: ${e.message}`, true)
    }
  }

  async function exportReproManifest() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    setBusy(true, 'Writing reproducibility manifest...')
    try {
      const result = await sendCommand('export_repro_manifest', {})
      const pane = el('structure-results')
      if (pane) {
        pane.innerHTML = `<div class="v2-summary">Reproducibility manifest written: ${result.artifact_count} artifact(s), ${result.job_count} job(s) recorded.</div>
          <div class="v2-placeholder">${escapeHtml(result.path)}\nsha256: ${escapeHtml(result.sha256 || '')}</div>`
      }
      showStatus('Reproducibility manifest written')
    } catch (e) {
      showStatus(`Manifest export failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  // ==================================================================
  // V2.2: PDBe EU-mirror annotations + RCSB ModelServer subsets
  // ==================================================================

  async function showPdbeSummary() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    setBusy(true, 'Querying PDBe (EU mirror)...')
    const pane = el('structure-results')
    try {
      const data = await sendCommand('pdbe_entry_summary', {})
      if (pane) {
        const s = data.summary?.[0] || {}
        pane.innerHTML = `<div class="v2-summary">PDBe entry summary (${escapeHtml(data.source || '')})</div>
          <table><tbody>
            <tr><td>Title</td><td>${escapeHtml(s.title || '—')}</td></tr>
            <tr><td>Method</td><td>${escapeHtml((s.experiment_type || []).join(', ') || '—')}</td></tr>
            <tr><td>Resolution</td><td>${s.resolution ? Number(s.resolution).toFixed(2) + ' Å' : '—'}</td></tr>
            <tr><td>Released</td><td>${escapeHtml(s.release_date || '—')}</td></tr>
            <tr><td>Chains</td><td>${escapeHtml(String((s.in_chain_ids || []).join(', ') || '—'))}</td></tr>
          </tbody></table>`
      }
      showStatus('PDBe summary fetched')
    } catch (e) {
      if (pane) pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(e.message)}</div>`
      showStatus(`PDBe lookup failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  async function fetchModelserverSubset() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    const chain = await appPrompt(
      'Fetch a chain subset via RCSB ModelServer (mmCIF)',
      `label_asym_id (blank = whole entry, e.g. A)`)
    if (chain === null) return
    setBusy(true, 'Fetching subset from ModelServer...')
    const pane = el('structure-results')
    try {
      const data = await sendCommand('modelserver_subset', {
        label_asym_id: chain || undefined,
      })
      // Overlay the subset as its own model (blue carbons) so the fetched
      // coordinates are visible without replacing the loaded structure.
      if (viewer) {
        viewer.addModel(data.content, 'cif')
        const modelId = (viewer.models?.length || 2) - 1
        viewer.setStyle({ model: modelId }, {
          cartoon: { color: '#58a6ff' },
          stick: { radius: 0.2, colorscheme: 'blueCarbon' },
        })
        viewer.zoomTo({ model: modelId })
        viewer.render()
      }
      if (pane) {
        pane.innerHTML = `<div class="v2-summary">ModelServer subset loaded as overlay model: ${escapeHtml(data.url)}</div>`
      }
      showStatus('ModelServer subset rendered (blue)')
    } catch (e) {
      if (pane) pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(e.message)}</div>`
      showStatus(`ModelServer fetch failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  // ==================================================================
  // V2.1: pockets (real fpocket) + secondary structure (real DSSP)
  // ==================================================================

  async function detectPockets() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    setBusy(true, 'Running fpocket...')
    try {
      const data = await sendCommand('detect_pockets', {})
      renderPockets(data)
    } catch (e) {
      showStatus(`Pocket detection failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  function renderPockets(data) {
    const pane = el('pockets-results')
    if (!pane) return
    if (!data.available) {
      pane.innerHTML = `<div class="v2-placeholder">fpocket is not installed — ${escapeHtml(data.error || '')}. Nothing is simulated.</div>`
      return
    }
    if (data.error) {
      pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(data.error)}</div>`
      return
    }
    const overlaps = data.ligand_overlap || []
    let html = `<div class="v2-summary">${data.pocket_count} pocket(s) detected by fpocket</div>`
    html += '<table><thead><tr><th>Pocket</th><th>Atoms</th><th>Score</th><th>Druggability</th><th>Dist. to ligand Å</th><th>Lig. atoms in 4 Å</th><th></th></tr></thead><tbody>'
    for (const p of data.pockets) {
      const ov = overlaps.find(o => o.pocket_id === p.pocket_id) || {}
      html += `<tr>
        <td>${escapeHtml(p.pocket_id)}</td>
        <td>${p.atom_count}</td>
        <td>${p.descriptors?.Score ?? '—'}</td>
        <td>${p.descriptors?.['Druggability Score'] ?? '—'}</td>
        <td>${ov.min_distance != null ? Number(ov.min_distance).toFixed(2) : '—'}</td>
        <td>${ov.ligand_atoms_within_4A != null ? `${ov.ligand_atoms_within_4A}/${ov.ligand_atom_count}` : '—'}</td>
        <td><button class="btn btn-secondary btn-sm pocket-zoom" data-id="${escapeHtml(p.pocket_id)}">Zoom</button></td>
      </tr>`
    }
    html += '</tbody></table>'
    pane.innerHTML = html
    pane.querySelectorAll('.pocket-zoom').forEach(btn => {
      btn.addEventListener('click', () => {
        const p = data.pockets.find(x => x.pocket_id === btn.dataset.id)
        if (!p || !viewer) return
        // Render fpocket's own pocket atoms as a sphere cloud model.
        const lines = ['MODEL        1']
        let serial = 1
        for (const c of p.coordinates) {
          lines.push(`HETATM${String(serial).padStart(5, ' ')}  C   PKT P   1    ` +
            `${c[0].toFixed(3).padStart(8)}${c[1].toFixed(3).padStart(8)}${c[2].toFixed(3).padStart(8)}  1.00  0.00           C`)
          serial++
        }
        lines.push('ENDMDL', 'END')
        viewer.addModel(lines.join('\n'), 'pdb')
        viewer.setStyle({ model: 1 }, { sphere: { scale: 0.25, color: '#a371f7', opacity: 0.6 } })
        viewer.zoomTo({ model: 1 })
        viewer.render()
      })
    })
  }

  async function runDssp() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    setBusy(true, 'Running DSSP...')
    try {
      const data = await sendCommand('run_dssp', {})
      const pane = el('pockets-results')
      if (!pane) return
      if (!data.available) {
        pane.innerHTML = `<div class="v2-placeholder">DSSP is not installed — ${escapeHtml(data.error || '')}. Nothing is simulated.</div>`
        return
      }
      if (data.error) {
        pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(data.error)}</div>`
        return
      }
      const counts = data.ss_counts || {}
      const countStr = Object.entries(counts).map(([k, v]) => `${k}: ${v}`).join(' · ')
      pane.innerHTML = `<div class="v2-summary">DSSP assigned ${data.residues.length} residues — ${escapeHtml(countStr)}</div>
        <div class="v2-placeholder">Open the Structure tab → Sequence View to see the assignment mapped per residue.</div>`
      showStatus('DSSP done')
    } catch (e) {
      showStatus(`DSSP failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  // ==================================================================
  // V2.1: PDB-wide discovery (RCSB search API, live)
  // ==================================================================

  async function runDiscovery() {
    const mode = el('discover-mode')?.value || 'text'
    const input = el('discover-input')?.value?.trim() || ''
    const status = el('discover-status')
    const pane = el('discover-results')
    if (!pane) return
    try {
      let data
      if (mode === 'text') {
        if (!input) { showStatus('Enter search text', true); return }
        if (status) status.textContent = 'searching the PDB…'
        data = await sendCommand('discover_text', { text: input, rows: 25 })
      } else if (mode === 'sequence') {
        if (!input) { showStatus('Paste a FASTA sequence', true); return }
        if (status) status.textContent = 'running MMseqs2…'
        data = await sendCommand('discover_sequence', { sequence: input, rows: 25 })
      } else if (mode === 'chemical') {
        if (!input) { showStatus('Enter a SMILES string', true); return }
        if (status) status.textContent = 'screening chemistry…'
        data = await sendCommand('discover_chemical', { smiles: input, rows: 25 })
      } else if (mode === 'same_ligand') {
        if (status) status.textContent = 'finding same-ligand entries…'
        data = await sendCommand('discover_same_ligand', { rows: 25 })
      } else if (mode === 'similar_components') {
        if (status) status.textContent = 'finding similar components…'
        data = await sendCommand('discover_similar_components', {
          smiles: input || undefined, rows: 25,
        })
      }
      if (status) status.textContent = ''
      renderDiscovery(data)
    } catch (e) {
      if (status) status.textContent = ''
      pane.innerHTML = `<div class="v2-error v2-summary">${escapeHtml(e.message)}</div>`
      showStatus(`Discovery failed: ${e.message}`, true)
    }
  }

  function renderDiscovery(data) {
    const pane = el('discover-results')
    if (!pane) return
    const hits = data.hits || []
    let html = `<div class="v2-summary">${data.total_count} total hits (${escapeHtml(data.service || '')} · ${escapeHtml(data.source || '')}) — showing ${hits.length}</div>`
    if (hits.length === 0) {
      html += '<div class="v2-placeholder">No hits</div>'
    } else {
      html += '<table><thead><tr><th>ID</th><th>Score</th><th></th></tr></thead><tbody>'
      for (const h of hits) {
        const isComponent = /^[0-9a-z]{1,3}$/i.test(h.id) && !/^[0-9][a-z0-9]{3}$/i.test(h.id)
        html += `<tr class="search-hit">
          <td>${escapeHtml(h.id)}</td>
          <td>${h.score != null ? Number(h.score).toFixed(3) : '—'}</td>
          <td>${isComponent
            ? '<button class="btn btn-secondary btn-sm discover-open" data-kind="component" data-id="' + escapeHtml(h.id) + '">Enrich</button>'
            : '<button class="btn btn-primary btn-sm discover-open" data-kind="entry" data-id="' + escapeHtml(h.id) + '">Open</button>'}</td>
        </tr>`
      }
      html += '</tbody></table>'
    }
    pane.innerHTML = html
    pane.querySelectorAll('.discover-open').forEach(btn => {
      btn.addEventListener('click', () => {
        const id = btn.dataset.id
        if (btn.dataset.kind === 'entry') {
          openPdbIdValue(id)
        } else {
          enrichComponent(id)
        }
      })
    })
  }

  /** Fetch a CCD component's real identity from the live resolver. */
  async function enrichComponent(compId) {
    setBusy(true, `Fetching ${compId} from the CCD...`)
    try {
      const data = await sendCommand('enrich_component', { comp_id: compId })
      switchTab('discover')
      const pane = el('discover-results')
      if (pane) {
        pane.innerHTML = `
          <div class="v2-summary">${escapeHtml(data.name || compId)} (${escapeHtml(compId)}) — ${escapeHtml(data.pdbx_type || data.type || '')}</div>
          <div class="ligand-properties">
            ${ligandRow('Formula', data.formula)}
            ${ligandRow('MW', data.molecular_weight != null ? Number(data.molecular_weight).toFixed(2) : null)}
            ${ligandRow('SMILES', data.smiles)}
            ${ligandRow('InChI Key', data.inchi_key)}
          </div>`
      }
    } catch (e) {
      showStatus(`Enrichment failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  // ==================================================================
  // V2.1: docking box editing + multi-ligand queue + redocking compare
  // ==================================================================

  async function applyDockingBox() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    const vals = ['box-cx', 'box-cy', 'box-cz', 'box-sx', 'box-sy', 'box-sz']
      .map(id => parseFloat(el(id)?.value))
    if (vals.some(v => Number.isNaN(v))) {
      showStatus('Fill all six box fields (center x/y/z, size x/y/z)', true)
      return
    }
    try {
      await sendCommand('set_docking_box', {
        center: [vals[0], vals[1], vals[2]],
        size: [vals[3], vals[4], vals[5]],
      })
      showStatus('Docking box override applied')
    } catch (e) {
      showStatus(`Box override failed: ${e.message}`, true)
    }
  }

  async function resetDockingBox() {
    try {
      // Fetch the current box (user-set or ligand extent) and display it.
      const data = await sendCommand('get_docking_box', {})
      const b = data.box || {}
      const cx = el('box-cx'), cy = el('box-cy'), cz = el('box-cz')
      const sx = el('box-sx'), sy = el('box-sy'), sz = el('box-sz')
      if (cx && b.center) cx.value = b.center[0].toFixed(1)
      if (cy && b.center) cy.value = b.center[1].toFixed(1)
      if (cz && b.center) cz.value = b.center[2].toFixed(1)
      if (sx && b.size) sx.value = b.size[0].toFixed(1)
      if (sy && b.size) sy.value = b.size[1].toFixed(1)
      if (sz && b.size) sz.value = b.size[2].toFixed(1)
      showStatus(`Current box (${b.source || 'ligand_extent'}): center ${b.center?.map(v => Number(v).toFixed(1)).join(', ') || '—'}`)
    } catch (e) {
      showStatus(`Box fetch failed: ${e.message}`, true)
    }
  }

  async function queueDocking() {
    if (!state.structure) {
      showStatus('Open a structure first', true)
      return
    }
    const exhaustiveness = Math.max(1, parseInt(el('dock-exhaustiveness')?.value || '8', 10))
    const numModes = Math.max(1, parseInt(el('dock-num-modes')?.value || '5', 10))
    const progress = el('dock-progress')
    if (progress) progress.style.display = 'flex'
    const results = el('dock-results')
    if (results) results.innerHTML = '<div class="v2-placeholder">Queueing one real Vina job per non-solvent ligand…</div>'
    try {
      const started = await sendCommand('queue_docking', { exhaustiveness, num_modes: numModes })
      const queued = started.queued || []
      const failed = started.failed || []
      showStatus(`Queued ${queued.length} docking job(s)`)  
      // Poll until all terminal.
      const ids = queued.map(q => q.job_id)
      let done = 0
      while (done < ids.length) {
        await new Promise(r => setTimeout(r, 3000))
        done = 0
        for (const id of ids) {
          const job = await sendCommand('get_job', { job_id: id })
          if (job.status === 'completed' || job.status === 'failed') done++
        }
      }
      await renderQueueResults(ids)
    } catch (e) {
      if (results) results.innerHTML = `<div class="v2-error v2-summary">Queue failed: ${escapeHtml(e.message)}</div>`
      showStatus(`Queue failed: ${e.message}`, true)
    } finally {
      if (progress) progress.style.display = 'none'
    }
  }

  async function renderQueueResults(ids) {
    const results = el('dock-results')
    if (!results) return
    let html = '<table><thead><tr><th>Ligand</th><th>Status</th><th>Best affinity</th><th>Poses</th><th></th></tr></thead><tbody>'
    state.queueJobs = []
    for (const id of ids) {
      const job = await sendCommand('get_job', { job_id: id })
      state.queueJobs.push(job)
      const poses = job.result?.poses || []
      const best = poses.length ? Number(poses[0].affinity).toFixed(2) : '—'
      html += `<tr>
        <td>${escapeHtml(job.params?.ligand_name || id.slice(0, 8))}</td>
        <td><span class="job-status ${job.status === 'completed' ? 'completed' : 'failed'}"></span> ${escapeHtml(job.status)}</td>
        <td>${best}</td>
        <td>${poses.length}</td>
        <td>${job.status === 'completed' && poses.length
          ? `<button class="btn btn-secondary btn-sm queue-pose" data-job="${escapeHtml(id)}">Show in 3D</button>`
          : `<span class="v2-error" title="${escapeHtml(job.error || '')}">${job.status === 'failed' ? 'failed' : ''}</span>`}</td>
      </tr>`
    }
    html += '</tbody></table>'
    results.innerHTML = html
    results.querySelectorAll('.queue-pose').forEach(btn => {
      btn.addEventListener('click', () => {
        const job = state.queueJobs.find(j => j.job_id === btn.dataset.job)
        if (job) {
          state.lastDockJob = job
          showPoseInViewer(job.result.poses[0].pose_id)
        }
      })
    })
    // Enable the redocking comparison now that jobs exist.
    const row = el('crystal-dock-row')
    if (row) row.style.display = 'flex'
  }

  async function compareCrystalDocked() {
    const job = state.lastDockJob || (state.queueJobs || [])[0]
    if (!job) {
      showStatus('Run docking first', true)
      return
    }
    try {
      const data = await sendCommand('compare_crystal_docked', {
        job_id: job.job_id,
        pose_id: job.result?.poses?.[0]?.pose_id,
      })
      const pane = el('crystal-dock-results')
      if (pane) {
        pane.innerHTML = data.error
          ? `<div class="v2-error v2-summary">${escapeHtml(data.error)}</div>`
          : `<div class="v2-summary">Redocking RMSD (heavy atoms, element-aware Hungarian + Kabsch): ${data.rmsd} Å over ${data.atom_count} atoms</div>
             <div class="v2-placeholder">${escapeHtml(data.method || '')}</div>`
      }
      showStatus(data.error ? 'Comparison refused' : `Redocking RMSD ${data.rmsd} Å`)
    } catch (e) {
      showStatus(`Redocking compare failed: ${e.message}`, true)
    }
  }

  // ==================================================================
  // V2.1: structure superposition
  // ==================================================================

  async function runSuperposition() {
    const chainA = el('superpose-chain-a')?.value?.trim()
    const ref = el('superpose-ref')?.value?.trim()
    const chainB = el('superpose-chain-b')?.value?.trim()
    if (!chainA || !ref || !chainB) {
      showStatus('Fill chain (this structure), PDB ID and chain (reference)', true)
      return
    }
    setBusy(true, `Fetching ${ref.toUpperCase()} and superposing...`)
    try {
      const data = await sendCommand('superpose_structures', {
        chain_a: chainA, reference_pdb_id: ref.toUpperCase(), chain_b: chainB,
      })
      const pane = el('superpose-results')
      if (pane) {
        pane.innerHTML = data.error
          ? `<div class="v2-error v2-summary">${escapeHtml(data.error)}</div>`
          : `<div class="v2-summary">RMSD ${data.rmsd} Å over ${data.residues_used} core residues (${data.nw_aligned_pairs} aligned, core fraction ${data.core_fraction})</div>
             <div class="v2-placeholder">Sequence identity ${data.sequence_identity ?? '—'} · method: NW alignment + seed-and-extend core fit, cutoff ${data.refine_cutoff} Å</div>`
      }
      showStatus(data.error ? 'Superposition refused' : `Superposition RMSD ${data.rmsd} Å`)
    } catch (e) {
      showStatus(`Superposition failed: ${e.message}`, true)
    } finally {
      setBusy(false)
    }
  }

  // ==================================================================
  // V2.1: keyboard ligand switching + recent ligands
  // ==================================================================

  const recentLigands = []

  function rememberLigand(ligandId) {
    if (!ligandId) return
    const i = recentLigands.indexOf(ligandId)
    if (i >= 0) recentLigands.splice(i, 1)
    recentLigands.unshift(ligandId)
    if (recentLigands.length > 5) recentLigands.pop()
  }

  function ligandIdsInOrder() {
    return (state.structure?.ligands || [])
      .filter(l => l.classification_hint && !['HETAS', 'HETAI'].includes(l.classification_hint))
      .map(l => l.id)
  }

  function switchLigand(direction) {
    const ids = ligandIdsInOrder()
    if (ids.length === 0) return
    const idx = ids.indexOf(state.selectedLigandId)
    const next = idx < 0 ? 0 : (idx + direction + ids.length) % ids.length
    selectLigand(ids[next])
  }

  function onGlobalKeydown(e) {
    // Ignore typing contexts.
    const tag = (e.target?.tagName || '').toLowerCase()
    if (tag === 'input' || tag === 'textarea' || tag === 'select' || e.target?.isContentEditable) return
    if (e.key === '[') switchLigand(-1)
    if (e.key === ']') switchLigand(1)
  }

  // ==================================================================
  // V2: tab switching
  // ==================================================================

  function initV2Tabs() {
    const dock = el('v2-dock')
    if (!dock) return
    dock.querySelectorAll('.v2-tab').forEach(tab => {
      tab.addEventListener('click', () => switchTab(tab.dataset.tab))
    })
    dock.querySelector('.v2-tab')?.classList.add('active')
    el('v2-panel-batch')?.classList.remove('v2-hidden')
  }

  function switchTab(name) {
    const dock = el('v2-dock')
    if (!dock) return
    const tab = dock.querySelector(`.v2-tab[data-tab="${name}"]`)
    if (!tab) return
    dock.querySelectorAll('.v2-tab').forEach(t => t.classList.remove('active'))
    dock.querySelectorAll('.v2-panel').forEach(p => p.classList.add('v2-hidden'))
    tab.classList.add('active')
    el(`v2-panel-${name}`)?.classList.remove('v2-hidden')
    if (name === 'compare') refreshComparisonJobs()
    if (name === 'editor') drawEditorCanvas()
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
    // Keyboard: [ and ] switch ligands (never while typing).
    document.addEventListener('keydown', onGlobalKeydown)
  }

  async function start() {
    initViewer()

    el('btn-open-file')?.addEventListener('click', openLocalFile)
    el('btn-open-pdb')?.addEventListener('click', openPdbId)
    el('btn-run-analysis')?.addEventListener('click', runAnalysis)
    el('btn-export')?.addEventListener('click', exportArtifacts)
    el('btn-measure')?.addEventListener('click', toggleMeasure)
    el('btn-get-status')?.addEventListener('click', refreshStatus)

    // Spec V1 §5.1 top bar: search, recent, settings.
    el('btn-topbar-search')?.addEventListener('click', topbarSearch)
    el('topbar-search')?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') topbarSearch()
    })
    el('select-recent')?.addEventListener('change', (e) => {
      if (e.target.value) openPdbIdValue(e.target.value)
    })
    el('btn-settings')?.addEventListener('click', showSettings)

    // Spec V1 §4.4/§5.2/§4.6: pocket surface, camera presets, scene PNG.
    el('btn-pocket-surface')?.addEventListener('click', () => {
      if (!state.pocketResidues?.length) {
        showStatus('Run contact analysis first (it computes the pocket)', true)
        return
      }
      state.pocketSurfaceOn = !state.pocketSurfaceOn
      el('btn-pocket-surface')?.classList.toggle('active', state.pocketSurfaceOn)
      applyView()
      showStatus(`Pocket surface ${state.pocketSurfaceOn ? 'on' : 'off'}`)
    })
    el('select-camera-preset')?.addEventListener('change', (e) => {
      if (e.target.value) {
        applyCameraPreset(e.target.value)
        e.target.value = ''
      }
    })
    el('btn-scene-image')?.addEventListener('click', exportSceneImage)
    el('btn-contacts-csv')?.addEventListener('click', exportContactsCsv)

    // Spec V2: editor sync, clustering, PDBe/ModelServer.
    el('btn-editor-sync')?.addEventListener('click', editorSyncTo3d)
    el('btn-cluster-run')?.addEventListener('click', runClustering)
    el('btn-pdbe-run')?.addEventListener('click', showPdbeSummary)
    el('btn-modelserver-run')?.addEventListener('click', fetchModelserverSubset)

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
    el('btn-dock-run')?.addEventListener('click', runDocking)
    el('btn-dock-queue')?.addEventListener('click', queueDocking)
    el('btn-box-apply')?.addEventListener('click', applyDockingBox)
    el('btn-box-reset')?.addEventListener('click', resetDockingBox)
    el('btn-crystal-dock')?.addEventListener('click', compareCrystalDocked)
    el('btn-superpose-run')?.addEventListener('click', runSuperposition)
    el('btn-discover-run')?.addEventListener('click', runDiscovery)
    el('discover-input')?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') runDiscovery()
    })
    el('btn-pockets-run')?.addEventListener('click', detectPockets)
    el('btn-dssp-run')?.addEventListener('click', runDssp)
    el('btn-diagram-run')?.addEventListener('click', showInteractionDiagram)
    el('btn-sequence-run')?.addEventListener('click', showSequenceView)
    el('btn-quality-run')?.addEventListener('click', showEntryQuality)
    el('btn-uniprot-run')?.addEventListener('click', showUniprotContext)
    el('btn-alphafold-run')?.addEventListener('click', showAlphafoldModel)
    el('btn-bindingdb-run')?.addEventListener('click', showBindingDB)
    el('btn-covale-run')?.addEventListener('click', showCovalentLinks)
    el('btn-md-run')?.addEventListener('click', runMd)
    el('btn-manifest-run')?.addEventListener('click', exportReproManifest)
    el('btn-search-run')?.addEventListener('click', runSearch)
    el('btn-search-refresh')?.addEventListener('click', refreshSearchIndex)
    el('search-input')?.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') runSearch()
    })

    await refreshStatus()
  }

  return { upgrade, start }
}
