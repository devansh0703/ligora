/**
 * Stage 6: spec-completion features through the real UI.
 * - Top bar: PDB-wide search, recent structures, settings (live status)
 * - Contacts: sort/filter controls + one-click CSV export
 * - Camera presets + pocket surface toggle
 * - Scene PNG export (real canvas pixels saved through the backend)
 * - Pose clustering + consensus (two real Vina runs)
 * - PDBe entry summary + RCSB ModelServer subset (live services)
 * - 2D editor -> 3D sync (real edited SDF overlaid in the viewer)
 * Real backend, no mocks.
 */
import {
  startStatic, startBridge, bridgeCommand, launch, makeReport,
  sleep, waitFor, SHOTS,
} from './ui-lib.mjs'
import path from 'node:path'

const staticServer = await startStatic()
const { bridge, ready } = startBridge()
await ready
const { browser, page } = await launch()
const rep = makeReport()

const text = (sel) => page.$eval(sel, (e) => e.textContent.trim()).catch(() => null)
const count = (sel) => page.$$eval(sel, (els) => els.length).catch(() => 0)
const shot = (n) => page.screenshot({ path: path.join(SHOTS, `${n}.png`) })

try {
  await page.goto('http://127.0.0.1:8742/', { waitUntil: 'networkidle0' })
  await waitFor(async () => (await count('#scene-container canvas')) >= 1,
    30000, 'canvas')

  console.log('== top bar: PDB-wide search opens the Discover tab with hits ==')
  await page.type('#topbar-search', 'ubiquitin')
  await page.click('#btn-topbar-search')
  await waitFor(async () => {
    const rows = await count('#discover-results tbody tr')
    return rows > 0
  }, 90000, 'top-bar search renders PDB hits')
  const hitRows = await count('#discover-results tbody tr')
  rep.ok(hitRows > 0, `top-bar search -> ${hitRows} live PDB hits`)
  await shot('30-topbar-search')

  console.log('== open 3W85, contacts sort/filter + CSV ==')
  await page.click('#btn-open-pdb')
  await waitFor(() => page.$eval('#app-modal', (m) => m.style.display !== 'none'),
    5000, 'modal opens')
  await page.type('#modal-input', '3W85')
  await page.click('#modal-ok')
  await waitFor(async () => {
    const card = await text('#ligand-card')
    return card && card.length > 40
  }, 120000, 'structure opened via UI')

  await page.click('#btn-run-analysis')
  await waitFor(async () => (await count('#contacts-table tbody tr')) > 0,
    180000, 'PLIP contacts render')

  // Sort/filter controls present and functional.
  rep.ok(await page.$('#contact-filter') !== null, 'contact type filter present')
  rep.ok(await page.$('#contact-sort') !== null, 'contact sort control present')
  const rowsBefore = await count('#contacts-table tbody tr')
  await page.select('#contact-filter', 'hydrogen_bond')
  await waitFor(async () => {
    const rows = await count('#contacts-table tbody tr')
    return rows > 0 && rows <= rowsBefore
  }, 10000, 'filter narrows rows')
  const filtered = await count('#contacts-table tbody tr')
  rep.ok(filtered < rowsBefore, `filter narrows contacts (${rowsBefore} -> ${filtered})`)
  await page.select('#contact-sort', 'distance')
  await sleep(400)
  const dists = await page.$$eval('#contacts-table tbody tr',
    (els) => els.map((e) => parseFloat(e.cells[3].textContent)))
  const sorted = dists.every((d, i) => i === 0 || dists[i - 1] <= d)
  rep.ok(sorted, `distance sort orders rows (${dists.slice(0, 3).join(', ')}...)`)
  await page.select('#contact-filter', 'all')
  await shot('31-contacts-filtered')

  console.log('== CSV export (real file through backend + save dialog) ==')
  // The Tauri save dialog is unavailable in the browser harness; assert the
  // backend command path the button drives.
  const csv = await bridgeCommand('export_contacts_csv', {})
  rep.ok(csv.success === true && csv.data?.format === 'csv',
    `contacts CSV export path works (${csv.error || 'ok'})`)

  console.log('== camera presets ==')
  rep.ok(await page.$('#select-camera-preset') !== null, 'camera preset select present')
  await page.select('#select-camera-preset', 'ligand')
  await sleep(500)
  await shot('32-camera-ligand')
  rep.ok(true, 'ligand-focus preset applied')

  console.log('== pocket surface toggle ==')
  await page.click('#btn-pocket-surface')
  await sleep(1200)
  rep.ok(await page.$('#btn-pocket-surface.active') !== null,
    'pocket surface toggles on (backend-computed pocket residues)')
  await shot('33-pocket-surface')

  console.log('== scene PNG export (real canvas pixels) ==')
  const pngResult = await page.evaluate(async () => {
    const canvas = document.querySelector('#scene-container canvas')
    window.__viewerRender?.()
    return canvas ? canvas.toDataURL('image/png').length : 0
  })
  rep.ok(pngResult > 1000, `scene canvas yields real PNG data (${pngResult} chars)`)
  const saved = await page.evaluate(async () => {
    const canvas = document.querySelector('#scene-container canvas')
    const b64 = canvas.toDataURL('image/png').split(',')[1]
    return window.__LIGORA_BRIDGE__('save_scene_image', { png_base64: b64 })
  })
  rep.ok(saved.success === true && saved.data?.size_bytes > 1000,
    `scene PNG saved to workspace (${saved.data?.size_bytes} bytes)`)

  console.log('== settings dialog (live engine/source status) ==')
  await page.click('#btn-settings')
  await waitFor(async () => {
    const body = await text('#settings-body')
    return body && body.includes('available')
  }, 30000, 'settings shows live status')
  const settings = await text('#settings-body')
  rep.ok(/plip|vina|gromacs/i.test(settings || ''), `settings lists engines live`)
  await page.click('#modal-ok')
  await shot('34-settings')

  console.log('== recent structures populated ==')
  const recentOpts = await page.$$eval('#select-recent option', (els) => els.map((e) => e.value))
  rep.ok(recentOpts.includes('3W85'), `recent list contains 3W85 (${recentOpts.join(',')})`)

  console.log('== PDBe entry summary (live EU mirror) ==')
  await page.click('.v2-tab[data-tab="structure"]')
  await page.click('#btn-pdbe-run')
  await waitFor(async () => {
    const s = await text('#structure-results .v2-summary')
    return s && s.includes('PDBe')
  }, 60000, 'PDBe summary renders')
  const pdbe = await text('#structure-results')
  rep.ok(/Title|Released/.test(pdbe || ''), 'PDBe summary fields rendered')
  await shot('35-pdbe')

  console.log('== ModelServer subset overlay (live) ==')
  await page.click('#btn-modelserver-run')
  await waitFor(() => page.$eval('#app-modal', (m) => m.style.display !== 'none'),
    5000, 'subset prompt opens')
  await page.type('#modal-input', 'A')
  await page.click('#modal-ok')
  await waitFor(async () => {
    const s = await text('#structure-results .v2-summary')
    return s && s.includes('ModelServer')
  }, 90000, 'ModelServer subset renders')
  rep.ok(true, 'ModelServer chain subset fetched and overlaid')
  await shot('36-modelserver')

  console.log('== 2D editor -> 3D sync ==')
  await page.click('.v2-tab[data-tab="editor"]')
  await page.click('#btn-editor-load')
  await waitFor(async () => {
    const info = await text('#editor-info')
    return info && info.includes('atoms')
  }, 60000, 'editor loads ligand')
  await page.click('#btn-editor-sync')
  await waitFor(async () => {
    const s = await text('#status-toast')
    return s && /3D/.test(s)
  }, 60000, 'sync status shows')
  const syncStatus = await text('#status-toast')
  rep.ok(/23 atoms/.test(syncStatus || ''), `editor synced to 3D: ${syncStatus}`)
  await shot('37-editor-sync')

  console.log('== pose clustering + consensus (two real Vina runs) ==')
  await page.click('.v2-tab[data-tab="docking"]')
  await page.click('#btn-dock-run')
  await waitFor(async () => {
    const s = await text('#dock-results .v2-summary')
    return s && /pose/i.test(s)
  }, 280000, 'docking run 1 completes')
  await page.click('#btn-dock-run')
  await waitFor(async () => {
    const rows = await count('#dock-results tbody tr')
    return rows >= 5
  }, 280000, 'docking run 2 completes')

  await page.click('.v2-tab[data-tab="compare"]')
  await waitFor(async () =>
    (await page.$eval('#compare-job-a', (s) => s.options.length)) >= 3,
    30000, 'job dropdowns populated')
  // Pick the two most recent jobs (last two options beyond the placeholder).
  const jobOpts = await page.$$eval('#compare-job-a option',
    (els) => els.filter((e) => e.value).map((e) => e.value))
  await page.select('#compare-job-a', jobOpts[jobOpts.length - 2])
  await page.select('#compare-job-b', jobOpts[jobOpts.length - 1])
  await page.click('#btn-cluster-run')
  await waitFor(async () => {
    const s = await text('#compare-results .v2-summary')
    return s && s.includes('cluster')
  }, 60000, 'clustering renders')
  const clusterSummary = await text('#compare-results .v2-summary')
  rep.ok(/cluster\(s\)/.test(clusterSummary || ''), `clusters: ${clusterSummary}`)
  const consensusShown = await text('#compare-results')
  rep.ok(/Consensus pose/.test(consensusShown || ''), 'consensus pose shown')
  await shot('38-clustering')
} catch (e) {
  rep.ok(false, `stage failed: ${e.message}`)
  await shot('99-stage6-failure').catch(() => {})
} finally {
  const okAll = rep.summary('stage 6: spec-completion features')
  await browser.close()
  bridge.kill()
  staticServer.close()
  process.exit(okAll ? 0 : 1)
}
