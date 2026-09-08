/**
 * Stage 4: v2.2 features through the real UI.
 * - Covalent links (1HRC: real deposited thioether links to HEC)
 * - BindingDB affinities (live REST service)
 * - Reproducibility manifest export
 * - Viewer-only file types (SDF/xyz) open in 3Dmol.js without faking chemistry
 * Real backend, no mocks.
 */
import {
  startStatic, startBridge, bridgeCommand, launch, makeReport,
  sleep, waitFor, SHOTS,
} from './ui-lib.mjs'
import path from 'node:path'
import { mkdtempSync, writeFileSync } from 'node:fs'
import os from 'node:os'

const staticServer = await startStatic()
const { bridge, ready } = startBridge()
await ready
const { browser, page } = await launch()
const rep = makeReport()

const text = (sel) => page.$eval(sel, (e) => e.textContent.trim()).catch(() => null)
const shot = (n) => page.screenshot({ path: path.join(SHOTS, `${n}.png`) })  // A real caffeine molecule in SDF and a real acetic-acid molecule in
  // xyz (converted from PubChem's own 3D SDF with Open Babel), both from
  // real sources — not fabricated files.
  async function fetchReal(name, url) {
    const res = await fetch(url)
    if (!res.ok) throw new Error(`fetch ${name}: ${res.status}`)
    return res.text()
  }

try {
  await page.goto('http://127.0.0.1:8742/', { waitUntil: 'networkidle0' })
  await waitFor(async () => {
    return (await page.$$eval('#scene-container canvas', (els) => els.length)) >= 1
  }, 30000, 'canvas')

  console.log('== open 1HRC (real covalent heme entry) ==')
  await page.click('#btn-open-pdb')
  await waitFor(() => page.$eval('#app-modal', (m) => m.style.display !== 'none'),
    5000, 'modal opens')
  await page.type('#modal-input', '1HRC')
  await page.click('#modal-ok')
  await waitFor(async () => {
    const card = await text('#ligand-card')
    return card && card.length > 40
  }, 120000, 'structure opened via UI')

  console.log('== Covalent links through the UI ==')
  await page.click('.v2-tab[data-tab="structure"]')
  await waitFor(async () => !(await page.$eval('#v2-panel-structure',
    (e) => e.classList.contains('v2-hidden'))), 5000, 'structure panel visible')
  await page.click('#btn-covale-run')
  await waitFor(async () => {
    const s = await text('#structure-results .v2-summary')
    return s && s.includes('covalent link')
  }, 60000, 'covalent links render')
  const cov = await text('#structure-results .v2-summary')
  rep.ok(/covalent link\(s\) declared/.test(cov || ''), `covalent summary: ${cov}`)
  const covRows = await page.$$eval('#structure-results tbody tr',
    (els) => els.map((e) => e.textContent))
  rep.ok(covRows.some((r) => r.includes('covale') && r.includes('HEC')),
    `real thioether link row to HEC present (${covRows.length} rows)`)
  await shot('20-covale')

  console.log('== BindingDB affinities through the UI (live service) ==')
  await page.click('#btn-bindingdb-run')
  await waitFor(async () => {
    const s = await text('#structure-results .v2-summary')
    return s && (s.includes('binding record') || s.includes('failed'))
  }, 90000, 'BindingDB pane renders')
  const bdb = await text('#structure-results .v2-summary')
  rep.ok(!/failed/i.test(bdb || ''), `BindingDB responded: ${bdb?.slice(0, 90)}`)
  const bdbRows = await page.$$eval('#structure-results tbody tr', (e) => e.length)
  rep.ok(bdbRows >= 1, `binding records rendered (${bdbRows} rows)`)
  await shot('21-bindingdb')

  console.log('== Reproducibility manifest through the UI ==')
  await page.click('#btn-manifest-run')
  await waitFor(async () => {
    const s = await text('#structure-results .v2-summary')
    return s && s.includes('Reproducibility manifest written')
  }, 60000, 'manifest summary renders')
  const man = await text('#structure-results .v2-summary')
  rep.ok(/artifact\(s\)/.test(man || ''), `manifest written: ${man}`)
  const manPath = await page.$eval('#structure-results .v2-placeholder',
    (e) => e.textContent.trim())
  rep.ok(/sha256: [0-9a-f]{64}/.test(manPath || ''), 'manifest sha256 shown')
  await shot('22-manifest')

  console.log('== Viewer-only file types (real PubChem SDF + xyz) ==')
  const tmp = mkdtempSync(path.join(os.tmpdir(), 'ligora-ui-files-'))
  const sdf = await fetchReal('caffeine.sdf',
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2519/SDF')
  writeFileSync(path.join(tmp, 'caffeine.sdf'), sdf)
  const sdf3d = await fetchReal('acetic.sdf',
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/176/SDF?record_type=3d')
  const xyz = await fetchReal('acetic.xyz',
    'https://files.rcsb.org/pubchem/ACM/ACM_3D.sdf')
    .catch(() => null)
  if (xyz == null) {
    // Convert PubChem's real 3D SDF to xyz with the app's own Open Babel
    // integration (real conversion, no fabrication).
    const { execFileSync } = await import('node:child_process')
    writeFileSync(path.join(tmp, 'acetic.sdf'), sdf3d)
    execFileSync('obabel', [path.join(tmp, 'acetic.sdf'), '-Oxyz',
      '-O', path.join(tmp, 'acetic.xyz')])
  } else {
    writeFileSync(path.join(tmp, 'acetic.xyz'), xyz)
  }

  // Open the SDF through the app's local-file command (the UI's open-local
  // flow uses the native file dialog, which puppeteer can't drive; the
  // backend command is exactly what the dialog handler invokes).
  const opened = await bridgeCommand('open_local_file',
    { file_path: path.join(tmp, 'caffeine.sdf') })
  rep.ok(opened.success === true && opened.data?.viewer_only === true,
    `SDF opens as viewer-only (${opened.error || 'ok'})`)
  await waitFor(async () =>
    (await page.$$eval('#scene-container canvas', (els) => els.length)) >= 1,
    30000, 'SDF renders in viewer')
  const st = await bridgeCommand('get_status', {})
  rep.ok(st.data?.session?.viewer_only === true,
    'status reports viewer_only session')
  // Analysis must honestly decline rather than approximate chemistry.
  const contacts = await bridgeCommand('run_contact_analysis', {})
  rep.ok(contacts.success === false || /no parsed structure|viewer-only/i
    .test(contacts.error || ''), 'analysis honestly declines for viewer-only file')
  await shot('23-sdf-viewer')

  const openedXyz = await bridgeCommand('open_local_file',
    { file_path: path.join(tmp, 'acetic.xyz') })
  rep.ok(openedXyz.success === true && openedXyz.data?.viewer_only === true,
    `xyz opens as viewer-only (${openedXyz.error || 'ok'})`)
  await shot('24-xyz-viewer')
} catch (e) {
  rep.ok(false, `stage failed: ${e.message}`)
  await shot('99-stage4-failure').catch(() => {})
} finally {
  const okAll = rep.summary('stage 4: v2.2 features + file types')
  await browser.close()
  bridge.kill()
  staticServer.close()
  process.exit(okAll ? 0 : 1)
}
