/**
 * Stage 4: different real data through the real UI - verifies the app
 * generalizes beyond 3W85. Real backend, real data, no mocks.
 *
 * - 1MBN (myoglobin): heme (HEM) must be the picked ligand, never water
 *   or the hydroxide ion - false-positive guard on different data.
 * - 1CRN (crambin, protein-only): honest "no drug-like ligand" path, no crash.
 * - Search UI: real BM25 queries over live RCSB/PubChem documents.
 */
import {
  startStatic, startBridge, launch, makeReport, sleep, waitFor, SHOTS,
} from './ui-lib.mjs'
import path from 'node:path'

const staticServer = await startStatic()
const { bridge, ready } = startBridge()
await ready
const { browser, page } = await launch()
const rep = makeReport()

const text = (sel) => page.$eval(sel, (e) => e.textContent.trim()).catch(() => null)
const exists = (sel) => page.$(sel).then((e) => !!e)
const count = (sel) => page.$$eval(sel, (els) => els.length).catch(() => 0)
const shot = (n) => page.screenshot({ path: path.join(SHOTS, `${n}.png`) })

async function openViaModal(pdbId) {
  await page.click('#btn-open-pdb')
  await waitFor(() => page.$eval('#app-modal', (m) => m.style.display !== 'none'),
    5000, 'modal opens')
  await page.evaluate(() => { document.getElementById('modal-input').value = '' })
  await page.type('#modal-input', pdbId)
  await page.click('#modal-ok')
}

try {
  console.log('== boot ==')
  await page.goto('http://127.0.0.1:8742/', { waitUntil: 'networkidle0' })
  await waitFor(async () => (await count('#scene-container canvas')) >= 1,
    30000, '3Dmol canvas')

  console.log('== 1MBN: heme picked, not water/ions ==')
  await openViaModal('1MBN')
  // Ligand card must populate from live CCD with the heme.
  await waitFor(async () => {
    const card = await text('#ligand-card')
    return card && /heme|protoporphyrin|HEM/i.test(card)
  }, 120000, '1MBN ligand card shows heme identity from live CCD')
  const mbnCard = await text('#ligand-card')
  rep.ok(/HEM|heme|protoporphyrin/i.test(mbnCard), '1MBN picks HEM (real cofactor)')
  rep.ok(!/\bHOH\b/.test(mbnCard), 'water never presented as the ligand (1MBN)')
  const mbnStatus = await text('#status-toast')
  rep.ok(!/no drug-like ligand/i.test(mbnStatus || ''), 'no false no-ligand message (1MBN)')
  await shot('41-1mbn-heme')

  console.log('== 1MBN: PLIP contacts against heme ==')
  await page.click('#btn-run-analysis')
  await waitFor(async () => (await count('#contacts-table tbody tr')) > 0,
    180000, 'PLIP returns contacts for heme')
  const mbnContacts = await count('#contacts-table tbody tr')
  rep.ok(mbnContacts > 0, `PLIP finds real heme contacts (${mbnContacts} rows)`)
  await shot('42-1mbn-contacts')

  console.log('== search UI: live BM25 over real documents ==')
  // Run while 1MBN is still the active session: "Index current session"
  // fetches real docs for the session's own structure + ligands.
  await page.click('.v2-tab[data-tab="search"]')
  await waitFor(async () => !(await page.$eval('#v2-panel-search',
    (e) => e.classList.contains('v2-hidden'))), 5000, 'search panel visible')

  await page.click('#btn-search-refresh')
  await waitFor(async () => {
    const s = await text('#search-status')
    return s && /docs indexed/i.test(s)
  }, 120000, 'session documents indexed from live sources')

  await page.evaluate(() => { document.getElementById('search-input').value = '' })
  await page.type('#search-input', 'myoglobin')
  await page.click('#btn-search-run')
  await waitFor(async () => (await count('#search-results tbody tr')) > 0,
    120000, 'search "myoglobin" returns live hits')
  const hits1 = await count('#search-results tbody tr')
  const first1 = await text('#search-results tbody tr')
  rep.ok(hits1 > 0, `BM25 search "myoglobin" -> ${hits1} hits`)
  rep.ok(/myoglobin|MYO|1MBN|MBN/i.test(first1), `top hit relevant: ${first1.slice(0, 60)}`)
  await shot('44-search-myoglobin')

  // Chemical-identity query: CCD name "PROTOPORPHYRIN IX CONTAINING FE"
  // must surface the HEM doc for this session's live-indexed chemicals.
  await page.evaluate(() => { document.getElementById('search-input').value = '' })
  await page.type('#search-input', 'protoporphyrin')
  await page.click('#btn-search-run')
  await waitFor(async () => {
    const rows = await page.$$eval('#search-results tbody tr',
      (els) => els.map((e) => e.textContent)).catch(() => [])
    return rows.some((r) => /HEM/i.test(r))
  }, 120000, 'chemical query resolves to HEM')
  const rows2 = await page.$$eval('#search-results tbody tr',
    (els) => els.map((e) => e.textContent))
  rep.ok(rows2.some((r) => /HEM/i.test(r)), 'chemical query "protoporphyrin" -> HEM doc')
  await shot('45-search-hem')

  console.log('== 1CRN: protein-only, honest no-ligand path ==')
  await openViaModal('1CRN')
  await waitFor(async () => {
    const s = await text('#status-toast')
    return s && /no drug-like ligand/i.test(s)
  }, 120000, 'honest no-ligand status for 1CRN')
  const crnStatus = await text('#status-toast')
  rep.ok(/no drug-like ligand/i.test(crnStatus), '1CRN reports no drug-like ligand honestly')
  const crnCard = await text('#ligand-card')
  rep.ok(!crnCard || !/HOH|heme/i.test(crnCard), 'no ligand card fabricated for 1CRN')
  await shot('43-1crn-no-ligand')

  const allPassed = rep.summary('stage 4: different data + search UI')
  if (!allPassed) process.exitCode = 1
} catch (e) {
  console.error('STAGE ERROR:', e.message)
  await shot('99-stage4-error').catch(() => {})
  process.exitCode = 1
} finally {
  await browser.close().catch(() => {})
  staticServer.close()
  bridge.kill()
}
