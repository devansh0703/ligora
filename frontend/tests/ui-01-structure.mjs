/**
 * Stage 1: boot, open real structure via modal, live enrichment, PLIP
 * contacts through the real UI. Real backend, real data, no mocks.
 */
import {
  startStatic, startBridge, launch, makeReport, sleep,
  waitFor, SHOTS,
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

try {
  console.log('== boot ==')
  await page.goto('http://127.0.0.1:8742/', { waitUntil: 'networkidle0' })
  await waitFor(async () => (await count('#scene-container canvas')) >= 1,
    30000, '3Dmol canvas')
  rep.ok(await exists('#v2-dock'), 'V2 dock present')
  rep.ok(await exists('#batch-sources'), 'batch input present')
  await shot('01-app-boot')

  console.log('== open 3W85 through the real modal ==')
  await page.click('#btn-open-pdb')
  await waitFor(() => page.$eval('#app-modal', (m) => m.style.display !== 'none'),
    5000, 'modal opens')
  rep.ok(await exists('#modal-input'), 'real modal dialog appears (no browser prompt)')
  await page.type('#modal-input', '3W85')
  await shot('02-modal')
  await page.click('#modal-ok')
  // The card renders the resolved ligand name (live CCD), plus a select of
  // the structure's ligands. Wait for either to appear.
  await waitFor(async () => {
    const card = await text('#ligand-card')
    const select = await page.$$eval('#ligand-select option',
      (os) => os.map((o) => o.textContent)).catch(() => [])
    return (card && card.length > 40) ||
      select.some((s) => s.includes('W85'))
  }, 120000, 'ligand card populated from live CCD')

  console.log('== live enrichment assertions ==')
  const name = await text('.ligand-name')
  rep.ok(name && name.length > 3, `ligand name from live CCD: "${name?.slice(0, 45)}"`)
  const badge = await text('.ligand-badge')
  rep.ok(!!badge, `classification badge: ${badge}`)
  const props = await page.$$eval('.ligand-property', (els) =>
    els.map((e) => e.textContent))
  rep.ok(props.some((p) => p.includes('73167555')),
    'PubChem CID 73167555 rendered (live PubChem)')
  rep.ok(props.some((p) => p.includes('InChI Key')) &&
    !props.some((p) => p.includes('—') && p.includes('InChI')),
    'InChI key populated from live sources')
  await shot('03-structure-enriched')

  console.log('== real PLIP contact analysis via UI button ==')
  await page.click('#btn-run-analysis')
  await waitFor(async () => (await count('#contacts-table tbody tr')) > 0,
    200000, 'contacts table populated')
  const rows = await count('#contacts-table tbody tr')
  rep.ok(rows >= 10, `PLIP contacts rendered (${rows} rows)`)
  const badge0 = await page.$eval('#contacts-table tbody tr td span',
    (e) => e.textContent)
  rep.ok(!!badge0, `typed contact badge present (${badge0})`)
  await shot('04-contacts')

  console.log('== evidence pane (live sources) ==')
  const ev = await count('#evidence-pane .evidence-item')
  rep.ok(ev > 0, `evidence items rendered (${ev})`)
  await shot('05-evidence')
} catch (e) {
  rep.ok(false, `stage failed: ${e.message}`)
  await shot('99-stage1-failure').catch(() => {})
} finally {
  const okAll = rep.summary('stage 1: structure + enrichment + contacts')
  await browser.close()
  bridge.kill()
  staticServer.close()
  process.exit(okAll ? 0 : 1)
}
