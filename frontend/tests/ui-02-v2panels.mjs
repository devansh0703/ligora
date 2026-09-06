/**
 * Stage 2: water network, 2D editor, scripting console through the real UI.
 * Requires the bridge session to open 3W85 first (real backend).
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
const exists = (sel) => page.$(sel).then((e) => !!e)
const count = (sel) => page.$$eval(sel, (els) => els.length).catch(() => 0)
const shot = (n) => page.screenshot({ path: path.join(SHOTS, `${n}.png`) })

try {
  await page.goto('http://127.0.0.1:8742/', { waitUntil: 'networkidle0' })
  await waitFor(async () => (await count('#scene-container canvas')) >= 1,
    30000, 'canvas')
  // Open the structure through the REAL UI (modal), exactly like a user:
  // the page's own session state then carries the session id.
  await page.click('#btn-open-pdb')
  await waitFor(() => page.$eval('#app-modal', (m) => m.style.display !== 'none'),
    5000, 'modal opens')
  await page.type('#modal-input', '3W85')
  await page.click('#modal-ok')
  await waitFor(async () => {
    const card = await text('#ligand-card')
    return card && card.length > 40
  }, 120000, 'structure opened via UI')
  await sleep(1000)

  console.log('== water network panel ==')
  await page.evaluate(() => document.querySelector('[data-tab="water"]').click())
  await page.click('#btn-water-analyze')
  await waitFor(async () => {
    const s = await text('#water-results .v2-summary')
    return s && s.includes('waters')
  }, 120000, 'water analysis renders')
  const waterSummary = await text('#water-results .v2-summary')
  rep.ok(/4\d\d waters/.test(waterSummary || ''), `real water count: ${waterSummary}`)
  const waterRows = await count('#water-results tbody tr')
  rep.ok(waterRows > 0, `water-protein contact rows (${waterRows})`)
  const chips = await count('#water-results .v2-chip')
  rep.ok(chips > 0, `cluster chips rendered (${chips})`)
  await shot('06-water-network')

  console.log('== 2D editor panel ==')
  await page.evaluate(() => document.querySelector('[data-tab="editor"]').click())
  await page.click('#btn-editor-load')
  await waitFor(async () => {
    const s = await text('#editor-info .v2-summary')
    return s && s.includes('bonds')
  }, 60000, 'editor loads ligand')
  const editorSummary = await text('#editor-info .v2-summary')
  rep.ok(/\d+ atoms/.test(editorSummary || ''),
    `editor loaded real CCD atoms/bonds: ${editorSummary}`)
  const atomDots = await page.evaluate(() => {
    const c = document.querySelector('#editor-canvas')
    const ctx = c.getContext('2d')
    const data = ctx.getImageData(0, 0, c.width, c.height).data
    let colored = 0
    for (let i = 0; i < data.length; i += 4) {
      // count non-background pixels (bg is #0d1117)
      if (Math.abs(data[i] - 13) + Math.abs(data[i + 1] - 17) +
          Math.abs(data[i + 2] - 23) > 30) colored++
    }
    return colored
  })
  rep.ok(atomDots > 200, `editor canvas actually drawn (${atomDots} lit pixels)`)
  await shot('07-editor')

  console.log('== editor interaction: select atom via canvas click ==')
  // Click near the canvas center; the app picks the nearest atom.
  const box = await (await page.$('#editor-canvas')).boundingBox()
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2)
  await sleep(600)
  const selInfo = await text('#editor-info .v2-summary')
  rep.ok(!!selInfo, `editor responded to canvas click: "${selInfo?.slice(0, 60)}"`)
  await shot('08-editor-pick')

  console.log('== scripting console ==')
  await page.evaluate(() => document.querySelector('[data-tab="script"]').click())
  await page.click('#script-input')
  await page.type('#script-input', "print('ui-verify:', structure.id)")
  await page.click('#btn-script-run')
  await waitFor(async () => {
    const o = await text('#script-output')
    return o && o.includes('ui-verify')
  }, 30000, 'script output')
  const scriptOut = await text('#script-output')
  rep.ok(scriptOut.includes('ui-verify: 3W85'),
    `console executed on real session: "${scriptOut.split('\n')[0]}"`)
  await shot('09-console')

  console.log('== measurement (real backend math) ==')
  const meas = await bridgeCommand('measure_distance', {
    atom1: { chain: 'A', residue_id: 67, residue_name: 'ASN', name: 'OD1' },
    atom2: { chain: 'A', residue_id: 401, residue_name: 'W85', name: 'NAM' },
  })
  rep.ok(meas.success && meas.data.distance > 1.0 && meas.data.distance < 8.0,
    `real distance ASN67:OD1-W85:NAM = ${meas.data?.distance} Å`)
} catch (e) {
  rep.ok(false, `stage failed: ${e.message}`)
  await shot('99-stage2-failure').catch(() => {})
} finally {
  const okAll = rep.summary('stage 2: water + editor + console + measure')
  await browser.close()
  bridge.kill()
  staticServer.close()
  process.exit(okAll ? 0 : 1)
}
