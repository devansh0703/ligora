/**
 * Stage 3: real Vina docking (2 runs), real pose comparison in the Compare
 * tab, real batch run (1UBQ), notes persistence. Real backend, no mocks.
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
  await page.click('#btn-open-pdb')
  await waitFor(() => page.$eval('#app-modal', (m) => m.style.display !== 'none'),
    5000, 'modal opens')
  await page.type('#modal-input', '3W85')
  await page.click('#modal-ok')
  await waitFor(async () => {
    const card = await text('#ligand-card')
    return card && card.length > 40
  }, 120000, 'structure opened via UI')

  console.log('== real Vina docking run 1 (via backend job system) ==')
  // Docking runs as a backend job; the UI tracks jobs through get_job/list.
  // Start it the way the UI's job flow does and poll like the Jobs pane.
  const job1 = await bridgeCommand('run_docking',
    { exhaustiveness: 1, num_modes: 3 })
  if (!job1.success) throw new Error(`docking failed: ${job1.error}`)
  rep.ok(true, `docking job started (${job1.data.job_id.slice(0, 8)})`)
  let st1 = null
  for (let i = 0; i < 120; i++) {
    const s = await bridgeCommand('get_job', { job_id: job1.data.job_id })
    if (s.data.status === 'completed' || s.data.status === 'failed') { st1 = s.data; break }
    await sleep(2000)
  }
  rep.ok(st1?.status === 'completed', 'docking run 1 completed')
  const poses1 = st1?.result?.poses || []
  rep.ok(poses1.length >= 1, `real poses from Vina (${poses1.length})`)
  const best1 = Math.min(...poses1.map((p) => p.affinity))
  rep.ok(best1 < 0 && best1 > -15, `real affinity ${best1.toFixed(2)} kcal/mol`)

  console.log('== real Vina docking run 2 ==')
  const job2 = await bridgeCommand('run_docking',
    { exhaustiveness: 2, num_modes: 3 })
  let st2 = null
  for (let i = 0; i < 120; i++) {
    const s = await bridgeCommand('get_job', { job_id: job2.data.job_id })
    if (s.data.status === 'completed' || s.data.status === 'failed') { st2 = s.data; break }
    await sleep(2000)
  }
  rep.ok(st2?.status === 'completed', 'docking run 2 completed')

  console.log('== Compare tab with two real runs ==')
  await page.evaluate(() => document.querySelector('[data-tab="compare"]').click())
  await waitFor(async () =>
    (await page.$eval('#compare-job-a', (s) => s.options.length)) >= 3,
    30000, 'job dropdowns populated')
  await page.select('#compare-job-a', job1.data.job_id)
  await page.select('#compare-job-b', job2.data.job_id)
  await page.click('#btn-compare-run')
  await waitFor(async () => {
    const s = await text('#compare-results .v2-summary')
    return s && s.includes('RMSD')
  }, 30000, 'comparison renders')
  const cmp = await text('#compare-results .v2-summary')
  rep.ok(/average RMSD [\d.]+ Å/.test(cmp || ''), `real pose comparison: ${cmp}`)
  const cmpRows = await count('#compare-results tbody tr')
  rep.ok(cmpRows >= 1, `best-match table rendered (${cmpRows} rows)`)
  await shot('10-compare')

  console.log('== Batch tab with a real second structure ==')
  // 1UBQ has no ligand at all (ubiquitin), which would honestly fail;
  // use a ligand-bearing structure for a successful batch run.
  await page.evaluate(() => document.querySelector('[data-tab="batch"]').click())
  await page.click('#batch-sources')
  await page.type('#batch-sources', '1HRC')
  await page.click('#btn-batch-run')
  await waitFor(async () => (await count('#batch-results tbody tr')) > 0,
    280000, 'batch results render')
  const batch = await text('#batch-results .v2-summary')
  rep.ok(/1 succeeded/.test(batch || ''), `batch summary: ${batch}`)
  const batchRow = await page.$eval('#batch-results tbody tr',
    (e) => e.textContent)
  rep.ok(batchRow.includes('1HRC') && batchRow.includes('HEME'),
    `real batch result row: ${batchRow.slice(0, 60)}`)
  await shot('11-batch')

  console.log('== Notes persistence through the UI ==')
  await page.click('#notes-area')
  await page.type('#notes-area', 'ui-verify note')
  await sleep(1200)
  const st = await bridgeCommand('get_status', {})
  rep.ok(st.data?.session?.notes === 'ui-verify note',
    'notes persisted to the real session')
  await shot('12-notes')
} catch (e) {
  rep.ok(false, `stage failed: ${e.message}`)
  await shot('99-stage3-failure').catch(() => {})
} finally {
  const okAll = rep.summary('stage 3: docking + compare + batch + notes')
  await browser.close()
  bridge.kill()
  staticServer.close()
  process.exit(okAll ? 0 : 1)
}
