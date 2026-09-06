/**
 * Shared harness for staged UI verification of Ligora.
 * Real app + real backend; no mocks. See ui-01*.mjs ... ui-04*.mjs.
 */
import puppeteer from 'puppeteer-core'
import { spawn } from 'node:child_process'
import { createServer } from 'node:http'
import { readFile } from 'node:fs/promises'
import { fileURLToPath } from 'node:url'
import path from 'node:path'

const ROOT = path.dirname(path.dirname(fileURLToPath(import.meta.url)))
export const DIST = path.join(ROOT, 'dist')
export const SHOTS = path.join(ROOT, 'tests', 'ui-screenshots')
export const BACKEND_ROOT = path.resolve(ROOT, '..', 'backend')

const mime = {
  '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css',
  '.svg': 'image/svg+xml', '.png': 'image/png', '.ico': 'image/x-icon',
}

export async function startStatic(port = 8742) {
  const server = createServer(async (req, res) => {
    try {
      const url = new URL(req.url, 'http://localhost')
      const file = path.join(DIST, url.pathname === '/' ? 'index.html' : url.pathname)
      const content = await readFile(file)
      res.writeHead(200, { 'Content-Type': mime[path.extname(file)] || 'application/octet-stream' })
      res.end(content)
    } catch {
      res.writeHead(404).end()
    }
  })
  await new Promise((r) => server.listen(port, '127.0.0.1', r))
  return server
}

export function startBridge(port = 8741) {
  const bridge = spawn('node', [path.join(ROOT, 'dev-bridge.mjs')], {
    env: { ...process.env, PYTHONPATH: BACKEND_ROOT, LIGORA_BRIDGE_PORT: String(port) },
    stdio: ['ignore', 'pipe', 'inherit'],
  })
  const ready = new Promise((resolve, reject) => {
    const t = setTimeout(() => reject(new Error('bridge start timeout')), 25000)
    bridge.stdout.on('data', (d) => {
      if (String(d).includes('listening')) { clearTimeout(t); resolve() }
    })
    bridge.on('exit', (code) => reject(new Error(`bridge exited early (${code})`)))
  })
  return { bridge, ready }
}

export async function bridgeCommand(type, payload = {}, port = 8741) {
  const res = await fetch(`http://127.0.0.1:${port}/command`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ type, payload }),
  })
  return res.json()
}

export async function launch() {
  const browser = await puppeteer.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: 'new',
    args: [
      '--no-sandbox', '--disable-dev-shm-usage',
      '--use-gl=swiftshader', '--enable-unsafe-swiftshader',
      '--window-size=1440,960',
    ],
    defaultViewport: { width: 1440, height: 960 },
  })
  const page = await browser.newPage()
  page.setDefaultTimeout(30000)
  page.on('pageerror', (e) => console.log('  [pageerror]', String(e).slice(0, 200)))
  await page.evaluateOnNewDocument(() => {
    window.__LIGORA_BRIDGE__ = async (type, payload) => {
      const res = await fetch('http://127.0.0.1:8741/command', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ type, payload }),
      })
      return res.json()
    }
  })
  return { browser, page }
}

export function makeReport() {
  let passed = 0
  let failed = 0
  const failures = []
  return {
    ok(cond, label) {
      if (cond) { passed++; console.log(`  ✔ ${label}`) }
      else { failed++; failures.push(label); console.log(`  ✘ ${label}`) }
    },
    summary(name) {
      console.log(`\n[${name}] ${passed} passed, ${failed} failed`)
      for (const f of failures) console.log('  ✘', f)
      return failed === 0
    },
  }
}

export const sleep = (ms) => new Promise((r) => setTimeout(r, ms))

export async function waitFor(fn, timeoutMs = 60000, label = 'condition') {
  const t0 = Date.now()
  while (Date.now() - t0 < timeoutMs) {
    try { if (await fn()) return true } catch { /* retry */ }
    await sleep(400)
  }
  throw new Error(`timeout waiting for ${label}`)
}
