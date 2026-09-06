#!/usr/bin/env node
/**
 * Dev-only bridge: runs the REAL Python backend (ligora_backend.server)
 * and exposes its newline-delimited JSON IPC over a local HTTP endpoint
 * so browser automation can drive the same commands the Tauri app uses.
 *
 * This is a development/test harness only - the shipped app uses the Rust
 * IPC bridge in src-tauri/src/lib.rs. No mock data anywhere: every command
 * is executed by the real backend against live sources.
 */
import { spawnSync, spawn } from 'node:child_process'
import http from 'node:http'
import fs from 'node:fs'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

// The Python package's PARENT directory: `python -m ligora_backend.server`
// must resolve the package from here (cwd and PYTHONPATH both).
const BACKEND_ROOT = path.resolve(
  path.dirname(fileURLToPath(import.meta.url)), '..', 'backend')
const PORT = Number(process.env.LIGORA_BRIDGE_PORT || 8741)

/** Resolve a usable Python interpreter (absolute path; no PATH games). */
function findPython() {
  if (process.env.LIGORA_PYTHON && fs.existsSync(process.env.LIGORA_PYTHON)) {
    return process.env.LIGORA_PYTHON
  }
  for (const candidate of ['/usr/bin/python3', '/usr/bin/python', '/usr/local/bin/python3']) {
    if (fs.existsSync(candidate)) return candidate
  }
  const probe = spawnSync('python3', ['--version'])
  return probe.error ? null : 'python3'
}

const pythonBin = findPython()
if (!pythonBin) {
  console.error('No Python interpreter found; cannot start backend')
  process.exit(1)
}

if (!fs.existsSync(path.join(BACKEND_ROOT, 'ligora_backend', 'server.py'))) {
  console.error(`backend package not found at ${BACKEND_ROOT}/ligora_backend`)
  process.exit(1)
}

const child = spawn(pythonBin, ['-u', '-m', 'ligora_backend.server'], {
  cwd: BACKEND_ROOT,
  env: { ...process.env, PYTHONPATH: BACKEND_ROOT },
  stdio: ['pipe', 'pipe', 'inherit'],
})
child.on('error', (e) => {
  console.error(`backend spawn failed: ${e.message}`)
  process.exit(1)
})

let sessionId = ''
let buffer = ''
let seq = 0
const pending = new Map()

child.stdout.on('data', (chunk) => {
  buffer += chunk.toString()
  let idx
  while ((idx = buffer.indexOf('\n')) >= 0) {
    const line = buffer.slice(0, idx).trim()
    buffer = buffer.slice(idx + 1)
    if (!line) continue
    let msg
    try {
      msg = JSON.parse(line)
    } catch {
      continue
    }
    if (msg.event === 'job') continue // events surfaced via /events
    if (msg.command_id && pending.has(msg.command_id)) {
      pending.get(msg.command_id)(msg)
      pending.delete(msg.command_id)
    }
  }
})

function sendCommand(type, payload) {
  const command_id = `ui-${++seq}`
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => {
      pending.delete(command_id)
      reject(new Error(`backend timeout: ${type}`))
    }, 280_000)
    pending.set(command_id, (msg) => {
      clearTimeout(timer)
      if (!msg.success) reject(new Error(msg.error || `${type} failed`))
      else resolve(msg.data || {})
    })
    child.stdin.write(
      JSON.stringify({ command_id, type, payload, session_id: sessionId }) + '\n',
    )
  }).then((data) => {
    if (data?.session_id) sessionId = data.session_id
    return data
  })
}

const server = http.createServer((req, res) => {
  res.setHeader('Access-Control-Allow-Origin', '*')
  res.setHeader('Access-Control-Allow-Headers', 'Content-Type')
  if (req.method === 'OPTIONS') {
    res.writeHead(204).end()
    return
  }
  if (req.method === 'GET' && req.url === '/health') {
    res.writeHead(200, { 'Content-Type': 'application/json' })
    res.end(JSON.stringify({ ok: child.pid > 0, pid: child.pid }))
    return
  }
  if (req.method === 'POST' && req.url === '/command') {
    let body = ''
    req.on('data', (c) => (body += c))
    req.on('end', () => {
      let parsed
      try {
        parsed = JSON.parse(body || '{}')
      } catch {
        res.writeHead(400).end(JSON.stringify({ error: 'invalid JSON' }))
        return
      }
      sendCommand(parsed.type, parsed.payload || {})
        .then((data) => {
          res.writeHead(200, { 'Content-Type': 'application/json' })
          res.end(JSON.stringify({ success: true, data }))
        })
        .catch((e) => {
          res.writeHead(200, { 'Content-Type': 'application/json' })
          res.end(JSON.stringify({ success: false, error: String(e.message || e) }))
        })
    })
    return
  }
  res.writeHead(404).end()
})

server.listen(PORT, '127.0.0.1', () => {
  console.log(`ligora dev bridge listening on http://127.0.0.1:${PORT}`)
})

function shutdown() {
  try { child.stdin.end() } catch {}
  try { child.kill() } catch {}
  server.close(() => process.exit(0))
  setTimeout(() => process.exit(0), 1500)
}
process.on('SIGINT', shutdown)
process.on('SIGTERM', shutdown)
