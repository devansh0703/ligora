/**
 * Ligora - Molecular Analysis Workstation
 * Entry point: wires the UI (app.js) to the Python backend.
 *
 * Two transports, one backend:
 * - Tauri app: commands go through the Rust bridge (`send_command` IPC),
 *   which spawns the Python backend and correlates responses by id.
 * - Browser (dev harness): `window.__LIGORA_BRIDGE__` is injected by the
 *   automation environment and posts to the dev-bridge HTTP endpoint,
 *   which drives the same real Python backend over its stdio IPC.
 */

import { invoke } from '@tauri-apps/api/core'
import { createApp } from './app'

// Session state shared by all commands.
const sessionState = { id: null }

const IN_TAURI = typeof window !== 'undefined' && '__TAURI_INTERNALS__' in window

async function tauriSendCommand(type, payload = {}) {
  const response = await invoke('send_command', {
    commandType: type,
    payload: payload,
    sessionId: sessionState.id || '',
  })
  if (!response || response.success !== true) {
    throw new Error(
      (response && response.error) || `backend command failed: ${type}`)
  }
  const data = response.data || {}
  if (data.session_id) {
    sessionState.id = data.session_id
  }
  return data
}

async function bridgeSendCommand(type, payload = {}) {
  const bridge = window.__LIGORA_BRIDGE__
  if (!bridge) throw new Error('dev bridge not available')
  const response = await bridge(type, payload)
  if (!response.success) throw new Error(response.error || `${type} failed`)
  const data = response.data || {}
  if (data.session_id) sessionState.id = data.session_id
  return data
}

const sendCommand = IN_TAURI ? tauriSendCommand : bridgeSendCommand

async function main() {
  const app = createApp({ sendCommand, sessionState })

  const container = document.getElementById('app')
  if (container) {
    // index.html already contains the full static layout; app.js upgrades
    // it in place (viewer + panel bindings) instead of duplicating it.
    app.upgrade(container)
  }

  await app.start()
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', main)
} else {
  main()
}
