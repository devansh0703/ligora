/**
 * Ligora - Molecular Analysis Workstation
 * Main application entry point
 */

import { invoke } from '@tauri-apps/api/core'
import { open } from '@tauri-apps/plugin-dialog'
import { readTextFile, writeTextFile } from '@tauri-apps/plugin-fs'
import { createApp } from './app'

// Initialize the app
async function main() {
  try {
    // Create and mount the app
    const app = createApp({
      onStructureLoaded: async (structure) => {
        // Send to backend
        await invoke('load_structure', { structure })
      },
      onLigandSelected: async (ligandId) => {
        await invoke('select_ligand', { ligandId })
      },
      onRunContactAnalysis: async () => {
        await invoke('run_contact_analysis')
      },
      onRunDocking: async (params) => {
        await invoke('run_docking', { params })
      },
      onSaveArtifact: async (path) => {
        await writeTextFile(path, JSON.stringify({
          type: 'artifact',
          timestamp: new Date().toISOString(),
        }))
      },
    })

    // Mount to DOM
    const container = document.getElementById('app')
    if (container) {
      container.appendChild(app)
    }

    // Start the app
    app.start()

  } catch (error) {
    console.error('Failed to initialize Ligora:', error)
    document.body.innerHTML = `
      <div style="padding: 20px; color: red;">
        <h2>Failed to initialize Ligora</h2>
        <p>${error.message || error}</p>
      </div>
    `
  }
}

// Run when DOM is ready
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', main)
} else {
  main()
}
