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
        if (structure.file_path) {
          const result = await invoke('open_local_file', { filePath: structure.file_path })
          if (!result || !result.success) {
            throw new Error(result && result.error ? result.error : 'open local file failed')
          }
          return result.data
        }
        if (structure.pdb_id) {
          const result = await invoke('open_pdb_id', { pdbId: structure.pdb_id })
          if (!result || !result.success) {
            throw new Error(result && result.error ? result.error : 'open PDB ID failed')
          }
          return result.data
        }
        throw new Error('structure must include file_path or pdb_id')
      },
      onLigandSelected: async (ligandId) => {
        const result = await invoke('select_ligand', { ligandId })
        if (!result || !result.success) {
          throw new Error(result && result.error ? result.error : 'ligand selection failed')
        }
        return result.data || {}
      },
      onRunContactAnalysis: async () => {
        const result = await invoke('run_contact_analysis')
        if (!result || !result.success) {
          throw new Error(result && result.error ? result.error : 'contact analysis failed')
        }
        const data = result.data || {}
        if (data.ligand_resolved) {
          data.ligand = data.ligand_resolved
        }
        return data
      },
      onRunDocking: async (params) => {
        const result = await invoke('run_docking', { params })
        if (!result || !result.success) {
          throw new Error(result && result.error ? result.error : 'docking failed')
        }
        return result.data
      },
      onGetStatus: async () => {
        const result = await invoke('get_status')
        if (!result || !result.success) {
          throw new Error(result && result.error ? result.error : 'status fetch failed')
        }
        return result.data
      },
      onOpenPdbId: async (pdbId) => {
        const result = await invoke('open_pdb_id', { pdbId })
        if (!result || !result.success) {
          throw new Error(result && result.error ? result.error : 'open PDB ID failed')
        }
        return result.data
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
