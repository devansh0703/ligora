import { defineConfig } from 'vite'

export default defineConfig({
  // Tauri expects a fixed dev port (see tauri.conf.json build.devUrl).
  server: {
    port: 1420,
    strictPort: true,
  },
  build: {
    outDir: 'dist',
    emptyOutDir: true,
    target: 'es2021',
  },
})
