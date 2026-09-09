# Ligora — Working UI Screenshots

Real screenshots of the packaged Ligora app (core26 snap) driving **live data** —
no mocks. Captured from the actual WebKit webview while running the real
Python backend (PLIP, Vina, RDKit) and live data sources (RCSB, PubChem, ChEMBL).

| Screenshot | What it shows |
|---|---|
| `01-app-start.png` | App on first launch — full workstation layout |
| `02-structure-3w85-loaded.png` | **3W85** fetched from RCSB and rendered by 3Dmol.js (WebGL) |
| `03-contacts-plip-3w85.png` | Real **PLIP contact analysis**: hydrogen bonds ASN67/ASN194/ASN127, distances in Å |
| `04-evidence-pane-3w85.png` | **Evidence pane** with live provenance: RCSB CCD, PubChem, UniChem endpoints |
| `05-search-myoglobin.png` | **PDB-wide search** for "myoglobin" — 15 live hits |
| `06-docking-panel-3w85.png` | Docking panel (Vina) with box controls, exhaustiveness, MD section |
| `07-docking-results-3w85.png` | **Real AutoDock Vina run**: 5 poses, affinity table (kcal/mol) |

## Reproduce

```bash
snap run ligora   # or: cd frontend && npm run tauri:dev
```

1. Enter a PDB ID (e.g. `3W85`) → *Open*
2. *Run Analysis* → contacts table populates from PLIP
3. *Docking* tab → *Run Docking* → Vina poses with affinities
