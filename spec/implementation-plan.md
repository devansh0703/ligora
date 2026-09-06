# Implementation Plan — Ligand Interface Analysis Workstation
# Working name: Ligora (provisional)
# Stack: Tauri (Rust webview shell) + Python backend (torch where useful, PLIP-backed contact analysis) + embedded WebGL molecular viewer
# Targets: V1 (shipable analysis workstation) + V2 (batch/2D-sync/water/scripting)
# Distribution: Snap (strict confinement with GPU), plus portable artifacts later.

## 1. What torch is and is not here
Torch is used as a GPU-backed numeric/learning component, NOT as a molecular dynamics or QM engine. This app is NOT a simulation package. Concretely:

### Torch IS used for
- **Ligand similarity/embeddings.** Fingerprint-based Tanimoto is fine for V1; torch can host learned molecular embeddings for nicer similarity/nearest-neighbor context if and when we have data and a model worth running locally. This is optional and pragmatic, not a selling point.
- **Optional learned confidence on contacts/pockets.** If at some point there is a small local model that improves contact/pocket confidence beyond geometry heuristics, torch is where it would live. Do not build this just to use torch. Start with geometry + PLIP, add ML only when it clearly helps and the model is small enough to bundle.
- **Numerical heavy lifting where GPU helps.** Anything that becomes a batch/array problem (similarity over a local ligand set, contact featurization, comparison math) can use torch for convenience and GPU when present. Again: only if it helps.

### Torch is NOT used for
- MD, minimization, docking, scoring, QM, dynamics, trajectory replay as a first-class feature. Those are out of scope (per spec non-goals). If people want them, they go to ChimeraX/PyMOL/Rosetta, not here.

### Practical decision
Start with PLIP-backed analysis and classic cheminformatics (fingerprints/Tanimoto) in V1. Add torch paths incrementally and only where they clearly improve the workflow and can be bundled small. Do not force a heavy DL stack into a snappable workstation that should feel local-first.

## 2. Architecture (Tauri + Python backend + embedded WebGL)
### 2.1 Host shell (Tauri, Rust)
Tauri owns:
- Window, menus, file handling, system integration, packaging/Snap.
- The webview that renders the UI.
- The bridge to the Python backend.

Pattern: Tauri spawns/manages a local Python process or uses Tauri's command bridge to call into Python via a small IPC layer (STDIO JSON, or a local socket, or Tauri invoke into a bundled Python runtime). Important properties:
- Python runs as the analysis backend; the frontend sends structured commands and gets back JSON + blobs (contact tables, SDF, images) + scene data.
- Scene geometry can be cooked by Python and streamed to the WebGL viewer, OR the frontend loads parsed structure data and the viewer renders it.
- Offline-first: everything meaningful works without network; network calls are enrichment via the Python backend.

### 2.2 Python backend (the analysis core)
Responsibilities:
- Parse/mmCIF/PDB load, ligand detection/resolution, coordinate handling.
- Contact analysis via PLIP (subprocess/ library use) plus app-level geometry as backup/extension.
- Cheminformatics: identity normalization, fingerprints, Tanimoto similarity, SDF/SMILES handling.
- Live enrichment clients: RCSB, ChEMBL, PubChem, PDBBind, with caching and honest status.
- Export: scene image (cooked or captured), contact CSV, ligand SDF/2D, analysis summary JSON.
- Optional: torch paths for similarity/embeddings/batch math.

Keep the backend as a clean library + CLI, so it can be tested headless and reused. The Tauri app is one consumer of it.

### 2.3 Frontend (Tauri webview, likely React/Vue/Solid + a WebGL molecular viewer)
- Embedded WebGL molecular viewer for the scene: NGL Viewer / 3Dmol.js / Mol* web components are the realistic choices here. Feeding parsed coordinates + simple representation directives is the main integration.
- UI panels: structure/ligand opener, ligand card, contact table, evidence pane, measurement tools, comparison, notes/annotations, export.
- All analysis state owned by the backend; frontend is the presentation + interaction layer.

This matches your choices: PLIP does the contact heavy lifting, the WebGL viewer shows the scene, torch is only where it helps.

## 3. Data-source map → concrete implementation
This is the live-data spine. Implement as backend clients with caching and honest status; network is optional enrichment.

### 3.1 Structure / coordinates
- **RCSB Search API v2**: resolve PDB ID / CSM-ID, find candidate ligands, related structures (same ligand, similar ligand, sequence/3D neighbors).
- **RCSB Data API / GraphQL**: per-entry annotations (chains, entities, ligands, chem_comp, observations).
- **RCSB ModelServer**: fetch coordinates efficiently (BinaryCIF, subsets by assembly/chain/ligand). This is the efficient path for on-demand loading.
- **RCSB File Download Services**: full mmCIF/PDB for local cache or when needed.
- **Local files**: mmCIF primary, PDB legacy fallback. App treats mmCIF as first-class.

Implementation notes:
- Prefer ModelServer for coordinate fetching in the live path; keep a local cache for repeat work.
- For offline, a local mmCIF/PDB is enough; all geometric analysis and export work from it.

### 3.2 Ligand identity / chemistry
- **wwPDB Chemical Component Dictionary (CCD)** / RCSB chem_comp + core_nonpolymer_entity: ligand identity from the entry (name, formula, etc.).
- **PubChem**: compound identity, synonyms, annotations (toxicity, solubility-type data where present) when the ligand maps.
- **ChEMBL molecule endpoint**: properties, structures (SMILES/InChI), 2D image, synonyms, xrefs when the ligand maps.
- **SDF/SMILES handling**: local cheminformatics for normalization, fingerprinting, export.

Status handling (part of the product, not a detail):
- Resolved / partial / not found — shown in the ligand card, not faked.

Mapping strategy:
- From a bound ligand in a PDB entry, start from the chemical component / residue name and formula; try to map to PubChem/ChEMBL where possible via identifiers and structure keys; accept partial mappings.

### 3.3 Evidence context
- **PDBBind**: curated binding affinities for protein-ligand complexes (affinity layer where present).
- **ChEMBL**: activity (binding/functional/ADMET-type), target, assay, binding_site; similarity/substructure search; molecule/target/activity/assay/binding_site endpoints.
- **UniProt**: target identity/function/features/disease context; available via RCSB annotations and UniProt APIs.
- **RCSB Search**: related structures as optional context.

Implementation notes:
- Enrichment is best-effort and clearly labeled; absence is fine.
- Cache what you fetch for the session and optionally beyond.

### 3.4 Owned / local computation
- Geometry/contact analysis: PLIP-backed + app geometry.
- Cheminformatics: fingerprints, Tanimoto, SDF/SMILES normalization and export.
- Torch paths: similarity/embeddings/batch math, only where they help.

## 4. Backend feature spec (V1 + V2)

### 4.1 Structure load
- Load mmCIF (primary) and PDB (legacy fallback) from local file or fetched via RCSB ModelServer/Download.
- Parse polymers and non-polymer entities; identify candidate small-molecule ligands.
- Expose: chains, residues, atoms, ligand entities, ligand atom sets.

### 4.2 Ligand detection & resolution (V1)
- Candidate ligand detection from non-polymer entities.
- Ligand card data: chemical component/residue name, formula, MW, 2D where available, SMILES/SDF where resolvable, classification hint, resolution status.
- Mapping to PubChem/ChEMBL where possible; honest status everywhere.

### 4.3 Contact analysis (V1, PLIP-backed)
- Run PLIP on the chosen structure/ligand to get interaction types and details.
- Present: contact table (chain, residue, atom, distance, contact type), pickable, with scene highlighting.
- App-level geometry as a fallback/extension and to fill gaps PLIP does not cover in a given case (e.g., custom distance/angle queries, pocket shell, contact-type heuristics where useful).
- Measurement: distance and angle for selected atom pairs.

### 4.4 Pocket (V1)
- Ligand-centric pocket shell on demand; toggle; widen/narrow option; surface coloring by contact status where relevant.

### 4.5 Live enrichment (V1, optional)
- PubChem ligand mapping/identity/annotations where available.
- ChEMBL bioactivity summary where the ligand maps.
- PDBBind affinity where the complex/ligand is covered.
- Related structures via RCSB Search where useful.
- Clear offline state when none of this is available.

### 4.6 Export (V1)
- Scene image: viewport screenshot or cooked render.
- Contact table CSV.
- Ligand SDF where resolvable; 2D image where available.
- Analysis summary JSON/text: ligand identity, contact summary, evidence summary, notes.
- Notes/annotations preserved per analysis; reloadable in V1.

### 4.7 Comparison (V1 light, V2 fuller)
- V1: open a second structure/ligand, side-by-side/overlay basics.
- V2: same-ligand across structures, pocket conservation summary, ligand similarity context via fingerprints/Tanimoto and/or learned embeddings.

### 4.8 V2 additions
- **Batch analysis**: analyze a list of PDB IDs / local set, produce a comparative summary (med-chem/series use case).
- **2D ligand edit + 3D sync**: edit ligand in 2D, reflect in 3D scene — only if it serves the analysis workflow, not a general editor.
- **Water network**: water-mediated contact analysis where the structure supports it.
- **Scripting/console**: small focused console for power users, below the main workflow.
- **Torch paths mature**: similarity/embeddings/batch math where useful; learned contact/pocket confidence only if a small bundleable model justifies it.

## 5. Frontend feature spec (V1 + V2)

### 5.1 UI skeleton (V1)
- Top: open (PDB ID / local file), search, recent, settings.
- Scene area: embedded WebGL viewer, viewer controls, pick/hover atom info.
- Side/dock: ligand card, contact table, evidence pane, measurement tools, export.
- Notes panel: free notes + structured annotations attached to the analysis.

### 5.2 Scene interaction (V1)
- Representations directed from app state: protein cartoon/ribbon by chain, ligand ball-and-stick with atom labels on pick/hover, stick for near-ligand sidechains, ligand-centric surface on demand.
- Coloring: by chain; by contact status; ligand-centric palette.
- Picking: atom/residue pick drives ligand card + contact table + scene highlight.
- Camera: rotate/zoom/pan; a few preset views for the interface/pocket/ligand.

### 5.3 Analysis interaction (V1)
- Contact table with pick-to-highlight, sort/filter basics.
- Distance/angle measurement with readout.
- Pocket surface toggle.
- Ligand card with honest resolution status + enrichment where available.

### 5.4 Export & artifact (V1)
- Export buttons: scene image, contact CSV, ligand SDF/2D, analysis summary.
- Save/load analysis artifact (state + notes) for reload in V1.

### 5.5 V2 UI additions
- Comparison workspace (side-by-side/overlay, two analyses).
- Batch results view (table + per-entry drill-in).
- 2D editor panel with 3D sync toggle.
- Water network layer toggle.
- Console panel for scripting.

## 6. Implementation order (so the product is real early)
### Phase 0 — scaffold
- Repo with Tauri app + Python backend package + shared types/IPC.
- Packaging layout for Snap and portable artifacts.
- Basic Tauri window + a placeholder scene + a "load local mmCIF" path that reaches the backend and shows parsed chains/ligands.

### Phase 1 — core offline analysis (MVP within V1)
- mmCIF/PDB load + ligand detection.
- PLIP-backed contact analysis + contact table + pick-to-highlight.
- Measurement: distance/angle.
- Ligand card with identity from archive data; honest status.
- Pocket surface toggle.
- Export: scene image + contact CSV + ligand SDF/2D where available + analysis summary.
- Save/load analysis artifact.

### Phase 2 — live enrichment (V1 enrichment layer)
- RCSB resolve/fetch (Search + Data/GraphQL + ModelServer/Download) with cache.
- PubChem/ChEMBL/PDBBind enrichment clients with honest status and offline fallback.
- Evidence pane in the UI.

### Phase 3 — polish + comparison (V1 completion)
- Representation/color controls, camera presets, picking ergonomics.
- Simple comparison (open second structure/ligand).
- Render/export quality improvements.

### Phase 4 — V2
- Batch analysis.
- 2D ligand edit + 3D sync.
- Water network.
- Scripting console.
- Torch paths matured only where they clearly help and stay small enough to bundle.

## 7. Packaging & Snap (strict confinement, GPU)
### 7.1 Confinement model
Strict confinement with:
- **opengl** plug for GPU rendering (both the WebGL viewer's GL context and any GPU work in the backend).
- Desktop plugs as needed by Tauri/webview (display, x11/wayland as appropriate, session integration).
- **home** or file-selector portal for user structure files.
- **network** for live enrichment (PubChem/ChEMBL/RCSB/PDBBind), optional, graceful offline.
- Optionally **session-lock** later if a presentation/kiosk mode is added.

### 7.2 Bundle shape
- Tauri app bundle with the webview frontend.
- Python backend bundled for the target (so the snap is self-contained). Options: ship a bundled Python runtime with the backend packaged inside the snap, or bundle a venv/conda-like tree; choose based on size and build-time constraints. The point: the snap should not depend on the host having a particular Python/chemlib stack installed.
- PLIP and cheminformatics deps available inside the bundle.
- torch present if/when a torch path is actually used; if a torch model is bundled, keep it small and optional.

### 7.3 GPU in strict snap
- opengl plug + desktop helpers as Tauri/webview needs; validate on X11 and Wayland early.
- If the WebGL viewer is the renderer, its GL context still needs opengl to work in the snap; verify early with a real GPU setup.
- Keep a software/worst-case path for scene display if GPU setup fails (even if slower), so the app still opens.

### 7.4 Channels & updates
- Snap Store publish with a channel strategy (e.g., edge/candidate/stable).
- Offline-first means auto-update is fine; the app still works if a refresh is delayed.

### 7.5 Secondary artifacts
- After Linux/Snap ships, add a portable tar/AppImage and consider Windows/macOS as "nice to have," not MVP blockers.

## 8. What "no similar app" means in implementation terms
The implementation must preserve the combined workflow, because that is the gap:
- Scene + PLIP-backed contact analysis + ligand identity resolution + open-live evidence context + structured export, all local-first with network as optional enrichment, in one snappable app.
- If we ever drift to "viewer with a PLIP button" or "data browser with a scene," we have lost the product. The integration and the artifact are the product.

## 9. Risks and mitigations
- **PLIP integration + licensing/distribution in a Snap.** Validate that PLIP can be bundled/run in the snap environment and that its dependencies fit; if PLIP is awkward to bundle, treat app geometry as the primary contact path and PLIP as enrichment/fallback.
- **WebGL viewer integration with Tauri.** Pick the viewer early and prove the feed path (coordinates + representation directives) before building UI; some viewers want specific input shapes.
- **mmCIF-first correctness.** Build on mmCIF with PDB legacy fallback; use RCSB ModelServer/BinaryCIF where possible instead of hand-rolling everything.
- **Ligand identity mapping is patchy.** Design the ligand card for honest statuses from the start.
- **Torch bloat.** Only add torch paths when they clearly help and stay small; do not let "torch for GPU" become a forced heavy DL stack.
- **Live-data dependence.** Core analysis must be fully usable offline from a local file; network enrichment is bonus.

## 10. Acceptance checks (V1)
- Open a representative bound PDB entry and a bound CSM entry; pick a ligand; see ligand-centered analysis (scene + ligand card + contact table + optional evidence where network allows).
- Contact table pickable and correct; distances/angles correct; contact types transparent.
- Ligand identity honest about resolved/partial/not-found; PubChem/ChEMBL/PDBBind enrichment appears where it exists.
- Export produces a coherent artifact (scene image + contact CSV + ligand SDF/2D where available + summary/notes), reloadable in V1.
- Works offline from a local file; network enrichment degrades gracefully.
- Snap installs, launches from app grid + CLI, renders with GPU on default Ubuntu (X11 + Wayland where feasible), strict confinement holds for the chosen plugs.
- No existing app ships this combined workflow as one snappable product; implementation preserves that.

## 11. Next concrete step
Scaffold the repo with the Tauri + Python backend split and the IPC/shared-types boundary, then implement Phase 1 (core offline analysis) end to end: load local mmCIF, detect ligand, run PLIP-backed contact analysis, show contact table + scene pick sync, measure, export. That proves the product spine before any live enrichment or V2.
