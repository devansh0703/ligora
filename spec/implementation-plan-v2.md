# Implementation Plan — Ligand Interface Analysis Workstation
# Working name: Ligora (provisional)
# Stack: Tauri (Rust webview shell) + Python backend (torch where useful, PLIP-backed contact analysis) +
#        embedded WebGL molecular viewer + orchestration of real external engines (docking, MD, QM, etc.)
# Distribution: Snap (strict confinement for the app core; engines handled with an explicit strategy) +
#               portable artifacts later.

## 1. What "real simulation engine" means here
This app does not reimplement docking/MD/QM. It integrates existing, proven local engines and makes them part of one coherent analysis workflow in one snappable app. The engines are the muscle; this app is the workstation, session manager, file/parameter plumbing, scene/visual context, and artifact assembly.

Realistic local engines to integrate (open or free-for-the-use-case), with the important caveats:

### Docking
- **AutoDock Vina** (open source, CLI, widely used, fast). Good baseline.
- **smina / gnina** (Vina-family; gnina adds CNN scoring and more options). Good for stronger scoring/pose prediction.
- These are CLI tools; the app wraps them, manages inputs/outputs, and visualizes results in-scene.

### Molecular dynamics
- **OpenMM** (open toolkit, Python API, GPU-capable). Best fit for "local MD integrated into a Python backend" because it has a Python API and GPU support; can be run as a library, not only a CLI.
- **GROMACS** (open, widely used, CLI, GPU-capable). Excellent but heavier to bundle; commonly user-provided or installed separately.
- **NAMD** (free but not fully OSS-ish depending on context/license). Be careful; treat as "user-provided" unless you verify redistribution terms.
- Production MD is a heavier ask than a snappable workstation usually wants to bundle. Treat MD as an engine path, with the realistic bundling model below.

### Quantum chemistry / electronic structure
- **PySCF** (open source, Python, Apache 2.0) — best fit for integration because it's a Python library, easy to call from the backend, and redistributable under a clean license.
- **Psi4** (open source, Python, GPL-ish — verify current license before bundling).
- **ORCA** (free for academic/personal use, but redistribution terms matter; do not assume it can be bundled in a freely distributed Snap without checking the license for the exact distribution model).
- QM is the heaviest and most license-sensitive; integrate carefully and prefer PySCF for the bundled path unless a user brings ORCA/Psi4 themselves.

### What is NOT in scope as a first-class engine integration
- Building/customizing force fields from scratch.
- Long production MD campaign management as a replacement for specialized MD tooling.
- Proprietary commercial suites (Schrödinger, MOE, etc.) — not open, not bundleable here.
- Anything that requires redistribution of a non-bundleable binary in a freely distributed Snap without checking terms.

## 2. Big packaging/licensing reality check (this is the product decision)
A "full V1+V2 with real engines" and a "clean strict-Snap that bundles all engines" are, in practice, in tension. This must be decided explicitly.

### 2.1 Snap strict confinement basics
- The app core can be strict-confined with opengl + desktop + home + network.
- At runtime, a strict snap generally **cannot just `pip install` arbitrary things from the internet into a system location or writable system Python** in the normal way; it has network only per its plugs, and the bundle is what you ship. You can have the app download/install things at runtime into its own writable area only if that is intentional and within confinement, but that is not the same as "bundle everything."
- So the clean model is: **bundle the app core + a chosen minimal engine set**, and treat heavier or license-sensitive engines as "user-provided / user-installed / classic-snap / external runtime" with a clean integration path.

### 2.2 Which engines can realistically be bundled in a strict Snap
- **Easy to bundle (Python-library engines):** PySCF (Apache 2.0), and any pure-Python/numerical backend that installs from PyPI and is license-clean. These are good bundled candidates.
- **Bundlable but needs care (CLI tools you can ship):** Vina/smina/gnina if licensing and binary redistribution are fine and you can ship the binaries for the target arches. Validate each.
- **Hard to bundle cleanly in a strict free Snap:** heavy MD binaries (GROMACS fully bundled for all arches is a lot), and especially license-sensitive binaries like ORCA if the redistribution terms do not permit it for your distribution model. For those, the realistic model is "user provides / user installs, app integrates via a known interface."

### 2.3 Two packaging strategies, pick based on your goals
#### Strategy A — strict Snap app core + bundled small engine set + "bring your own engine" for the rest
- Bundle: app core, PLIP-backed analysis, cheminformatics, PySCF (if used), and possibly Vina/gnina if redistributable and you ship binaries for amd64/arm64.
- MD/QM heavier engines: integrate by discovering them on the user's system (PATH, known install locations, or a user-configured engine path), or via a separate classic-snap / separate runtime the user installs.
- This keeps the strict-Snap story clean and the app distributable, while still being a real workstation that runs engines the user has.

#### Strategy B — classic confinement Snap for a "batteries-included workstation"
- Classic confinement can make it easier to bundle more engines and to install/enrich at runtime, at the cost of the strict-confinement story.
- If "one snap with everything inside" is more important than strict confinement, this is the tradeoff. It also raises more review/questions for the Snap Store.

### 2.4 Recommendation (to keep the product coherent)
- Design the app core as strict-confinement-ready, with engine integration as a first-class but pluggable layer.
- Bundle the small clean engine set you can legitimately redistribute.
- Support "user-provided engines" as a supported integration path, with a clean UI for engine discovery/selection and clear messaging when an engine is missing.
- Be explicit and honest in the product about which engines are bundled vs user-provided. That honesty is part of the product's credibility.

## 3. Product model with engines (what the user actually gets)
The app is still the analysis workstation from before, plus an **engine orchestration layer** that lets the user do real computational work on the same structure/ligand they are analyzing, with results flowing back into the scene and the artifact.

### Core workflow with engines
1. Open a structure / pick a ligand / analyze the interface (scene + contacts + ligand identity + evidence).
2. Choose an engine task relevant to the ligand/target/pocket:
   - Dock a ligand (or a set of ligands) into the pocket.
   - Run a short local MD / minimization to relax geometry or probe stability (where appropriate and where the user has the engine).
   - Run a quantum/cheap electronic-structure computation on the ligand or a fragment for properties/charges/energy where the user has the engine.
3. The app prepares inputs from the app state (structure, ligand, pocket selection, parameters), invokes the engine, monitors the run, parses outputs, and presents results in-scene and in the artifact.
4. Results become part of the analysis: pose(s) visualized, contact changes after relaxation visible, energies/properties shown where relevant, and everything exportable into the artifact.

This is the "deep analysis" with real engines — not just looking, but doing real work in one place.

## 4. Architecture (Tauri + Python backend + embedded WebGL + engine adapters)
### 4.1 Host shell (Tauri, Rust)
Same as before: window, menus, file handling, system integration, Snap/packaging, webview, IPC to Python backend. Adds: engine status/launch UI wiring, job monitoring, and result handling via the backend.

### 4.2 Python backend (analysis + orchestration)
Responsibilities:
- Analysis spine: structure load (mmCIF/PDB), ligand detection/resolution, contact analysis (PLIP + app geometry), cheminformatics, live enrichment, export.
- Engine orchestration layer: adapters per engine, input preparation, job launch/monitoring, output parsing, result normalization into a common result model.
- Session/job management: queue, status, logs, cancel, result lifecycle.
- Optional torch paths: similarity/embeddings/batch math, and possibly engine-result featurization if it helps.

Keep engine adapters as a pluggable subsystem. Common result model is the key: different engines produce different outputs; the app normalizes what it needs (pose, contacts, energy/properties, logs) into one shape the UI and artifact can use.

### 4.3 Frontend (Tauri webview + embedded WebGL viewer)
- Scene + ligand card + contact table + evidence pane + measurement + notes + export, as before.
- Adds: engine panel (engine selection, parameters, job control, logs, results list), result visualization in-scene (poses, post-relaxation contacts, etc.), result diff/compare where useful.

## 5. Backend feature spec (V1 + V2)

### 5.1 Analysis spine (V1, unchanged core)
- Structure load (mmCIF primary, PDB legacy fallback), local or fetched via RCSB.
- Ligand detection/resolution + ligand card with honest status.
- Contact analysis via PLIP + app geometry; contact table + pick sync + measurement.
- Pocket shell option.
- Live enrichment: RCSB, PubChem, ChEMBL, PDBBind, UniProt-adjacent context; honest offline fallback.
- Export: scene image + contact CSV + ligand SDF/2D where available + analysis summary; save/load analysis artifact.

### 5.2 Engine orchestration (V1 — start small and real)
Pick a minimal real engine set for V1 that is both useful and realistically bundleable/integrable:
- **Docking (V1):** at least one CLI docking adapter (e.g., Vina or gnina) wired in; adapter model is pluggable so more can be added.
  - Input prep: protein + ligand + pocket/box definition derived from the app state (ligand-centric pocket or user-defined box).
  - Run + monitor + parse outputs into poses + scores.
  - Visualize poses in-scene; compare poses; include in artifact.
- **Local geometry cleanup / quick minimization (V1-lite):** a small local step (backend geometry or a simple engine call) to check ligand strain / relax the picked ligand's geometry for analysis context — not production MD. This is the "simulation-lite" that belongs in V1.

Keep V1 engine scope tight: one docking adapter + a small geometry/relaxation helper. That is already unusually capable for a single snappable workstation and is real. Everything bigger is V2.

### 5.3 Engine adapters (plugable, V1 foundation, V2 expansion)
Common adapter interface:
- Discover/engine path (bundled, or user-provided via PATH/config).
- Input preparation from app state.
- Run + status + logs + cancel.
- Output parsing into the common result model.
- Result visualization hooks.

Engines to add in V2 (as bundling/user-provided strategy allows):
- **Docking:** add gnina and/or smina as additional adapters; support scoring comparisons.
- **MD:** OpenMM as a library adapter (Python, GPU-capable) is the best V2 fit; GROMACS as a CLI adapter where user-provided/bundleable; support short relaxations and simple trajectory handling for analysis (not full MD campaign management).
- **QM/electronic structure:** PySCF as the clean bundled Python path; user-provided ORCA/Psi4 as adapters where the user has them and redistribution terms allow; support small property/charge/energy computations relevant to the analysis.

### 5.4 Job/session model (V1 foundation, V2 richer)
- Jobs have: engine, inputs, parameters, status, logs, results, timestamps.
- Results attach to the analysis artifact.
- V2: batch jobs, job templates, history, compare across jobs/poses/engines.

### 5.5 V2 additions (beyond engine expansion)
- Batch analysis of a list of structures/ligands with engine tasks where appropriate.
- 2D ligand edit + 3D sync (only if it serves the workflow).
- Water network analysis where the structure supports it.
- Scripting/console: a focused console that can drive the backend, including engine-adjacent scripting.
- Result diff/compare across poses, engines, and conditions.

## 6. Frontend feature spec (V1 + V2)

### 6.1 V1 UI
- Existing analysis panels (scene, ligand card, contact table, evidence, measurement, notes, export).
- Engine panel:
  - Engine selection (bundled vs discovered user engines, clearly labeled).
  - Task selection (dock; local geometry cleanup; later MD/QM).
  - Parameter input tied to the app state (pocket/box from ligand, ligand selection, simple parameter controls).
  - Job control: run, cancel, status, logs.
  - Results: list of jobs/poses/results; pick to visualize in-scene; include in artifact.
- Result visualization: poses in-scene; post-relaxation contacts; relevant numeric readouts.

### 6.2 V2 UI
- Comparison workspace for results across poses/engines/conditions.
- Batch results view with drill-in.
- 2D editor panel with 3D sync toggle.
- Water network layer toggle.
- Console panel for scripting.
- Job history/templates.

## 7. Engine integration details (honest, implementable)

### 7.1 Docking adapter pattern (V1 example)
- Choose the pocket/box from the app state: ligand-centric pocket as the default, with a user-adjustable box where needed.
- Prepare protein + ligand inputs in the formats the engine expects.
- Run the CLI engine in the backend; capture stdout/stderr/logs; parse poses/scores.
- Normalize into: pose coordinates, score, metadata, logs.
- Show poses in-scene; let the user pick one for further analysis/export.

### 7.2 Local geometry cleanup (V1)
- A small backend step (or a tiny engine call) that checks/relaxes the picked ligand geometry for analysis context, with clear labeling that this is not production dynamics.
- Useful to flag strain and to give a cleaner ligand for downstream visualization.

### 7.3 MD adapter pattern (V2)
- Prefer OpenMM as a library call from Python for the bundled path (GPU-capable, Python API, clean integration story).
- Support short simulations/relaxations and simple trajectory import for analysis (not full campaign management).
- For GROMACS or other CLI MD engines, treat as user-provided/bundleable per the packaging strategy.

### 7.4 QM adapter pattern (V2)
- PySCF as the clean bundled Python path for small property/charge/energy computations relevant to the analysis.
- User-provided ORCA/Psi4 where the user has them and redistribution terms permit; do not bundle license-sensitive binaries in a freely distributed Snap without verifying terms for the exact model.

## 8. Packaging & Snap (the real constraint)
### 8.1 Core Snap (strict confinement)
- App core: Tauri app + Python backend + analysis libs + PLIP-backed contact analysis + cheminformatics + embedded WebGL frontend.
- Plugs: opengl, desktop, home, network. GPU works for the scene and any GPU-backed backend work (torch, OpenMM if used locally, etc.).
- Strict confinement keeps the distributable story clean.

### 8.2 Engine bundling strategy
- Bundle the clean, redistributable engine set you can ship for amd64/arm64: PySCF (Python lib), and CLI docking tools only if their licensing and binary redistribution are acceptable for your distribution model.
- For heavier or license-sensitive engines (GROMACS fully bundled for all arches; ORCA if redistribution terms do not permit), support user-provided/discovered engines and/or a classic-snap/separate runtime path.
- Make engine discovery explicit in the UI: "bundled," "found on system," "not found — add path / install."

### 8.3 If you want "one snap with more inside"
- Classic confinement is the lever, with the tradeoff that it is not strict confinement and may face more review/questions. Choose this only if "everything inside one snap" matters more than the strict story.
- Do not assume you can cleanly pip-install big engine stacks at runtime inside a strict snap as if it were a normal environment; model the bundle + discovery + optional runtime install explicitly.

### 8.4 Licensing discipline (part of the product)
- Be accurate about each engine's license and redistribution terms in the product docs and in the Snap; do not bundle binaries whose terms you have not checked for your distribution model.
- PySCF = clean for bundling (Apache 2.0). Vina/smina/gnina = check binary redistribution terms for your model. ORCA/Psi4 = do not assume bundleability; verify per engine and per distribution model.

## 9. Data-source map (unchanged core; engines are separate local tools)
Same as before for the open-live analysis spine:
- RCSB Search/Data/GraphQL/ModelServer/File Download for structures and annotations.
- wwPDB CCD + RCSB chem_comp for ligand identity.
- PubChem + ChEMBL + PDBBind + UniProt-adjacent context for evidence.
- Local cheminformatics + geometry/contact analysis + PLIP.
- Engines are local; they read/write files in a workspace the app manages.

## 10. Implementation order (so the product is real early)
### Phase 0 — scaffold
- Repo: Tauri app + Python backend package + shared types/IPC + packaging layout (Snap + portable).
- Basic Tauri window + embedded WebGL viewer placeholder + "load local mmCIF reaches backend shows chains/ligands."

### Phase 1 — offline analysis spine (V1 core)
- mmCIF/PDB load + ligand detection.
- PLIP-backed contact analysis + contact table + pick sync + measurement + pocket shell.
- Ligand card with identity from archive data + honest status.
- Live enrichment clients (RCSB, PubChem, ChEMBL, PDBBind) + cache + offline fallback.
- Export + save/load analysis artifact.

### Phase 2 — engine orchestration V1 (small and real)
- Engine adapter framework + common result model + job/session model.
- Docking adapter (choose one redistributable CLI docking engine) + pocket/box from app state + run/monitor/parse/visualize.
- Local geometry cleanup/relaxation helper.
- Engine panel UI + result visualization in-scene + results in artifact.

### Phase 3 — polish + integration
- Representation/color/camera ergonomics, picking, export quality.
- Clean messaging for bundled vs discovered vs missing engines.
- Result diff/compare basics.

### Phase 4 — V2
- Additional docking adapters (gnina/smina where appropriate).
- OpenMM MD adapter (library path) + short simulations/relaxations + simple trajectory import for analysis.
- PySCF QM adapter (bundled) + user-provided ORCA/Psi4 discovery where terms allow.
- Batch analysis + jobs/templates/history.
- 2D edit/sync, water network, scripting console.

## 11. Acceptance checks (V1)
- Opens a representative bound PDB/CSM entry; picks ligand; shows ligand-centered analysis (scene + ligand card + contact table + optional evidence).
- Contact table pickable/correct; distances/angles correct; contact types transparent.
- Ligand identity honest; PubChem/ChEMBL/PDBBind enrichment where it exists.
- Export produces coherent artifact (scene image + contact CSV + ligand SDF/2D where available + summary/notes), reloadable in V1.
- Docking adapter runs a real dock from the app state into the pocket; results visualized in-scene and included in artifact; job status/logs/cancel work.
- Local geometry cleanup is available and clearly labeled as an analysis aid, not production dynamics.
- Works offline from a local file; network enrichment degrades gracefully.
- Bundle/discovery model for engines is clear and honest in the UI and docs.
- Snap installs, launches from app grid + CLI, renders with GPU on default Ubuntu (X11 + Wayland where feasible).
- No existing app ships this combined workflow (analysis spine + real engine orchestration + structured artifact) as one snappable product; implementation preserves that.

## 12. What "no similar app" means with engines
The gap is still the join + artifact, now with real engines in the loop:
- Viewers + separate CLI engines + separate data browsers exist.
- A snappable workstation that lets you analyze a bound structure, then run a real dock (or later MD/QM) from the same state, visualize the results in-scene, and fold them into one structured artifact — with a clear bundled vs user-provided engine model — is not a common single product.
- The product is the orchestration + integration + artifact, not the engines themselves.

## 13. Risks and mitigations
- **Engine bundling/licensing.** Verify each engine's license and redistribution terms for your Snap distribution model before bundling binaries; do not bundle what you have not cleared. Use PySCF/Python-lib engines as the clean bundled path where possible.
- **Snap + big engine stacks.** Do not pretend a strict snap can pip-install large engine stacks at runtime like a normal environment; model bundle + discovery + optional runtime install explicitly; classic-snap is the tradeoff if you really want more inside one snap.
- **MD/QM weight.** Keep V1 to one docking adapter + small geometry helper; push MD/QM to V2 and only with a realistic bundling/user-provided strategy.
- **Engine result normalization.** Different engines produce different outputs; invest early in a common result model and parsers, or the UI/artifact promise breaks.
- **Pocket/box definition for docking.** Derive a reasonable default from the ligand/pocket state; let the user adjust; do not pretend the app solves box selection perfectly.
- **Performance/runtime control.** Engine runs can be long; the app must handle job monitoring, cancel, logs, and not block the UI.

## 14. Next concrete step
Scaffold the repo (Tauri + Python backend + shared types/IPC + packaging layout), then implement the offline analysis spine end to end (Phase 1), then add the engine adapter framework + one real docking adapter + pocket-from-state + run/visualize (Phase 2). That proves the product spine and the engine integration before any V2 MD/QM.

## 15. Explicit decision points still needed
- Which docking engine(s) for V1, and are their binaries redistributable for your Snap model on amd64/arm64?
- Strict-Snap core + bundled small engine set + user-provided engines (Strategy A), or classic-snap "more inside one snap" (Strategy B)?
- Whether OpenMM counts as "bundled MD" in your model (Python lib, GPU-capable — good fit), and whether you want GROMACS too and how you will source it.
- Whether to bundle PySCF now (clean) and treat ORCA/Psi4 as user-provided only.
