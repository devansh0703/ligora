# Full Implementation Plan — Ligora
# x86_64 only. One install = fully working UI. Open-source engines only, used (not reimplemented).
# Stack: Tauri (Rust shell) + Python backend (torch where useful, PLIP, cheminformatics, enrichment, engine adapters)
#        + embedded 3Dmol.js WebGL viewer.
# Distribution: one Snap install on x86_64 that ships a usable workstation, with engines available in the same
#               install experience (bundled where open-source license allows; otherwise part of the same install flow).

## 1. Product in one line
A snappable x86_64 workstation that, in one install, opens a bound structure, lets you pick a ligand, analyze the
interface (scene + contacts + ligand identity + open-live evidence), and run real open-source engines (docking in V1;
MD/QM in V2) from the same state, with results visualized in-scene and folded into one structured artifact.

No existing app ships this combined workflow (analysis spine + real open-source engine orchestration + structured artifact)
as one snappable product. That is the gap.

## 2. Engine set (all open source; bundlability noted)
This app does NOT reimplement docking/MD/QM. It orchestrates existing open-source engines through adapters.

### 2.1 V1 engine set
- **Docking:** AutoDock Vina (Apache 2.0). Clean to bundle/redistribute. This is the V1 docking anchor.
- **Local geometry cleanup:** backend geometry / simple minimizer (owned code, open license), used as an analysis aid,
  not production dynamics.
- **Contact analysis:** PLIP (GPL) + app geometry. GPL is fine to use; just don't pretend PLIP is the only contact path
  and be aware of copyleft if you bundle/depend on it in a way that triggers GPL obligations for the whole snap.
- Embedded viewer: **3Dmol.js** (BSD-3-Clause) for the WebGL scene.

### 2.2 V2 engine set (open source)
- **Docking expansion:** gnina (dual GPL/Apache) as an additional adapter where the license model is accepted; otherwise
  keep Vina as the primary and add features around it.
- **MD:** OpenMM (Python, open/available; LGPL-ish contingent on exact packaging) as the library adapter; alternatively
  a CLI MD engine if you want one, but OpenMM is the best Python-integrated open path.
- **QM/electronic structure:** PySCF (Apache 2.0) as the clean bundled Python path; other open engines only if their
  license model is acceptable for bundling.

### 2.3 Perpetual rule
Every engine in the snap must be open source and its license must be satisfied for the way you ship it. This plan assumes
you will actually verify each engine's license terms for the exact distribution model (bundle vs adapter vs process)
before shipping. The plan does not hardcode a license conclusion beyond the facts above; it builds the adapter framework
so you can plug engines in and out as license/bundle decisions are finalized.

## 3. Packaging model (one install = fully working)
### 3.1 The user-facing requirement
User installs the app snap, opens it, and can work immediately. All features the snap ships are usable without a separate
engine hunt.

### 3.2 How engines get there
Two legitimate open-source-compatible ways, and you can mix them:

- **Bundle into the app snap** what is cleanly redistributable under its license and acceptable for your overall snap
  license/structure. Vina (Apache 2.0) and PySCF (Apache 2.0) and 3Dmol.js (BSD) are examples of code you can ship inside
  the snap under the appropriate notices. OpenSSL/numpy/scipy-type Python deps are ordinary dependencies, handled like any
  Python app.
- **Bring in via the same install flow** for engines that are open source but where bundling has packaging/license
  implications you want to handle separately. In Snap terms this can mean provider snaps consumed via content interface, or
  classic confinement for the app snap if that simplifies engine access, or a bundled-but-separate runtime. The point is the
  user still ends up with a working install, not a "go install this yourself" step.

### 3.3 Recommended default for this project
Make the app snap **classically confined on x86_64** if that is what lets you ship one install = fully working with the
engine set you want, without fighting strict-confinement limitations on engine access and runtime layout. Classic confinement
is the tradeoff that buys you the user experience you asked for. Keep strict confinement as an option you revisit if/when you
want a stricter store story.

Why this is defensible: the product is the orchestration + analysis + artifact, not the engines themselves. The engines are
open source. Classic confinement on x86_64 is a reasonable way to ship one working install when the alternative is forcing a
strict model that complicates engine integration for no user benefit.

### 3.4 Licensing hygiene for the whole snap
- Include proper notices for each bundled open-source component (Vina, PySCF, 3Dmol.js, PLIP if bundled, OpenMM if bundled,
  Python deps, etc.).
- Be accurate in the snap metadata/docs about which engines are included and under what conditions.
- If you bundle GPL code (PLIP, or gnina with its GPL side), understand the implications for the rest of the snap and decide
  the overall license posture accordingly. The plan keeps the backend/adapter boundary clean so you can make those decisions
  without rewriting the whole product.

## 4. Architecture (real, no stubs)
### 4.1 Repo layout
- `frontend/` — Tauri webview app (3Dmol.js scene + analysis/evidence/engine UI).
- `backend/` — Python package: analysis core, enrichment clients, engine adapters, job/session, export, CLI.
- `snap/` — snapcraft project for x86_64.
- `shared/` — schemas/types used by frontend and backend over IPC (JSON for commands/results; file artifacts for blobs).
- `ci/` — CI config, basic tests, and acceptance helpers.

### 4.2 Runtime shape
- Tauri app launches and hosts the webview.
- Python backend runs as a local service/process owned by the app; Tauri talks to it over a local JSON IPC channel.
- Backend manages:
  - workspace files (input files, engine inputs/outputs),
  - analysis state,
  - jobs,
  - results normalization,
  - enrichment cache.
- Engines are invoked by the backend from the installed/bundle layout; outputs are parsed into the common result model.

### 4.3 IPC model
- Frontend sends JSON commands (open structure, pick ligand, run contact analysis, run docking, set parameters, start job,
  cancel, export, save/load artifact, etc.).
- Backend returns JSON results and events (job status, logs, result ready, errors).
- Larger artifacts (SDF, CSV, images) are saved to the workspace and referenced by path/URL over the IPC, not payload-inlined
  into every JSON message.

## 5. Backend feature implementation (real, no stubs)

### 5.1 Workspace + file management
- Each analysis session has a workspace directory with:
  - the loaded structure file(s) (local or fetched),
  - the chosen ligand file(s),
  - engine input/output files,
  - analysis artifacts (contact CSV, ligand SDF, scene image, analysis summary JSON),
  - job logs.
- Backend owns workspace creation, path management, cleanup, and artifact export.

### 5.2 Structure parsing and ligand detection
- Accept mmCIF (primary) and PDB (legacy fallback) from local file or fetched source.
- Parse polymers and non-polymer entities; identify candidate small-molecule ligands by entity type and chemical component
  data where available.
- Provide downstream consumers with: chains, residues, atoms, ligand entities, ligand atom sets, coordinates.

This is real parsing code, not a stub. For mmCIF, use a real mmCIF parser path (a Python mmCIF/PDBx parser) rather than
hand-splitting text. For PDB legacy, handle the fixed-column format correctly. In both cases, normalize to a common internal
structure model.

### 5.3 Ligand identity resolution
- From the archive chemical component data for the entry (chem_comp / non-polymer entity data), build the ligand card:
  component name/residue name, formula, molecular weight, atom counts, and any 2D/SMILES/SDF where available/resolvable.
- Attempt mapping to PubChem/ChEMBL where possible; report honest status (resolved / partial / not found).

### 5.4 Contact analysis (PLIP + owned geometry)
- Run PLIP on the structure+ligand to classify non-covalent interactions.
- Present a contact table: chain, residue, atom, distance, contact type, ligand atom pair, and any PLIP-internal detail you
  choose to expose.
- Supplement with owned geometry: distance/angle queries, contact selection, pocket shell, and any contact-type heuristics
  you want beyond PLIP.
- Pickable: contact row selection highlights the corresponding atoms/residues in the scene.

PLIP is invoked from the backend against the workspace structure/ligand files. The backend parses PLIP output and normalizes
it into the contact table. If PLIP is unavailable for a reason, the app falls back to owned geometry and clearly reports it.

### 5.5 Cheminformatics
- Normalize ligand identity where possible (canonical SMILES/InChI handling, SDF read/write).
- Compute fingerprints and Tanimoto similarity for ligand similarity context (V1-lite, V2 fuller).
- SDF export of the resolved ligand where available.

### 5.6 Live enrichment clients (optional, network, cache, honest)
- **RCSB:** Search API v2 to resolve PDB/CSM-ID and find candidates/related structures; Data API/GraphQL for annotations;
  ModelServer/Download for coordinates.
- **PubChem:** compound identity/synonyms/annotations where the ligand maps.
- **ChEMBL:** molecule/target/activity/binding_site context where the ligand/target maps; similarity/substructure where useful.
- **PDBBind:** affinity where the complex/ligand is covered.
- Cache fetched enrichments per session and optionally beyond; clearly mark what is live vs cached vs unavailable.

### 5.7 Enrichment client implementation notes (real)
- Each client is a small bounded module with retry/backoff and rate-limit awareness.
- Clients return structured data and a status; the backend merges them into the evidence pane with honest labeling.
- Network is optional; offline is a first-class mode.

### 5.8 Engine adapter framework
- Common adapter contract:
  - discover/availability (bundled path, or system path, or provider snap path),
  - prepare inputs from app state (structure, ligand, pocket/box, parameters),
  - run (blocking or async with monitoring),
  - parse outputs into a common result model,
  - provide logs/status/cancel support.
- Common result model:
  - job id, engine, parameters, status, started/finished/errors,
  - primary results: e.g. poses+scores for docking, energies/properties for QM, relaxation/trajectory info for MD,
  - files: input files, output files, parsed artifacts,
  - logs.

### 5.9 V1 docking adapter (Vina, real)
- Detect Vina availability from the installed/bundle layout.
- Prepare protein + ligand inputs in Vina-ready form from the workspace structure/ligand files.
- Derive a docking box from the app state: ligand-centric pocket as the default, with a user-adjustable box where the UI
  exposes it.
- Run Vina with the chosen parameters; capture stdout/stderr; parse poses and scores.
- Normalize into: pose coordinates, affinity/score, metadata, logs, output files.
- Visualize poses in-scene; let the user pick one for further analysis/export.

### 5.10 V1 local geometry cleanup
- A small owned backend step that checks/relaxes the picked ligand geometry enough to support analysis context (strain flag,
  cleaner ligand for visualization).
- Clearly labeled as an analysis aid, not production dynamics.

### 5.11 Job/session model
- Sessions own analysis state + workspace.
- Jobs belong to a session; they have engine, inputs, parameters, status, logs, results, result files.
- UI shows job list, progress, logs, result ready events, cancel.
- Results attach to the analysis artifact.

## 6. Frontend feature implementation (real, no stubs)

### 6.1 Scene (3Dmol.js)
- Embed 3Dmol.js in the Tauri webview.
- Feed parsed structure data (coordinates + representation directives) from the backend or from parsed local state.
- Representations: protein cartoon/ribbon by chain; ligand ball-and-stick with atom labels on pick/hover; stick for
  near-ligand sidechains; ligand-centric surface on demand.
- Coloring: by chain; by contact status; ligand-centric palette.
- Camera: rotate/zoom/pan; a few preset views for ligand/pocket/interface.

### 6.2 Interaction
- Atom/residue picking drives the ligand card, contact table, and scene highlight.
- Contact row selection highlights atoms/residues in the scene.
- Distance/angle measurement with readout.

### 6.3 Panels
- **Structure/ligand opener:** PDB ID / local file / recent; ligand candidate selection.
- **Ligand card:** identity, formula, MW, 2D where available, SMILES/SDF where available, mapping status, classification
  hint.
- **Contact table:** chain, residue, atom, distance, contact type, ligand atom pair; sort/filter basics; pickable.
- **Evidence pane:** PubChem/ChEMBL/PDBBind/RCSB context where available; honest offline labeling.
- **Measurement tools:** distance/angle.
- **Engine panel:** engine selection (bundled/available), task selection (dock; local geometry cleanup; later MD/QM),
  parameters, job control, logs, results list.
- **Notes/annotations:** free notes + structured annotations attached to the analysis.
- **Export:** scene image, contact CSV, ligand SDF/2D where available, analysis summary; save/load analysis artifact.

### 6.4 Job/results UX
- Start job from the engine panel; show status and logs; result-ready events update the results list and scene.
- Results visualized in-scene where applicable (e.g., docking poses).
- Selected results become part of the artifact.

### 6.5 Artifact
- Analysis summary plus scene images, contact CSV, ligand export where available, job results where present, notes.
- Save/load artifact to resume a session.

## 7. V1 scope (shipable)
- Open PDB ID / local mmCIF/PDB; ligand candidate selection.
- Scene (3Dmol.js) with representations/coloring/picking.
- Contact analysis (PLIP + geometry), contact table, pick sync, measurement, pocket shell.
- Ligand card with honest identity/mapping status.
- Live enrichment (RCSB/PubChem/ChEMBL/PDBBind) with offline fallback.
- Vina docking adapter with box-from-state + run/monitor/parse/visualize.
- Local geometry cleanup helper.
- Export + save/load artifact.
- Job/session model with logs/status/cancel.
- x86_64 Snap that installs and opens into a working UI with the above usable.

## 8. V2 scope (after V1 is real)
- gnina docking adapter (where license model accepted).
- OpenMM MD adapter (library path) with short simulations/relaxations and simple trajectory import for analysis.
- PySCF QM adapter for small property/charge/energy computations relevant to analysis.
- Batch analysis of a list of structures/ligands with engine tasks where appropriate.
- Job templates/history; result diff/compare across poses/engines/conditions.
- 2D ligand editor with 3D sync toggle (only if it serves the workflow).
- Water network layer toggle where the structure supports it.
- Scripting/console panel for backend-driven automation.

## 9. V2 MD/QM implementation notes
- MD via OpenMM as a Python library call from the backend is the most natural open integration here; keep it to short
  simulations/relaxations relevant to analysis, not campaign management.
- QM via PySCF as a Python library call for small computations; be explicit that this is small-scale analysis support.
- For any CLI engine you add, reuse the adapter contract: discover, prepare, run, parse, normalize.

## 10. Snap build (x86_64, one install = workable)
### 10.1 Project shape
- `snap/snapcraft.yaml` + parts for the Tauri app, the Python backend, and the bundled open-source engine components.
- Build for amd64.

### 10.2 Parts
- App part: Tauri frontend + Rust shell, producing the app binary.
- Backend part: the Python backend package with its dependencies.
- Engine parts (as needed): Vina and/or gnina binaries and/or PySCF/OpenMM and/or other open components, each only if their
  license supports the way you bundle them and you want them inside the one install.
- Shared/lib parts as needed (e.g., system libs the backend/engines need on x86_64).

### 10.3 Confinement
- Classic confinement as the default for this project if it is what lets one install = fully working with the engine set you
  ship. Use strict only if/when you decide the store/confinement story matters more than that.

### 10.4 Apps entry
- One or more app entries in the snap: the main GUI app, and optionally a CLI entrypoint that opens the same backend/workspace
  for headless use.

### 10.5 Engine layout inside the snap
- The backend knows where engines live in the snap installation and uses those paths; the adapter discovery layer is
  explicit about bundled paths vs external discovery.
- Do not require the user to install engines separately after the app install.

### 10.6 Licensing in the snap
- Include a notices/licenses bundle for the open-source components in the snap.
- Keep the backend/adapter boundary so you can change bundle decisions without rewriting the product.

## 11. CI + tests + acceptance (real)
### 11.1 CI
- Build the frontend and backend on x86_64.
- Build the snap on x86_64.
- Run backend unit tests (parsing, ligand detection, cheminformatics, enrichment client unit behavior, adapter input prep,
  output parsing on toy data).

### 11.2 Tests (real, not stub)
- Structure parsing tests on small mmCIF/PDB examples.
- Ligand detection tests on entries with multiple entities.
- Contact analysis tests using PLIP output parsing on sample structures.
- Cheminformatics tests on SMILES/SDF round-trips and fingerprint/Tanimoto on known examples.
- Enrichment client tests with recorded responses / mocked HTTP for offline CI.
- Engine adapter tests: input preparation and output parsing on sample Vina outputs.
- Export tests: contact CSV, SDF, analysis summary JSON written and reloaded.

### 11.3 Acceptance checks (what "done" means)
- Install the snap on x86_64; open it; open a representative bound PDB entry; pick a ligand; see the ligand-centered analysis.
- Contact table correct and pickable; distances/angles correct; contact types transparent.
- Ligand identity honest; enrichment present where available.
- Vina dock runs from the app state into the pocket; poses visualized and included in the artifact.
- Export produces a coherent artifact; save/load works.
- Offline mode from a local file works for the core analysis.
- No engine hunt after install; the UI is immediately usable.

## 12. Implementation order (so the product is real early)
1. Repo scaffold + shared types + backend package skeleton + Tauri app shell with 3Dmol.js embedded.
2. Workspace + file management + IPC.
3. Structure parsing + ligand detection + ligand card.
4. Contact analysis (PLIP + geometry) + contact table + pick sync + measurement + pocket shell.
5. Cheminformatics + ligand export.
6. Enrichment clients (RCSB/PubChem/ChEMBL/PDBBind) + cache + honest offline labeling.
7. Engine adapter framework + Vina docking adapter + job/session model + engine panel UI.
8. Local geometry cleanup helper.
9. Export + save/load artifact + notes.
10. Polish: representations/color/camera, comparison-lite, export quality.
11. Snap build for x86_64; CI; tests; acceptance.
12. V2: gnina, OpenMM, PySCF, batch, templates/history, diff/compare, 2D edit/sync, water network, console.

## 13. Lic — and what you must still decide before build
- Confirm Vina is acceptable to bundle (Apache 2.0 says it is, but you still include proper notices).
- Decide on gnina: useful but GPL-coupled; include only if you accept the license consequences.
- Decide on PLIP bundling vs use-as-external-tool; GPL applies if you bundle/depend in a triggering way.
- Decide on OpenMM and PySCF for V2 bundling; PySCF is clean (Apache 2.0), OpenMM needs careful packaging/license handling.
- Decide overall snap license posture based on the mix, especially if GPL components are bundled.
- Decide classic vs strict for the app snap on x86_64; default recommendation here is classic if it preserves one-install =
  fully working.

## 14. Next concrete step
Scaffold the repo with the real layout, shared IPC types, the Tauri shell with 3Dmol.js embedded, and the Python backend
package skeleton with workspace management and IPC. Then implement the parsing/ligand/contact/enrichment spine end to end
before the engine adapter, so the analysis workstation is real first. After that, add Vina + job model + engine panel, then
the snap build and CI/acceptance.

This plan goes all the way from repo to snap to CI/acceptance and explicitly keeps engines as used, open-source, and not
reimplemented, with one install = fully working UI on x86_64.
