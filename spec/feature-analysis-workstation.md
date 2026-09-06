# Product Spec — Ligand Interface Analysis Workstation
# Working name: "Ligora"  (pick later; not final)
# Version: 1.0 (spec draft)
# Constraint this spec is built around:
#   No existing app ships THIS workflow as one snappable, local-first product.

## 1. Product thesis
A snappable desktop workstation for one question:

> "Given a protein structure that has a bound small molecule, what is actually happening at the interface — chemically and structurally — and how does this fit into the open evidence around the ligand and target?"

This is NOT:
- a general molecular viewer (PyMOL, ChimeraX, Mol*, jlmol-style shells),
- a 2D interaction-diagram generator alone (LigPlot+, PLIP web),
- a cheminformatics/data browser alone (ChEMBL web, PubChem, PDBbind web),
- a docking/scoring/minimization suite.

It IS: a focused local-first analysis workflow that takes one bound structure + one ligand of interest and produces a **structured, shareable interface analysis** that unifies:
- the 3D scene and contact analysis,
- the ligand identity and chemistry,
- the target/evidence context from open live data,
- a comparison layer where it's useful,
- and a real exportable analysis artifact (image + structured data + notes), not just a screenshot.

Visualization is a component; the deliverable is the analysis.

## 2. Fields this serves (and why they are underserved)
- Structural biologists annotating a ligand site in a new/revised structure.
- Med-chem / comp-chem people reasoning about a binding mode from a PDB entry and wanting a coherent local write-up, not tab-switching between viewer + PDB pages + ChEMBL + notes.
- Students and teaching labs where the task is "analyze this interface," not "look at this molecule."
- Anyone who currently does this in pieces and then pastes screenshots + a hand-built table into a document.

Why this field has enough data to be real:
- PDB archive: 200k+ experimental structures + 1M+ CSMs (AlphaFold/ModelArchive), growing ~10k+/month in 2025-2026.
- RCSB Data API + Search API v2 + GraphQL + ModelServer (BinaryCIF, on-demand subsets) + File Download Services are open and queryable; ligand chemical component data is in the archive and the wwPDB Chemical Component Dictionary; RCSB exposes small-molecule ligand search.
- ChEMBL: >13M activity records, molecule/target/activity/assay/binding_site endpoints, documented REST API, substructure/similarity search.
- PubChem: public compound space with APIs and toxicity/solubility-type annotations.
- PDBBind: curated binding affinities for protein-ligand complexes (an affinity layer, not universal, but real where present).
- UniProt: target identity/function/features/disease context, tied in via RCSB annotations and its own APIs.

All are open and live-queryable. The app can work offline from a local file and use the network as an optional enrichment layer with graceful degradation.

## 3. What makes this NOT "a similar app already"
Be explicit about the boundaries and the gap.

### What exists today (and why it's not this)
- **Viewers with scripting/plugins (PyMOL, ChimeraX, Chimera, SAMSON):** can be coerced into contact analysis and some reporting, but not a coherent ligand-interface analysis flow with ligand identity + evidence context + structured export as a default artifact. Heavy to ship and distribute; not a snappable product for this workflow.
- **Web/browser viewers (Mol*, NGL in browser):** great at viewing and selecting; no offline local-file-first analysis workflow, no structured analysis export, no integrated ChEMBL/PDBBind/PubChem evidence pane in one app.
- **2D interaction-diagram tools (LigPlot+, PLIP web):** produce contact diagrams/tables; valuable, but separate from the interactive 3D scene and from ligand identity and evidence context, and not a unified snappable desktop analysis app.
- **Ligand/binding-data web resources (PDBbind web viewer, PDBeChem/PDB ligand dictionary, ChEMBL web, UniChem):** give you data about ligands and affinities, but not the local 3D scene + contact analysis + write-up.
- **Electron+JSmol/JSME desktop shells (jlmol-style):** wrap a viewer/editor; not this analysis workflow, and heavy for the Snap story.

### The gap this product owns
No existing product to our knowledge:
1. opens a bound structure (local file or PDB ID / CSM-ID) **and** immediately frames the analysis around a chosen ligand,
2. computes and presents interface/contact analysis **in the 3D scene AND as an inspectable structured table**,
3. resolves the ligand's identity and chemistry **from the archive/CCC dictionary and live compound sources** into one ligand card,
4. pulls the **evidence context** (affinity where available via PDBBind; bioactivity/target context via ChEMBL; compound identity/synonyms/annotations via PubChem) into the same view,
5. gives a **structured analysis export** (scene image(s) + contact table + ligand identity/export + summary/notes) as a first-class artifact,
6. does all of the above as a **snappable, local-first, offline-capable desktop app** with GPU rendering in strict confinement.

That combination is the product. Each piece exists somewhere; nobody ships them as one workflow for this use case.

### Honest edge cases where "no similar app" is weakest
- Heavy users of PyMOL/ChimeraX + scripting can approximate parts of this; we are not targeting them first. We target the broad user base that wants one coherent local analysis flow, not a modeling environment.
- PLIP and LigPlot+ do 2D contact analysis well; we complement, not replace — and we put the analysis next to the scene and the evidence, which they do not.
- ChEMBL/PubChem/PDBbind web tools give ligand/target/affinity context; we bring that context to the structure, and allow it to be optional and local-cached.

## 4. Core user workflow (the spine)
### 4.1 Inputs
- **Live:** a PDB ID or CSM-ID (e.g. a real PDB code or an AF-model ID). Resolved via RCSB Search/Data API.
- **Local:** an existing structure file (PDB, mmCIF preferred; PDB legacy as fallback). Open via file picker / drag-and-drop / command-line argument.
- **Ligand choice:** when a structure has multiple non-polymer entities, the app identifies candidate small-molecule ligands and lets the user pick the one to analyze. This is the analysis center.

### 4.2 First-view analysis
For the chosen ligand, the app immediately presents:
- The 3D scene: protein (cartoon/ribbon by chain as default), ligand (ball-and-stick, atom-label on pick/hover), the pocket on demand (surface), and contact highlights.
- A **ligand card**: resolved identity — chemical component / residue name, formula, molecular weight, 2D representation where available, SMILES/SDF where resolvable, stereochemistry/flag notes, classification hint (drug-like small molecule vs cofactor/ion/solvent — auto-classify with user override).
- A **contact table**: residues in contact with the ligand, by chain, with distance and a first-pass contact-type label; pickable rows highlight the residue in the scene.
- An **evidence summary pane**: what is known about this ligand and/or target from open sources (PubChem identity/synonyms/annotations where mapped; ChEMBL bioactivity where the ligand maps; PDBBind affinity where the complex/ligand is covered; related structures via RCSB Search where useful).

### 4.3 Analysis operations (the "deep analysis" part)
- **Contact inspection:** distance, angle for chosen contacts; contact-type classification heuristics (hydrogen bond, hydrophobic, pi-stacking, salt bridge, halogen bond, metal coordination, water-mediated where the data/structure supports it). Transparent, inspectable, not a black box.
- **Pocket definition:** a ligand-centric pocket shell on demand; option to widen/narrow; surface coloring by contact type or by property where relevant.
- **Measurement & picking:** pick atoms/residues; measure distances/angles; highlight and annotate selections.
- **Ligand identity resolution:** from the archive chemical component data and, where possible, live compound sources. Honest statuses: resolved / partial / not found — no fake confidence.
- **Evidence context lookup (live, optional):** 
  - PubChem mapping for the ligand when available (name, synonyms, annotations like toxicity/solubility where present).
  - ChEMBL mapping for the ligand where available (bioactivity summary, target context).
  - PDBBind affinity where the complex/ligand is covered.
  - Related structures via RCSB Search (same ligand, similar ligand, sequence/3D neighbors) as optional context, not clutter.
- **Comparison (small but high-value):**
  - Open a second structure/ligand and overlay/side-by-side for same-ligand or similar-ligand comparison.
  - Simple ligand similarity context where the app has a cheminformatics path (fingerprint + Tanimoto) for "how does this compare to known binders."
- **Annotation & notes:** user notes attached to the analysis; selection/annotation persistence for a session and into the export.

### 4.4 Output — the analysis artifact
This is what makes it analysis, not viewer.
- **Scene images:** screenshot at viewport resolution; optional higher-quality render of the interface/pocket.
- **Contact table export:** CSV/structured of residues, chains, distances, contact types, atom pairs.
- **Ligand export:** SDF of the resolved ligand where available; 2D image where available.
- **Analysis summary:** a structured summary (JSON/text) combining ligand identity, contact summary, evidence summary, and user notes — reloadable and shareable.
- **Notes/annotations** preserved so a user can return and continue.

### 4.5 Offline behavior
- With no network: open local files; do geometric contact analysis; render; export scene images + contact table + ligand export from the file. No live PubChem/ChEMBL/PDBBind/PDB enrichment. Clearly indicate what is unavailable.
- With network: enrich as above. Rate-limit/backoff for the APIs; do not assume unlimited connectivity.

## 5. Feature spec (MVP → V1 → V2)

### MVP (shipable first version)
**Opening & scene**
- Open PDB ID / CSM-ID (live) or local PDB/mmCIF file.
- Auto-detect candidate non-polymer ligands; user picks the one to analyze.
- Default scene: protein cartoon/ribbon by chain; ligand ball-and-stick with atom labels on pick/hover.
- Mouse/trackpad navigation: rotate, zoom, pan; atom/residue picking.
- Representations: cartoon/ribbon, stick (for sidechains near ligand), ball-and-stick (ligand), optional ligand-centric surface on demand.
- Coloring: by chain; by contact status (contacting residues highlighted); ligand-centric palette.

**Analysis**
- Contact detection: residues within contact distance of the chosen ligand; per-contact distance; first-pass contact-type heuristics.
- Contact table: chain, residue, atom, distance, contact type; row pick highlights in scene.
- Measurement: distance and angle for selected atom pairs; clear readout.
- Pocket surface option: ligand-centric shell, toggle on/off.
- Ligand card: chemical component/residue name, formula, MW, 2D where available, SMILES/SDF where resolvable, classification hint, resolution status.

**Evidence (live, optional)**
- PubChem mapping for the ligand where available (identity/synonyms/annotations where present).
- ChEMBL bioactivity summary where the ligand maps.
- PDBBind affinity where the complex/ligand is covered.
- Clear "offline" state when none of this is available.

**Export**
- Scene screenshot (PNG).
- Contact table (CSV) with the data behind the view.
- Ligand SDF where resolvable; 2D image where available.
- Analysis summary (JSON/text) with ligand identity, contact summary, evidence summary, notes.

**Session & UI**
- Window chrome: open/search, representation & color controls, ligand card, contact table, evidence pane, measurement tools, export.
- Notes/annotations for the analysis.
- Preferences: theme (dark default, light option), UI scale, default representation.
- Offline-first posture; live lookups are enrichment.

**Snap/package**
- Strict confinement with GPU (opengl), desktop, home (for user files), network (optional enrichment).
- Works offline from a local file; network lookups degrade gracefully.
- Target amd64 + arm64; publish to Snap Store; secondary portable artifact later.

### V1
- Extended contact classification: HB, hydrophobic, pi-stacking, salt bridge, halogen bond, metal coordination, water-mediated where sensible.
- Comparison views: same ligand across structures; side-by-side/overlay; pocket conservation summary.
- Ligand similarity context via fingerprint + Tanimoto where the cheminformatics path exists.
- Variant/mutation annotation on pocket residues where data is available (UniProt/PDB annotations).
- Reloadable analysis artifact: re-open a saved analysis and restore scene state, selections, notes.
- Better render options: ambient occlusion toggle, AA quality, better export resolution.

### V2 (only after the workflow is proven)
- Batch analysis of a list of PDB IDs / local set with a comparative summary (med-chem/series use case).
- 2D ligand editing/building with 3D sync ONLY if it serves the analysis workflow (not a general editor).
- Water-mediated network analysis where the structure supports it.
- Scripting/console for power users (small, focused), clearly below the main workflow.

## 6. Data-source map (live + offline)

### 6.1 Structure & coordinates (the scene)
- **RCSB PDB / wwPDB archive** — experimental structures + CSMs.
  - Search API v2 to resolve PDB ID / CSM-ID and find candidate ligands / related structures.
  - Data API / GraphQL for per-entry annotations (chains, entities, ligands, chem_comp, observations).
  - ModelServer for coordinates (BinaryCIF, subsets by assembly/chain/ligand) — preferred for on-demand, efficient fetching.
  - File Download Services for full mmCIF/PDB files when needed or offline cache.
- **Local file** — PDB (legacy) and mmCIF as primary; the app treats mmCIF as first-class now.

### 6.2 Ligand identity & chemistry
- **wwPDB Chemical Component Dictionary (CCD)** — standard chemical descriptions for ligands/small molecules in the archive.
- **RCSB chem_comp / core_nonpolymer_entity** — ligand identity from the entry: residue/name, formula, etc.
- **PubChem** — compound identity, synonyms, annotations (toxicity, solubility-type data where present) when the ligand maps.
- **ChEMBL molecule endpoint** — molecule properties, structures (SMILES/InChI), 2D image, synonyms, xrefs where the ligand maps.

Status handling is part of the product: "resolved", "partial", "not found" — not fake confidence.

### 6.3 Evidence context (binding activity, target, affinity)
- **PDBBind** — curated binding affinities for protein-ligand complexes (affinity layer where present; not universal).
- **ChEMBL** — activity (binding/functional/ADMET-type), target, assay, binding_site; similarity/substructure search; molecule/target/activity/assay endpoints.
- **UniProt** — target identity/function/features/disease context; tied via RCSB annotations and its own APIs.
- **RCSB Search API** — related structures (same ligand, similar ligand, sequence/3D neighbors) as optional context.

### 6.4 Local/computational (owned by the app)
- **Geometry/contact analysis** — contact detection, distances, angles, contact-type heuristics from coordinates. This is app logic; it is the "analysis" core.
- **Cheminformatics** — for ligand identity normalization, simple fingerprints/Tanimoto for similarity context (V1), SDF/SMILES handling where possible.
- **Caching** — optional local cache of enrichments so repeated lookups are fast and the app feels local-first.

### 6.5 Offline-first reality
- Core analysis (open local file, contact geometry, scene, export of scene image + contact table + ligand export) works fully offline.
- Live enrichment (PubChem/ChEMBL/PDBBind/related-structures) is an optional network layer; absent network = app still useful, just without enrichment.
- API rate limits/backoff handled; do not treat the live layer as required for the product to function.

## 7. Non-goals (explicit, to hold the line against "similar app" creep)
- Not a general viewer replacement for PyMOL/ChimeraX.
- Not a docking/scoring/minimization/QM suite.
- Not a full cheminformatics platform or a drug-discovery data warehouse.
- Not a 2D-only interaction-diagram tool (we complement LigPlot+/PLIP, we don't duplicate them as the main event).
- Not VR/AR/collaboration for now.
- Not a batch data-mining platform in V1 (batch analysis is V2 and only if the single-structure workflow proves out).

## 8. Acceptance criteria (first shipped version)
- Opens a representative bound PDB structure and a bound CSM entry; lets the user pick the ligand and immediately see a focused ligand-centered analysis (scene + ligand card + contact table + optional evidence where network allows).
- Contact table is inspectable and pickable; distances/angles read correctly; contact-type labels are transparent and editable in interpretation (heuristics, clearly labeled).
- Ligand identity resolution is honest about what is and isn't available; PubChem/ChEMBL/PDBBind enrichment appears where it exists and is clearly absent otherwise.
- Export produces a coherent analysis artifact: scene image + contact table CSV + ligand SDF/2D where available + summary/notes; the artifact is reloadable in V1.
- The app is useful offline from a local file; network lookups are enrichment and degrade gracefully.
- `snap install` installs, launches from the app grid and CLI, renders with GPU on default Ubuntu (X11 and Wayland where feasible), strict confinement works for the chosen plugs.
- No existing app is doing this exact combined workflow as one snappable product; the spec's boundaries (scene + contact analysis + ligand identity + evidence context + structured export, local-first) are real and hold.

## 9. Why this is defensible as a product (not a feature)
Because the value is the **join and the artifact**, not any single capability:
- A viewer shows the scene; a contact tool gives a table; a ligand dictionary gives chemistry; ChEMBL/PDBbind/PubChem give evidence; a browser viewer gives some of this in one web page but not a local offline-first analysis app with a structured export.
- This product takes one structure + one ligand, runs the geometric analysis locally, enriches with open live data where available, and hands the user a shareable analysis object. That join + artifact is what doesn't exist as a snappable product today.

## 10. Next (out of scope for this doc)
- Stack/lane decision (Rust GPU renderer + owned geometry/cheminformatics core vs Electron+JSmol/JSME lane) and architecture sketch.
- Snap shape for the chosen lane (GPU + desktop + home + network, strict confinement).
- Concrete API usage plan: exact RCSB Search/Data/ModelServer calls, ChEMBL endpoints to hit first, PubChem mapping strategy, PDBBind usage, and caching strategy.
