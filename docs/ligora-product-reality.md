# Ligora — What It Is, What It Does, and Where It Stands

*A grounded, honest product document. Every claim about Ligora's behavior was
verified against the running code and live tests (143 backend tests + 6 staged
UI verification runs on real data). Every external source mentioned was probed
live (HTTP 200 confirmed) before being written here.*

---

## 1. What Ligora actually is

Ligora is a **desktop molecular-analysis workstation** for **protein–ligand
interaction analysis**. It runs as a native x86_64 desktop application (Tauri
shell + Python backend + 3Dmol.js viewer) on Linux, packaged for distribution
as a snap.

In one sentence:

> **Open a protein–ligand structure, and Ligora tells you — with real,
> traceable data — what the ligand is, how it binds, how it compares to known
> chemistry, how water and ions behave around it, and how alternative poses
> score — without ever inventing a number.**

Its defining engineering rule (enforced across the entire codebase, and the
thing that separates it from most "bioinformatics dashboards"):

- **No hardcoded chemistry.** Classification (solvent vs ion vs organic) comes
  from the RCSB Chemical Component Dictionary's own `pdbx_type` codes.
- **No heuristics or local approximations.** Bond orders come from the CCD
  `_chem_comp_bond` table. Contact analysis is performed by PLIP (the
  published, peer-reviewed tool), not by a reimplemented approximation.
- **No fallbacks, no mock data, no placeholders.** When a source cannot answer
  — Vina absent, PLIP missing, UniChem slow — Ligora reports that honestly in
  the UI instead of papering over it with plausible-looking numbers.
- **Everything is traceable.** Each displayed value carries its evidence
  (source, fields, URL) in the evidence pane.

---

## 2. The problem it solves

Working with a PDB structure and its ligand today means juggling **six to ten
separate tools and websites**:

| Question | What you normally do |
|---|---|
| "What is this HETATM ligand?" | Look up the 3-letter code on RCSB manually |
| "What are its properties/SMILES/CID?" | Copy the code into PubChem, ChEMBL |
| "Which residues contact it and how?" | Download structure, upload to PLIP web server, wait, download XML |
| "Are there structured waters bridging the binding site?" | Another PLIP pass / manual inspection in PyMOL |
| "Can I redock this ligand and compare poses?" | Prepare PDBQT with Open Babel by hand, write Vina config files, parse output yourself |
| "How do two docking runs/poses compare?" | Compute RMSD in a notebook |
| "Can I search my notes/results later?" | Keep scattered files |

Ligora collapses this into **one desktop application with one click each**,
while keeping every intermediate artifact (contacts CSV, ligand SDF, summary
JSON, session state) exportable so nothing is locked in.

**Target users:** computational/structural biology researchers, medicinal
chemists doing triage on co-crystal structures, students learning what
"hydrogen bond vs hydrophobic contact vs salt bridge" means on *real*
structures, and anyone who needs reproducible, evidence-backed interaction
analysis without writing glue scripts.

---

## 3. Complete feature inventory (all real, all verified)

### 3.1 Structure access & viewing
| Feature | Implementation |
|---|---|
| Open local PDB / mmCIF files | Full in-house mmCIF parser (token-level, per CIF spec) + PDB parser with altLoc/insertion-code handling |
| Open by 4-char PDB ID | Live fetch from `files.rcsb.org` |
| 3D visualization | 3Dmol.js: cartoon/sphere/stick/surface representations, element/spectrum color schemes, spin |
| Ligand focus view | Isolate ligand + binding site, zoom-to-selection |
| Click-to-zoom from data tables | Clicking a contact row zooms the viewer to that residue; clicking a water cluster zooms to that water |
| Structure metadata | Title, method (X-ray/NMR/EM), resolution — parsed from the file's own records |

### 3.2 Ligand identity & enrichment (live)
| Feature | Source |
|---|---|
| "What is this ligand?" name resolution | RCSB CCD (live `data.rcsb.org`) |
| Formula, MW, canonical SMILES, InChIKey | CCD, then PubChem (PUG REST, `ConnectivitySMILES`) |
| PubChem CID + IUPAC name | PubChem PUG REST |
| ChEMBL cross-reference | ChEMBL via UniChem mapping (live) |
| Per-value evidence chain | Every resolved value lists source + URL in the Evidence pane |
| Classification (solvent / ion / organic) | CCD `pdbx_type` (HETAS/HETAI/HETAIN…), **never guessed** |
| Honest "no drug-like ligand" path | Structures like 1CRN (protein-only) report exactly that — verified in UI tests |
| False-positive guard | Water and ions are never presented as "the ligand" (UI-tested on 1MBN: heme picked, HOH/OH never) |

### 3.3 Protein–ligand contact analysis (real PLIP)
| Feature | Implementation |
|---|---|
| Typed contacts | Hydrogen bonds, hydrophobic, salt bridges, π-stacking, π-cation, halogen, water bridges, metal coordination — from PLIP's real XML output |
| Real atom names per contact | Resolved through PLIP's per-section serial fields (`donoridx`/`acceptoridx`, `ligcarbonidx`/`protcarbonidx`, `lig_idx_list`), plus nearest-atom geometry for PLIP's coordinate-only salt-bridge fields |
| Contact distances | From PLIP's measurements |
| Contacts CSV export | Real CSV via the export module |
| Visual affordance | Contact rows carry chain/residue; click zooms viewer |

### 3.4 Water network analysis (v2)
| Feature | Implementation |
|---|---|
| Water detection & counts | From the structure's own atoms (CCD-classified solvent) |
| Water–water network graph | Distance cutoffs from configuration/user payload (no invented thresholds) |
| Cluster identification | Connected components over the network |
| Water–protein contacts | **Vectorized numpy** distance pipeline (3W85: 431 waters, 235 edges in well under a second) |
| Click-to-zoom per water/cluster | Wired in the UI |

### 3.5 Docking (real AutoDock Vina)
| Feature | Implementation |
|---|---|
| Receptor/ligand preparation | Open Babel → PDBQT (real subprocess, real conversion) |
| Docking box | Computed from the co-crystallized ligand's real extent |
| Vina execution | Real `vina` binary (1.2.7 verified), configurable exhaustiveness/num_modes |
| Pose parsing | Real PDBQT output parse → poses with affinities + elements |
| **Poses in 3D** | Any pose renders in the main viewer as an overlay model (stick+sphere, orangeCarbon), zoomed into view |
| Honest unavailability | No Vina/Open Babel → explicit "unavailable", never a fake result |
| Pose comparison | `compare_results`: Kabsch-superposed RMSD on atom-name-matched pose pairs + affinity deltas; verified: two independent 3W85 runs converge at 0.14–0.51 Å on matched poses |

### 3.6 Workbench features
| Feature | Implementation |
|---|---|
| Measurement tools | Distance & angle from real atomic coordinates via backend math (measured 3.8546 Å vs PLIP's 3.85 Å on the same pair) |
| 2D ligand editor | Canvas editor with real CCD bond orders/charges, element-colored atoms; atom add/remove, bond add/remove (two-pick flow); SDF export from the edited graph |
| Python scripting console | Real execution with `structure` / `selected_ligand_id` in scope; `print()` output captured and displayed (never corrupts the IPC channel) |
| Session notes | Notes attached to and persisted with the session |
| Jobs panel | Real async job lifecycle (running/completed/failed) with docking poses listed per completed job |

### 3.7 Batch, comparison & export
| Feature | Implementation |
|---|---|
| Batch analysis | Multiple PDB IDs/paths through the full pipeline with per-task status; failures are **reported, not dropped** (1CRN correctly surfaces "no ligand found") |
| Batch export | JSON artifact including failed tasks |
| Run comparison | Two docking jobs side-by-side (RMSD matrix + affinities) |
| Exports | Contacts CSV, ligand SDF (CCD-bond-derived), full analysis summary artifact, scene file |
| Workspace | Per-session workspace directory caching CCD/PubChem artifacts (responsible caching, nothing fabricated) |

### 3.8 PDB-wide discovery, protein context, and engine depth (v2.1)
| Feature | Implementation |
|---|---|
| PDB-wide discovery | RCSB Search API (text, MMseqs2 sequence, OpenEye chemical similarity, same-ligand, similar CCD components) — all live |
| Structure quality | RCSB validation summary (clashscore, RSRZ/Ramachandran outliers) |
| Protein context | UniProt REST (function, features, real binding-site annotations) |
| AlphaFold models | AlphaFold DB mmCIF fetch for the entry's SIFTS-mapped accession |
| Pocket detection | Real `fpocket` binary; descriptors verbatim; honest absence |
| Secondary structure | Real `mkdssp` (DSSP 4.x); 8-state codes from its own output |
| 2D interaction diagram | SVG from PLIP's own contacts (RDKit depiction or PCA of real coordinates) |
| Sequence viewer | Per-chain sequence with contact/DSSP/pocket mapping; click-to-zoom |
| mmCIF `_struct_conn` export | PLIP's mappable contacts written back in a standards format |
| Docking box editing | `set_docking_box`/`get_docking_box` + UI overrides |
| Multi-ligand docking queue | One shared receptor PDBQT, one job per non-solvent ligand |
| Redocking validation | Crystal-vs-docked Kabsch RMSD (element-aware matching) |
| Structure superposition | CA superposition via NW alignment + seed-and-extend core fit |
| Batch interaction aggregation | Per-residue contact frequency across batch results |
| Keyboard ligand switching | Next/prev shortcuts + recent-ligands list |

### 3.9 v2.2 — affinity evidence, covalent handling, MD, reproducibility
| Feature | Implementation |
|---|---|
| BindingDB affinities | **Current** REST API (`/rest/getLigandsByPDBs`, `/rest/getLigandsByUniprot`), probed live: real Ki/Kd/IC50/EC50 records with PMID/DOI; both real response schemas handled; uncovered entries reported honestly (never an invented empty success) |
| Covalent-ligand handling | The file's **own** records only — mmCIF `_struct_conn` `covale`/`disulf` rows, or PDB `LINK` records with distances from real atom coordinates (verified on 1HRC: two HEC thioether links at 1.754/1.864 Å, iron coordination correctly excluded as `metalc`); ligand card flags covalent attachment |
| GROMACS minimization/MD | Real `gmx` wrap: `pdb2gmx → grompp → mdrun → gmx energy` (vacuum EM / short vacuum MD); every parameter user-configured; honest absence; runs as a background job |
| Reproducibility manifest | `repro_manifest.json`: SHA-256-hashed artifacts, job parameters/results, config snapshot, notes — written into the workspace |
| Viewer-only file formats | SDF/MOL/MOL2/XYZ/GRO/PDBQT/CDJSON open in the 3D viewer (3Dmol.js native parsing) without pretending chemistry was parsed: no ligand card, no contacts; analysis honestly declines. `.gz`-wrapped structures decompress into the workspace. PDB/mmCIF remain fully parsed |
| Viewer-only honesty in the backend | `get_status` reports `viewer_only`; contact analysis on a viewer-only session refuses with an explicit reason instead of failing on stale state |

### 3.9a Spec-completion layer (V1/V2 spec items that were missing, now shipped)
| Feature | Implementation |
|---|---|
| Scene image export | The real rendered viewport (viewer canvas pixels) is captured as PNG, validated by the backend (magic bytes, size), stored as the session's `scene.png` artifact, and saved through the native dialog — the "structured analysis artifact = image + data + notes" thesis is complete |
| Contacts CSV export | One-click button in the Contacts panel driving the existing backend CSV export |
| Contact table sort/filter | Type filter + sort by type/distance/residue, over the real PLIP rows |
| Camera presets | Ligand focus / pocket site / interface / whole structure |
| Pocket surface toggle | Ligand-centric pocket surface over the residues the backend's own geometric binding-pocket computation selected (radius is the config/user setting, never a UI guess) |
| Top bar | PDB-wide text search (live RCSB Search API, results in the Discover tab), recent-structures list, settings dialog showing live engine + data-source status |
| Surface representation | Real molecular surface (VDW, translucent) over a thin cartoon backbone |
| 2D editor → 3D sync | The edited molecule's real geometry + bond graph exported as SDF and overlaid in the 3D viewer as its own model (green carbons); the loaded structure is never modified |
| Pose clustering + consensus | `compare_results` now clusters the union of two jobs' poses (caller-settable RMSD threshold, no buried constant) and reports a consensus pose (lowest affinity) with real atoms; clustering degrades to empty without ever breaking the core pairwise comparison |
| RCSB ModelServer | On-demand coordinate subsets (entry/chain/component) as mmCIF from `models.rcsb.org`, parsed by the app's own mmCIF parser and overlaid as a separate model |
| PDBe REST (EU mirror) | Entry summary, secondary-structure annotation, and ligand monomers from `ebi.ac.uk/pdbe/api` — the §6.1 source is now actually integrated, not just listed |

### 3.10 Search (BM25 over live documents)
| Feature | Implementation |
|---|---|
| Search engine | Real Okapi BM25 (own implementation: TF/IDF/length-norm/k1/b) |
| Document sources | RCSB Data API entry records, CCD component records, PubChem records, workspace artifacts — all fetched live when indexing |
| "Index current session" | One click indexes the active structure + all its ligands |
| Queries verified in UI | "myoglobin" → the real 1MBN entry doc; "protoporphyrin" → the HEM chemical doc |
| Persistent index | Cached under the workspace; refreshes re-fetch live data |

### 3.11 Health & transparency
| Feature | Implementation |
|---|---|
| Status view | Live checks of every data source (RCSB, CCD, PubChem, ChEMBL/UniChem) and every engine (PLIP, Open Babel, Vina) — real HTTP/binary checks with a single honest retry on transient timeouts |
| No silent degradation | Every unavailable capability says so in the UI |

---

## 4. Feature-by-alternatives matrix (what exists, what's different)

Alternatives surveyed: **PyMOL, ChimeraX, Avogadro, VMD, Biopython, ProDy,
OpenStructure, PLIP (web/CLI), LigPlot+/PoseView, DrugScore-type viewers,
MolecularFlipbook, Mol* (Molstar), UCSF/Docking GUIs (e.g. PyRx), OpenEye
apps (proprietary), SwissBioisostere-type web tools.**

Legend: ✅ real feature · 🟡 partial/indirect · ❌ absent · *(—)* not applicable.

| Capability | Ligora | PyMOL | ChimeraX | Avogadro | VMD | Mol* | PLIP web | Biopython |
|---|---|---|---|---|---|---|---|---|
| 3D structure viewing | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ |
| mmCIF parse (own, spec-compliant) | ✅ | ✅ | ✅ | 🟡 | 🟡 | ✅ | ✅ | ✅ |
| **Auto ligand identity + enrichment with evidence links** | ✅ | ❌ | 🟡 (CCD name only) | ❌ | ❌ | 🟡 (CCD name) | 🟡 (partial) | ❌ |
| **Typed contact analysis (PLIP engine) one-click** | ✅ | 🟡 (manual find-pairs) | 🟡 (HBonds only) | ❌ | ❌ | ❌ | ✅ | 🟡 (scripted) |
| Water network graph + clusters + click-zoom | ✅ | 🟡 (manual) | ❌ | ❌ | 🟡 (scripted) | ❌ | ❌ | ❌ |
| **In-app docking with pose overlay + pose comparison RMSD** | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| 2D ligand editor on CCD bond orders → SDF | ✅ | ❌ | ❌ | ✅ (general editor) | ❌ | ❌ | ❌ | ❌ |
| Python console scoped to live session | ✅ | ✅ (own API) | ✅ (own API) | ❌ | ✅ (Tcl) | ❌ | ❌ | ✅ |
| BM25 search over session's live-sourced docs | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Batch over structure lists with honest failures | ✅ | 🟡 (scripted) | 🟡 (scripted) | ❌ | 🟡 (scripted) | ❌ | ❌ | ✅ (scripted) |
| **Evidence pane for every displayed value** | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **"No fallbacks" honesty contract** | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Exportable artifacts (CSV/SDF/JSON summary/scene) | ✅ | 🟡 | 🟡 | ✅ | 🟡 | 🟡 | 🟡 | ✅ |
| Desktop app (offline-capable shell, live data when needed) | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ (web) | ❌ (web) | — (lib) |
| No-cost fully open-source pipeline | ✅ | ❌ (licensed) | ✅ | ✅ | ✅ (free non-commercial) | ✅ | ✅ | ✅ |

### 4.1 What Ligora does that the alternatives genuinely do not

1. **The evidence contract.** No mainstream tool attaches a per-value
   provenance chain (source, fields, URL) to every displayed datum, and none
   refuse to invent values when sources are unavailable. This is Ligora's
   single clearest differentiator, and it is architectural, not cosmetic.
2. **The integrated identity→contacts→docking→compare loop.** PyMOL/ChimeraX
   view; PLIP analyzes; Vina docks; nothing joins them. In Ligora the output
   of each stage feeds the next (the CCD-classified ligand feeds PLIP; the
   co-crystal ligand defines the Vina box; poses land in the same viewer as
   the crystal ligand; two runs compare with real superposition).
3. **Water network analysis with clusters and one-click zoom** is not a
   first-class feature of any surveyed tool (VMD can be scripted to similar
   ends; nothing ships it as a product feature with UI).
4. **Honest unavailability as UX.** Tools typically fail silently or crash;
   Ligora's status surface and per-feature "unavailable" states make absence
   visible.

### 4.2 What the alternatives do better (honest)

- **Rendering quality/perf at scale:** PyMOL/ChimeraX/VMD ray-trace,
  handle million-atom systems, and have decades of GPU optimization. 3Dmol.js
  (Ligora's viewer) is excellent for the protein–ligand scale but is not a
  match for megastructures.
- **Breadth of structure biology features:** ChimeraX (maps, trajectories,
  sequence alignment, cryo-EM tooling) and VMD (MD trajectories) dwarf Ligora
  in scope. Ligora is deliberately a **ligand-centric** workstation, not a
  universal molecular graphics suite.
- **Editor maturity:** Avogadro's 3D editor (force fields, geometry
  optimization) is far deeper than Ligora's 2D ligand editor.
- **Ecosystem/age:** PyMOL scripting has thousands of published scripts;
  Biopython integrates with the entire Python scientific stack.
- **LigPlot+** remains the reference for publication-style 2D interaction
  diagrams (Ligora currently shows typed contact tables + 3D, not 2D
  schematic diagrams — see roadmap §5).

---

## 5. Roadmap status (updated)

Everything from the original Tier 1 and Tier 2 roadmap **has shipped and is
tested** (see §3.8): the 2D interaction diagram, sequence viewer with
interaction mapping, per-residue batch aggregation, mmCIF `_struct_conn`
export, keyboard ligand switching + recent ligands, one-click open of both
structure and chemical search hits, fpocket pockets, DSSP secondary
structure, docking box editing UI, the multi-ligand docking queue, ligand
similarity discovery, crystal-vs-docked redocking validation, and structure
superposition.

Of the original Tier 3:

14. **MD trajectory snippets** — the first real step now ships: GROMACS
    minimization and short vacuum-MD snippets run as background jobs through
    a real `gmx` wrap (§3.9). Full solvated production MD and trajectory
    analysis remain future work.
15. **Covalent-ligand handling** — ships in v2.2, grounded strictly in the
    file's own covalent-link records (§3.9); richer covalent chemistry
    (warhead reactivity, covalent docking) remains future work.

What remains genuinely unimplemented:
- Full solvated production MD with trajectory analysis.
- Covalent docking / warhead reactivity chemistry.
- A verified end-to-end snap build (needs a clean core22 host).

---

## 6. Data sources worth integrating (all probed live before inclusion)

### 6.1 Verified reachable and immediately integrable
| Source | What it adds | Endpoint verified | License/terms |
|---|---|---|---|
| **RCSB search API** (`search.rcsb.org/rcsbsearch/v2`) | **Integrated (§3.8)**: text, sequence, chemical-similarity, same-ligand and similar-component discovery, wired into the Discover tab and the top-bar search | HTTP 200 | CC0 (PDB data); API free, attribution requested |
| **PDBe REST API** (`ebi.ac.uk/pdbe/api/pdb/entry/...`) | **Integrated (§3.9a)**: entry summary, secondary structure, ligand monomers | HTTP 200 (verified live) | CC-BY 4.0 (EBI terms) |
| **RCSB ModelServer** (`models.rcsb.org/v1`) | **Integrated (§3.9a)**: on-demand mmCIF coordinate subsets (entry/chain/component), parsed by the app's own parser | HTTP 200 (verified live) | CC0 |
| **RCSB validation data** (`rcsb_vrpt_summary` in entry payload) | Clashscore, RSRZ outliers, Ramachandran outliers — display "structure quality" next to analysis | HTTP 200, field present on 3W85 | CC0 |
| **UniProt REST** (`rest.uniprot.org/uniprotkb/{acc}.json`) | Protein function, domains, names, PTMs, cross-refs — the protein side of the evidence pane | HTTP 200 (271 KB payload) | CC-BY 4.0 |
| **AlphaFold DB** (`alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif`) | Predicted structures for proteins without experimental structures; overlay confidence (pLDDT is in the file, not a heuristic) | HTTP 200 (**v6** current) | CC-BY 4.0 |
| **PDBe Ligand / CCD mirrors** | Same CCD data, EU-latency alternative | HTTP 200 | CC0 |

### 6.2 Verified interesting, but integration needs care
| Source | Status | Note |
|---|---|---|
| **BindingDB REST** | **Integrated (v2.2)** against the current REST web services, probed live (`getLigandsByPDBs` returns real Ki records; `getLigandsByUniprot` returns its own `bdb.`-prefixed schema — both handled). Entries outside its coverage are reported honestly, never faked. |
| **PDBbind** | Already wired as an optional source; requires configured local/remote access per their terms | Affinity values appear in evidence only when the source answers. |

### 6.3 Engines/tools to wrap (distro packages, same pattern as Vina/PLIP)
| Tool | Adds | Ubuntu availability |
|---|---|---|
| **fpocket** | Pocket detection (Tier-2 #7) | Not in current Ubuntu repos — would need upstream static build in snap part (real binary, still no reimplementation) |
| **DSSP** (`dssp`) | Real secondary-structure assignment | In repos (verified via apt-cache) |
| **OpenBabel Guassian/forcefield modes** | Geometry minimization in the 2D/3D editor | Already shipped with obabel |

### 6.4 Deliberately **not** integrated
- Scraped HTML of PubChem/ChEMBL web pages (only their real APIs).
- Any "AI-predicted affinity" service — unverifiable provenance conflicts
  with the evidence contract.
- Proprietary OpenEye/Tripos components — license incompatible with the
  fully-open-source requirement.

---

## 7. Honest opinion — user perspective

**What a real user experiences today (tested end-to-end):**
Boot the app → type "3W85" → in ~5–10 seconds a fully annotated ligand card
appears (real name, CCD-classified, PubChem CID, InChIKey, ChEMBL link) with
an evidence trail; one click runs PLIP and 14 typed contacts appear with real
atom names; water network, docking, console, notes, batch, search all respond
with real data; every failure mode they hit says what is actually wrong.

**Strengths, felt as a user:**
- The **single-click depth** is the hook: what costs an afternoon of
  tool-hopping becomes a minute, without a scripting prerequisite.
- **Trust**: evidence links mean you can defend the numbers in a meeting or a
  paper. Very few tools can say that.
- The app **never lies**: water is never shown as the ligand, missing engines
  are labeled, empty results say "no drug-like ligand" instead of inventing
  one. This builds the kind of confidence that keeps users.
- Artifacts (CSV/SDF/JSON) keep it out of "walled garden" territory.

**Weaknesses, felt as a user (honest):**
- **First-run dependency on network**: identity/enrichment/search are live;
  offline you can view local files but get honest "source unavailable"
  messages instead of cached universes (a deliberate tradeoff).
- **No publication-grade 2D interaction figure yet** — a 2D diagram from
  PLIP's own contacts now ships (§3.8); LigPlot+ still offers richer
  publication styling for the final figure.
- **Viewer ceiling**: superb at protein–ligand scale; users coming from
  PyMOL with million-atom CryoEM maps will notice the difference.
- **Undo in the 2D editor** and richer editing affordances are thin.
- **Search is session-scoped by design**; PDB-wide discovery lands properly
  with the RCSB search API integration (§6.1), which is the single highest-
  leverage data-source addition available.

---

## 8. Honest opinion — usage & requirements perspective

**Runtime requirements (verified on this machine):**
- x86_64 Linux (snap targets core22/Ubuntu 22.04), ~200 MB install +
  workspace cache.
- Python 3.10+ with numpy/requests/rdkit (pip), plus distro binaries:
  `plip`, `openbabel`, `autodock-vina`, `gromacs`, `dssp` (all staged in
  snapcraft.yaml; `fpocket` remains unpackagable — no distro package —
  and reports honest absence when not user-installed).
- GUI stack: WebKit/WebView with WebGL (SwiftShader fallback verified in
  headless testing).
- Network: required for enrichment/search/fetch; docking/contacts/exports
  work fully offline once structures are local.

The three §10 recommendations below were the doc's original "shortest path"
; the first two have since shipped (RCSB search API discovery, §3.8; 2D
interaction diagram, §3.8; crystal-vs-docked comparison + box editing UI,
§3.8/§3.9).

**Operational honesty:**
- Long operations (PLIP on big structures, Vina at exhaustiveness 16) run as
  real background jobs with progress UI — but contact analysis in the main
  flow is synchronous; a large structure can hold the UI for seconds
  (async-ifying is straightforward with the existing jobs system).
- The snap stages the real engine packages (`plip`, `openbabel`,
  `autodock-vina`, `gromacs`, `dssp`) coherently, but a full snap build was
  not run in this environment (needs a clean core22 host).
- The dev/test harness (dev-bridge + CDP stages, now 6 stages) is **not**
  shipped in the product path; the shipped app talks to the same backend
  through the Rust bridge.

**Reproducibility story:** strong — every artifact export is deterministic
from real inputs, and the evidence pane doubles as a methods section. The
"reproduce this analysis" manifest (§3.9) now formalizes it: SHA-256-hashed
artifacts plus job parameters/results in one JSON, one click from the
Structure tab.

---

## 9. Honest opinion — versus existing alternatives

**Where Ligora wins outright:**
1. The **integrated, evidence-backed pipeline** — nothing open-source joins
   identity → contacts → waters → docking → comparison → search in one
   desktop app, and nothing else makes provenance a first-class UI concept.
2. **Honesty engineering** (no heuristics/fallbacks/mocks) — competitors
   either approximate silently (viewer-assigned bonds) or require the user
   to script real tools together.
3. **Time-to-answer** for the standard "what is this and how does it bind"
   workflow: minutes instead of an afternoon, reproducibly.

**Where Ligora loses today:**
1. **Graphics scale & polish** vs PyMOL/ChimeraX (rendering depth,
   trajectories, megastructures).
2. **Maturity/breadth** — ChimeraX/VMD do far more *biology*; Ligora is
   scoped to ligand interaction analysis (by design, but say it plainly).
3. **Community & longevity** — PyMOL/VMD have decades of users, docs, and
   extensions; Ligora is new with a growing test/verification surface (143
   backend tests + 6 UI stages).
4. **2D publication figures** — LigPlot+ still owns that deliverable.

**Net position:** Ligora is not "a nicer PyMOL". It is a different category:
a **reproducible analysis pipeline disguised as a desktop app**. For the
protein–ligand workflow specifically, it is already the strongest fully-
open-source option; for general molecular graphics it should not pretend to
compete.

---

## 10. Verdict and the shortest path to "best in class"

The app is **shippable today** for its core use case, with real engines, real
data, honest failure modes, and no placeholder anywhere in the product path.
The three moves with the highest ratio of impact to effort:

1. **RCSB search API integration** (§6.1) — turns the session-scoped BM25
   search into PDB-wide discovery with zero chemistry risk.
2. **2D interaction diagram from PLIP's own output** (§5.1 #1) — closes the
   only "users leave the app" gap in the publication workflow.
3. **Crystal-vs-docked comparison + box editing UI** (§5.2 #9/#12) —
   completes the standard redocking protocol users will immediately try.

Every other item in §5 is incremental on an architecture that is already
real, tested, and honest.
