# Molara — Product Spec (v1.1)

## Intent
Molara is an industrial-grade open-source genomics workstation UI for X86_64 desktops. It does **not** try to re-implement 공유 open-source analysis engines. It makes existing analyses and datasets easier to reach, easier to compare, easier to reproduce, and easier to export. The UI is the product; the backend is a typed integration layer over live open data sources plus local computation where it is genuinely useful.

## Rules (binding)
1. **Everything visible must work end-to-end.** If a feature exists in code, it gets a real UI path. If a UI path exists, it must run against real data, not sample stubs.
2. **Do not reimplement shared open-source tools.** Use them, wrap them, or link out to them. Molara is a UI and integration layer, not a new engine library.
3. **Data sources only, not hardcoded chemistry/biology intuition.** Classification, thresholds, weights, formulas, mapping rules, significance defaults, and similar numeric or categorical knowledge must come from datasets, API responses, or explicit user settings — not from buried constants in app code.
4. **Prefer authoritative live data.** Use RCSB, PubChem, ChEMBL, PDBBind, genome annotations, GWAS catalog, reference variant sources, and other open endpoints where they exist. Cache responsibly, but do not invent values when a source can answer.
5. **Make rare and missing analyses usable.** If a dataset or analysis type already exists in the code but has no UI, finish the UI. If a useful analysis does not exist, add the smallest sensible layer that connects real data to the user.
6. **X86_64 desktop first, snap-ready.** Target desktop deployment and make packaging real, including data, desktop integration, and confinement semantics. No throwaway build scripts.
7. **No placeholder data in the product.** Sample data may exist only in tests, demos, or documentation. The app itself should operate on user data or live sources.

## Problem framing
Existing molecular-biology GUIs tend to be either:
- narrow viewers that do not move beyond one structure or one track, or
- engine CLIs dressed up lightly, with weak dataset integration and weak comparison/export workflows.

Molara targets the gap between those: a desktop UI that lets a user inspect a structure or variant context, bring in multiple open datasets, run or attach analyses without hand-editing files, and keep the result set comparable and exportable.

## Core jobs
1. Load a structure or genomic context and see it in a real viewer.
2. Inspect ligands, contacts, evidence, and derived properties in one place.
3. Enrich identity and metadata from open sources without manual copy-paste.
4. Run or attach analysis workflows through the UI, not only through scripts.
5. Compare results across structures, ligands, runs, or cohorts where applicable.
6. Export artifacts that are actually usable outside Molara.

## Feature matrix

### V1 — Core workstation
- Open local PDB/mmCIF and open PDB IDs from RCSB.
- 3D scene with representation and coloring controls that actually map to scene state.
- Ligand card grounded in enriched ligand data.
- Contacts panel with typed contact rows and export.
- Evidence pane connected to enrichment sources.
- Job panel for running tasks and viewing status.
- Measurement tools with real coordinate math.
- Notes attached to session/analysis.
- Export paths for contacts, ligand, summary, and scene artifacts.
- Health/status view for data sources and engines.

### V2 — Analysis depth
- Batch analysis for multiple structures/IDs with progress and result set export.
- 2D ligand editor with coordinate-aware editing and sync path to 3D.
- Water network analysis with clustering and contact mapping.
- Scripting console for custom workflows with access to Molara API objects.
- Result comparison with pose/aggregate comparison and clustering.

### V3+ — Industrial expansibility
- Additional open-data enrichments as datasets become usable.
- Additional comparison and cohort views as analytics are added.
- Packaging, metadata, and CI for repeatable desktop shipping.

## Data sources map

### Structural biology
- RCSB PDB / mmCIF download and metadata endpoints.
- RCSB chemical component dictionary for identity and classification.
- PubChem for compound identity, properties, and cross-references.
- ChEMBL for bioactivity and target context where available.
- PDBBind for affinity context where available.

### Genomics / molecular biology (where applicable to workflow)
- Reference annotations and public variant/gene data via open endpoints or downloadable open datasets.
- Public catalog-style sources for variant/disease associations when the workflow needs them.
- Local open datasets supported where live lookups are not appropriate.

## UI surfaces to complete
- Scene viewer controls that mutate scene representation and coloring.
- Ligand card that reflects enrichment state.
- Contacts table that is filterable/viewable and exportable.
- Evidence panel that shows source, field, value, and link.
- Job panel that shows running/complete/failed work.
- Measurement UI with selection and results.
- Notes/editor affordance attached to session.
- Batch UI for adding sources and launching batch runs.
- 2D editor UI synchronized with 3D where applicable.
- Water network UI for clusters and contacts.
- Scripting console UI with history and variable view.
- Comparison UI for result sets.

## Non-goals
- Re-implement Vina, gnina, PLIP, OpenMM, PySCF, RDKit-class cheminformatics, or other established open tools as primary Molara logic.
- Ship fabricated sample data as if it were real product data.
- Harden numeric/biological rules inside the app when a dataset or source can supply them.

## Acceptance
- A user can open a real structure from file or PDB ID and see it in the viewer.
- A user can select a ligand and see enriched identity and evidence from live sources.
- A user can run contact analysis and export contacts.
- A user can run V2 workflows from the UI with real progress and results.
- Packaging builds and runs on the intended desktop target.
