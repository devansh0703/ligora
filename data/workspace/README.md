# Workspace Data Directory

User data lives under `~/.ligora/` (override with `LIGORA_BASE_DIR`):

```
~/.ligora/
├── workspaces/<session-id>/
│   ├── structure.mmcif        # the exact file fetched/opened
│   ├── contacts.csv           # PLIP results (see ../contacts)
│   ├── analysis_summary.json  # saved artifact (summary + contacts + evidence)
│   ├── ligand_<ID>.sdf        # exported ligand with real CCD-derived bonds
│   ├── edited_ligand.sdf      # 2D editor export (RDKit-validated)
│   └── job_vina_<ts>/         # docking jobs (see ../docking)
├── batch_<batch-id>/          # batch analysis exports
└── ccd_cache/                 # CCD component CIF cache (fetched documents,
                               # keyed by component ID; TTL-limited)
```

## What is cached vs fetched

- **Cached**: CCD component files (they are immutable snapshots keyed by
  component ID) — this is transport caching of real documents, not
  fabrication.
- **Always fetched live**: entry metadata, PubChem properties, ChEMBL
  records, UniChem mappings (they change over time).
- **Never stored**: invented values. If a field has no source answer, it is
  absent in files and shown as "not available" in the UI.
