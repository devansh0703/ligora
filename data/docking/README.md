# Docking Data Directory

Each docking run writes a self-contained, reproducible job directory under
the session workspace (`~/.ligora/workspaces/<session>/job_vina_<ts>/`):

```
job_vina_<ts>/
├── output/
│   ├── receptor.pdb        # polymer chains, written from the parsed structure
│   ├── receptor.pdbqt      # Open Babel conversion (-xr rigid)
│   ├── ligand.sdf          # RDKit 3D SDF (CCD bonds + charges)
│   ├── ligand.pdbqt        # Open Babel conversion (--gen3d -h)
│   ├── out.pdbqt           # Vina poses (REMARK VINA RESULT per mode)
│   └── vina.log            # full Vina stdout (scoring params, timings)
└── job.json                # parameters, box, engine version, status
```

## Provenance

- **Docking box**: computed from the co-crystallized ligand's own extent
  (centroid + 8 Å padding), overridable per-run via `set_docking_box` or
  `run_docking {box_center, box_size}`. The box is geometry from the user's
  structure — never a fixed default box.
- **Affinities**: produced by the `vina` binary (version recorded in
  `vina.log`). Ligora never estimates or fabricates affinities; when Vina or
  Open Babel is absent the command fails with an install hint.
- **Pose atoms**: parsed from `out.pdbqt` fixed columns (PDBQT spec),
  including per-atom AD4 types → elements, so poses render in 3Dmol.js.

## Reproducing a run

The exact command line is reconstructed in `vina.log`; parameters
(exhaustiveness, num_modes, energy_range, seed) come from the run request
payload stored in `job.json`.
