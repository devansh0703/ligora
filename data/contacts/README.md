# Contacts Data Directory

Contact analysis writes `contacts.csv` in the session workspace
(`~/.ligora/workspaces/<session>/contacts.csv`) on every analysis run.

## Provenance

- Rows come from **PLIP's XML report** executed on the complex (protein +
  selected ligand) rendered to PDB by Ligora.
- Contact types are PLIP's own classification (hydrogen bond, hydrophobic,
  salt bridge, π-stacking, halogen bond, metal complex, water bridge).
- Participating atoms are resolved from PLIP's serial fields
  (`donoridx/acceptoridx`, `ligcarbonidx/protcarbonidx`, index lists) against
  the exact PDB file that was analyzed; for sections that report group
  coordinates instead of serials (salt bridges, π interactions), the atom is
  the nearest atom of the interacting residue to PLIP's reported coordinate —
  geometry from the file, not a heuristic.
- Distances/angles are PLIP's values.

## CSV schema

```
id, ligand_atom, ligand_res, ligand_res_id, ligand_chain,
protein_atom, protein_res, protein_res_id, protein_chain,
distance, contact_type, angle, water_mediated, description
```
