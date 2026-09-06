# Ligora Data Sources

Every fact Ligora displays comes from one of these sources. Nothing is
invented in application code: classifications, thresholds and chemistry come
from these datasets and APIs; when a source cannot answer, Ligora says so.

Metadata below was verified against each provider's published terms in
September 2026 (verification notes per source).

## Structural data

### RCSB PDB (wwPDB archive)
- **What**: experimental 3D structures (mmCIF/PDB files), entry metadata,
  primary citations, and the Chemical Component Dictionary (CCD) with
  per-component atoms, bonds, ideal coordinates, descriptors and the
  `pdbx_type` classification used for solvent/ion/organic decisions.
- **Endpoints used**:
  - `https://files.rcsb.org/download/<ID>.cif` (structure files)
  - `https://data.rcsb.org/rest/v1/core/entry/<ID>` (metadata + citation)
  - `https://data.rcsb.org/rest/v1/core/chemcomp/<ID>` (CCD component)
  - `https://files.rcsb.org/ligands/view/<ID>.cif` (CCD atom/bond loops)
  - `https://search.rcsb.org/rcsbsearch/v2/query` (text search API)
- **License**: PDB archive data are CC0 1.0 (wwPDB usage policy,
  https://www.wwpdb.org/about/usage-policies).
- **Citation**: available per-entry via `rcsb_primary_citation`
  (`data.rcsb.org/rest/v1/core/entry/<ID>`).
- **Verified**: 2026-09; endpoints exercised by the test suite with live
  requests (3W85, 1CRN, 1UBQ, ...).

## Chemical identity & properties

### PubChem (NCBI)
- **What**: compound resolution by name/CID, IUPAC name, connectivity SMILES,
  molecular formula, exact mass, TPSA, XLogP, synonyms.
- **Endpoints used**:
  - `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{cid|name}/...`
  - `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fast/similarity/...`
- **License**: PubChem is an open archive; data are public domain
  (CC0-equivalent; NCBI asks users to acknowledge PubChem,
  https://pubchem.ncbi.nlm.nih.gov/docs/data-sources).
- **Note**: `ConnectivitySMILES` is the current property name; the legacy
  `CanonicalSMILES` field is deprecated and no longer returned.
- **Verified**: 2026-09; live in `tests/test_enrichment.py`.

## Bioactivity & targets

### ChEMBL (EMBL-EBI)
- **What**: bioactivity records (assay type, target, pChEMBL value) and
  compound metadata by ChEMBL ID.
- **Endpoints used**:
  - `https://www.ebi.ac.uk/chembl/api/data/molecule/<ID>.json`
  - `.../activity.json?molecule_chembl_id=<ID>`
  - `.../chembl_id_lookup/search?q=<term>` (text search)
- **License**: CC-BY-SA 3.0 Unported
  (https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions).
- **Verified**: 2026-09; live in `tests/test_enrichment.py`.

### UniChem (EMBL-EBI)
- **What**: cross-references between chemical identifiers (InChIKey →
  ChEMBL ID and other source databases).
- **Endpoint used**: `https://www.ebi.ac.uk/unichem/rest/inchikey/<key>`
- **License**: EMBL-EBI Terms of Use impose no additional restrictions on
  contributed data (https://www.ebi.ac.uk/about/terms-of-use/); source
  databases retain their own licenses.
- **Verified**: 2026-09; live in `tests/test_enrichment.py`.

## Analysis engines (executed locally; never re-implemented)

### PLIP — Protein-Ligand Interaction Profiler
- **Role**: all contact classification (H-bonds, hydrophobic, salt bridges,
  π-stacking, halogen, metal, water bridges). Ligora parses PLIP's XML
  report; it never invents contact types.
- **License**: GPL-2.0-or-later (pharmai/plip).
- **Citation**: Adasme et al., Nucleic Acids Res. 2021 (PLIP 2021);
  Salentin et al., Nucleic Acids Res. 2015 (PLIP 2015).

### AutoDock Vina
- **Role**: molecular docking. Input prep via Open Babel; box from the
  co-crystallized ligand unless user-specified.
- **License**: Apache-2.0 (ccsb-scripps/AutoDock-Vina).
- **Citation**: Eberhardt et al., J. Cheminform. 2021 (Vina 1.2);
  Trott & Olson, J. Comput. Chem. 2010 (Vina 1.0).

### Open Babel
- **Role**: PDBQT preparation, geometry cleanup/minimization.
- **License**: GPL-2.0 (openbabel).
- **Citation**: O'Boyle et al., J. Cheminform. 2011.

### RDKit
- **Role**: bond perception (`rdDetermineBonds`), SDF generation,
  fingerprints/similarity, descriptors, sanitized 2D editing.
- **License**: BSD-3-Clause (rdkit).
- **Citation**: RDKit, https://www.rdkit.org (open-source cheminformatics).

### 3Dmol.js
- **Role**: WebGL molecular viewer in the UI.
- **License**: BSD-3-Clause (3dmol/3Dmol.js).
- **Citation**: Mitchell et al., Bioinformatics 2021.

## Provenance rules (binding, from RULES.md)

1. Classification (solvent/ion/organic) comes from CCD `pdbx_type`, never
   from name lists in code.
2. Bond orders, formal charges and elements come from CCD
   `_chem_comp_atom`/`_chem_comp_bond`, never guessed from geometry alone.
3. Contact types come from PLIP; affinities from Vina; both report absence
   honestly when the engine is unavailable.
4. When a source cannot answer (e.g. PDBBind has no open API), the UI shows
   "not available" — no synthesized values.
