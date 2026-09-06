# Ligora

**Molecular Analysis Workstation**

Ligora is an open-source desktop application for protein-ligand interaction analysis, molecular docking, and computational chemistry workflows. It combines a 3D molecular viewer with analysis tools and integration with external open-source engines.

## Features

### V1 (Core)
- **Structure Visualization**: 3D viewer for PDB/mmCIF structures with multiple representations (cartoon, stick, sphere, surface)
- **Ligand Detection**: Automatic detection and classification of ligands in structures
- **Contact Analysis**: Protein-ligand interaction analysis using PLIP and geometric methods
  - Hydrogen bonds
  - Hydrophobic contacts
  - Pi-stacking
  - Salt bridges
  - Halogen bonds
  - Metal coordination
  - Water-mediated contacts
- **Ligand Identity Resolution**: Automatic lookup in PubChem, ChEMBL, and PDBBind
- **Live Data Enrichment**: Fetch metadata from RCSB PDB, PubChem, ChEMBL, PDBBind
- **Measurement Tools**: Distance and angle measurements between atoms
- **Export**: Analysis summaries (JSON), contacts tables (CSV), ligand structures (SDF)

### V2 (Engines)
- **Molecular Docking**: AutoDock Vina integration for protein-ligand docking
- **Geometry Cleanup**: Local ligand geometry optimization
- **Batch Analysis**: Process multiple structures
- **2D Ligand Editor**: Sync between 2D and 3D views
- **Water Network Analysis**: Solvent molecule analysis

### Future (V3+)
- **MD Simulation**: OpenMM integration for molecular dynamics
- **QM Calculations**: PySCF integration for quantum chemistry

### Rules
- Everything visible in the app runs end-to-end against real data or user input.
- Existing open-source engines are used, wrapped, or linked out; they are not reimplemented.
- Numeric and categorical chemical/biological values come from data sources (RCSB, PubChem, ChEMBL, PDBBind, CCD) or user settings, not from hardcoded app constants.
- Atomic weights, formulas, classification hints, thresholds, and similar knowledge are not embedded in app code.
- There are no placeholder data or sample values in the product itself; sample data exists only in tests, demos, or docs.
- Engines are integrated under their own licenses; bundling is only done after licensing and redistribution terms are verified.

Full rules are in `RULES.md`.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                      Tauri Frontend (Rust + Web)            │
│  ┌─────────────────────────────────────────────────────────┐ │
│  │  3D Scene Viewer (Three.js/3Dmol.js)                   │ │
│  │  Ligand Card        │  Contacts Table                   │ │
│  │  Evidence Pane      │  Jobs Panel                      │ │
│  │  Notes              │  Engine Controls                 │ │
│  └─────────────────────────────────────────────────────────┘ │
│                              │                               │
│                    IPC (JSON-RPC)                          │
│                              ▼                               │
├─────────────────────────────────────────────────────────────┤
│                    Python Backend                           │
│  ┌─────────────┬─────────────┬─────────────┬─────────────┐ │
│  │  Parser     │  Ligand     │  Contacts   │ Cheminform  │ │
│  │  (mmCIF/   │  Resolver   │  Analyzer   │ (fingerprints│
│  │   PDB)     │  (PubChem/  │  (PLIP/     │ /similarity)│
│  │            │   ChEMBL)   │   geometry) │              │ │
│  └─────────────┴─────────────┴─────────────┴─────────────┘ │
│  ┌─────────────┬─────────────┬─────────────┬─────────────┐ │
│  │  Enrichment │  Engines    │  Jobs       │  Export     │ │
│  │  (RCSB/     │  (Vina/     │  Manager    │  (JSON/CSV/ │
│  │   PubChem/  │   gnina/    │             │   SDF)      │ │
│  │   ChEMBL/   │   OpenMM/   │             │              │ │
│  │   PDBBind)  │   PySCF)   │             │              │ │
│  └─────────────┴─────────────┴─────────────┴─────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

## Installation

### From Snap (recommended for Ubuntu/Debian)

```bash
# Install from Snap Store (when available)
sudo snap install ligora

# Or build from source
cd snap
snapcraft
sudo snap install ligora_*.snap --dangerous
```

### From Source

```bash
# Clone the repository
git clone https://github.com/ligora/ligora.git
cd ligora

# Install Python backend dependencies
cd backend
pip install numpy requests

# Install Tauri frontend dependencies
cd frontend
npm install

# Build and run
cd frontend
npm run tauri dev
```

## Usage

### Opening a Structure

1. Click **Open File** to load a local PDB/mmCIF file
2. Click **Open PDB ID** and enter a PDB ID (e.g., `1ABC`)

### Running Contact Analysis

1. Select a ligand from the Ligand Card panel
2. Click **Run Analysis** to detect protein-ligand contacts
3. View contacts in the Contacts Table panel
4. Check Evidence pane for external data enrichment

### Running Docking (V2)

1. Select a ligand
2. Configure the docking box (center and size)
3. Click **Run Docking**
4. View poses in the scene and results in the Jobs panel

## Engine Integration

Ligora integrates with existing open-source engines rather than implementing them:

| Engine | Purpose | License |
|--------|---------|---------|
| AutoDock Vina | Molecular docking | Apache 2.0 |
| gnina | Deep-learning docking | GPL/Apache |
| PLIP | Contact analysis | GPL |
| OpenMM | MD simulations (V2) | LGPL |
| PySCF | QM calculations (V2) | Apache 2.0 |
| 3Dmol.js | WebGL viewer | BSD-3 |

## Data Sources

Ligora enriches analysis with live data from:

- **RCSB PDB**: Structure metadata, annotations, related structures
- **PubChem**: Compound identity, properties, bioactivity
- **ChEMBL**: Bioactivity data, targets, binding data
- **PDBBind**: Binding affinity data for protein-ligand complexes

## Packaging

### Snap (x86_64)

The snap package includes:
- Tauri frontend (Rust + Web)
- Python backend with all dependencies
- Pre-installed engines (Vina, PLIP)
- Access to user's home directory and network

```yaml
# snap/snapcraft.yaml
name: ligora
version: '0.1.0'
base: core22
confinement: strict
```

### Building

```bash
# Build snap package
cd snap
snapcraft

# Install locally
sudo snap install ligora_*.snap --dangerous
```

## Development

### Project Structure

```
ligora/
├── backend/
│   ├── ligora_backend/
│   │   ├── __init__.py
│   │   ├── config.py           # Configuration
│   │   ├── workspace.py        # Session/workspace management
│   │   ├── parser.py           # Structure parsing (mmCIF/PDB)
│   │   ├── ligand.py           # Ligand resolution
│   │   ├── contacts.py         # Contact analysis (PLIP + geometry)
│   │   ├── cheminformatics.py  # Fingerprints and similarity
│   │   ├── enrichment.py       # Live data clients
│   │   ├── engines.py          # Engine adapters (Vina, gnina, etc.)
│   │   ├── jobs.py             # Job management
│   │   ├── export.py           # Artifact export
│   │   └── server.py           # IPC server
│   ├── tests/
│   └── requirements.txt
├── frontend/
│   ├── src/
│   │   ├── main.js            # App entry point
│   │   ├── app.js             # Main app component
│   │   └── styles.css         # Application styles
│   ├── src-tauri/
│   │   ├── src/
│   │   │   ├── lib.rs         # Tauri backend
│   │   │   └── main.rs        # Rust entry point
│   │   └── Cargo.toml
│   ├── package.json
│   └── tauri.conf.json
├── snap/
│   └── snapcraft.yaml
├── .github/
│   └── workflows/
│       └── ci.yml
└── README.md
```

### Running Tests

```bash
cd backend
pip install -r requirements.txt
pip install pytest
python -m pytest tests/ -v
```

### Code Style

```bash
cd backend
pip install ruff
ruff check . --select E,F,W
ruff format .
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

### Third-party Licenses

Components used by Ligora have their own licenses:

- **AutoDock Vina**: Apache 2.0
- **PLIP**: GPL v2
- **gnina**: GPL v2 / Apache 2.0 (dual)
- **OpenMM**: LGPL
- **PySCF**: Apache 2.0
- **3Dmol.js**: BSD-3-Clause

See the THIRD_PARTY_LICENSES file for details.

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

## Acknowledgments

- [RCSB PDB](https://www.rcsb.org/) for structure data
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) for compound data
- [ChEMBL](https://www.ebi.ac.uk/chembl/) for bioactivity data
- [PDBBind](http://www.pdbbind.org.cn/) for binding affinity data
- [AutoDock Vina](https://github.com/ccsb-scripps/autoDock-Vina) for docking
- [PLIP](https://github.com/pharmai/plip) for contact analysis

## Contact

- GitHub: https://github.com/ligora/ligora
- Issues: https://github.com/ligora/ligora/issues
