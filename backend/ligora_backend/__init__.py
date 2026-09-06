"""
Ligora Backend - Analysis and Engine Orchestration

This backend provides:
- Structure parsing (mmCIF/PDB)
- Ligand detection and identity resolution
- Contact analysis (PLIP + geometry)
- Cheminformatics (fingerprints, similarity, SDF/SMILES)
- Live enrichment (RCSB, PubChem, ChEMBL, PDBBind)
- Engine orchestration (Vina, gnina, OpenMM, PySCF)
- Job/session management
- Artifact export
"""

from .workspace import WorkspaceManager, Session
from .parser import StructureParser
from .ligand import LigandResolver
from .contacts import ContactAnalyzer
from .cheminformatics import Cheminformatics
from .enrichment import EnrichmentClient
from .engines import EngineAdapter, VinaAdapter, GninaAdapter, EngineRegistry
from .jobs import JobManager
from .export import ArtifactExporter
from .server import BackendServer, start_server
from .config import Config, get_config, reset_config

__version__ = "0.1.0"

__all__ = [
    "WorkspaceManager",
    "Session",
    "StructureParser",
    "LigandResolver",
    "ContactAnalyzer",
    "Cheminformatics",
    "EnrichmentClient",
    "EngineAdapter",
    "VinaAdapter",
    "GninaAdapter",
    "EngineRegistry",
    "JobManager",
    "ArtifactExporter",
    "BackendServer",
    "start_server",
    "Config",
    "get_config",
    "reset_config",
]
