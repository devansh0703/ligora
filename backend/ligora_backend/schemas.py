"""
Shared schemas and types for Ligora IPC.
These are used by both the Tauri frontend (via JSON) and the Python backend.
"""
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any, Literal
from enum import Enum
import json


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class EngineType(str, Enum):
    VINA = "vina"
    GNINA = "gnina"
    OPENMM = "openmm"
    PYSCF = "pyscf"
    GEOMETRY_CLEANUP = "geometry_cleanup"


class ContactType(str, Enum):
    HYDROGEN_BOND = "hydrogen_bond"
    HYDROPHOBIC = "hydrophobic"
    PI_STACKING = "pi_stacking"
    SALT_BRIDGE = "salt_bridge"
    HALOGEN_BOND = "halogen_bond"
    METAL_COORDINATION = "metal_coordination"
    WATER_MEDIATED = "water_mediated"
    VAN_DER_WAALS = "van_der_waals"
    UNKNOWN = "unknown"


class LigandResolutionStatus(str, Enum):
    RESOLVED = "resolved"
    PARTIAL = "partial"
    NOT_FOUND = "not_found"


@dataclass
class Atom:
    """Represents a single atom in a molecular structure."""
    id: int
    name: str
    residue_name: str
    residue_id: int
    chain_id: str
    x: float
    y: float
    z: float
    element: Optional[str] = None
    b_factor: float = 0.0
    occupancy: float = 1.0


@dataclass
class Residue:
    """Represents a residue in a molecular structure."""
    id: int
    name: str
    chain_id: str
    atoms: List[Atom] = field(default_factory=list)
    residue_number: int = 0


@dataclass
class Chain:
    """Represents a polymer chain."""
    id: str
    name: str
    residues: List[Residue] = field(default_factory=list)
    is_polymer: bool = True


@dataclass
class Ligand:
    """Represents a non-polymer ligand entity."""
    id: str
    name: str
    residue_name: str
    formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    atom_count: int = 0
    smiles: Optional[str] = None
    inchi_key: Optional[str] = None
    atoms: List[Atom] = field(default_factory=list)
    resolution_status: LigandResolutionStatus = LigandResolutionStatus.NOT_FOUND
    pubchem_cid: Optional[int] = None
    pubchem_name: Optional[str] = None
    chembl_id: Optional[str] = None
    pdbbind_affinity: Optional[float] = None
    classification_hint: str = "unknown"
    has_2d_structure: bool = False
    iupac_name: Optional[str] = None


@dataclass
class Structure:
    """Represents a complete molecular structure (PDB/mmCIF)."""
    id: str
    title: str = "Unknown structure"
    source: str = "local"  # "local", "rcsb", "pdb", "csm"
    chains: List[Chain] = field(default_factory=list)
    ligands: List[Ligand] = field(default_factory=list)
    has_density: bool = False
    resolution: Optional[float] = None
    experiment_type: Optional[str] = None
    file_path: Optional[str] = None
    file_format: str = "mmcif"


@dataclass
class Contact:
    """Represents a protein-ligand interaction contact."""
    id: int
    ligand_atom: str
    ligand_residue_name: str
    ligand_residue_id: int
    ligand_chain_id: str
    protein_residue_name: str
    protein_residue_id: int
    protein_chain_id: str
    protein_atom: str
    distance: float
    contact_type: ContactType
    angle: Optional[float] = None
    is_water_mediated: bool = False
    description: str = ""


@dataclass
class DockingPose:
    """Represents a docking pose result."""
    pose_id: int
    affinity: float  # kcal/mol
    rmsd_bound: Optional[float] = None
    rmsd_unbound: Optional[float] = None
    atoms: List[Atom] = field(default_factory=list)
    ligand_name: str = ""


@dataclass
class DockingResult:
    """Result from a docking run."""
    job_id: str
    engine: EngineType
    poses: List[DockingPose] = field(default_factory=list)
    box_center: Optional[List[float]] = None
    box_size: Optional[List[float]] = None
    runtime_seconds: float = 0.0
    output_file: Optional[str] = None


@dataclass
class GeometryCleanupResult:
    """Result from local geometry cleanup."""
    job_id: str
    original_energy: Optional[float] = None
    cleaned_energy: Optional[float] = None
    rmsd: Optional[float] = None
    strain_flag: bool = False
    atoms: List[Atom] = field(default_factory=list)


@dataclass
class EvidenceItem:
    """A piece of evidence from external sources."""
    source: str  # "pubchem", "chembl", "pdbbind", "rcsb"
    field: str
    value: Any
    score: Optional[float] = None
    url: Optional[str] = None
    status: str = "available"


@dataclass
class AnalysisSummary:
    """Complete analysis summary for export."""
    structure_id: str
    structure_title: str
    source: str
    ligand_id: str
    ligand_name: str
    ligand_formula: Optional[str]
    ligand_smiles: Optional[str]
    resolution_status: str
    contact_count: int
    contacts: List[Dict[str, Any]]
    evidence: List[Dict[str, Any]]
    job_results: List[Dict[str, Any]]
    notes: str
    created_at: str
    exported_at: str
    scene_image_path: Optional[str] = None
    ligand_export_path: Optional[str] = None
    contact_export_path: Optional[str] = None


# IPC command types - what frontend can send to backend
class CommandType(str, Enum):
    OPEN_LOCAL_FILE = "open_local_file"
    OPEN_PDB_ID = "open_pdb_id"
    SELECT_LIGAND = "select_ligand"
    RUN_CONTACT_ANALYSIS = "run_contact_analysis"
    RUN_DOCKING = "run_docking"
    RUN_GEOMETRY_CLEANUP = "run_geometry_cleanup"
    SET_DOCKING_BOX = "set_docking_box"
    CANCEL_JOB = "cancel_job"
    MEASURE_DISTANCE = "measure_distance"
    MEASURE_ANGLE = "measure_angle"
    SAVE_ARTIFACT = "save_artifact"
    LOAD_ARTIFACT = "load_artifact"
    GET_STATUS = "get_status"
    EXPORT_SCENE_IMAGE = "export_scene_image"


@dataclass
class Command:
    """A command from frontend to backend."""
    type: CommandType
    payload: Dict[str, Any]
    session_id: str
    command_id: str


@dataclass
class CommandResponse:
    """Response from backend to frontend."""
    command_id: str
    success: bool
    data: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


@dataclass
class Job:
    """A running or completed job."""
    id: str
    type: EngineType
    status: JobStatus
    parameters: Dict[str, Any]
    result: Optional[Any] = None
    progress: float = 0.0
    log: str = ""
    started_at: Optional[str] = None
    finished_at: Optional[str] = None
    error: Optional[str] = None


@dataclass
class Session:
    """Analysis session state."""
    id: str
    structure: Optional[Structure] = None
    selected_ligand_id: Optional[str] = None
    jobs: Dict[str, Job] = field(default_factory=dict)
    active_job_id: Optional[str] = None
    notes: str = ""
    created_at: str = ""
    workspace_path: str = ""


# Serialization helpers
def command_to_json(cmd: Command) -> str:
    """Serialize a command to JSON for IPC."""
    d = asdict(cmd)
    d["type"] = cmd.type.value
    return json.dumps(d)


def command_from_json(data: str) -> Command:
    """Deserialize a command from JSON."""
    d = json.loads(data)
    return Command(
        type=CommandType(d["type"]),
        payload=d.get("payload", {}),
        session_id=d.get("session_id", ""),
        command_id=d.get("command_id", "")
    )


def response_to_json(resp: CommandResponse) -> str:
    """Serialize a response to JSON."""
    d = asdict(resp)
    return json.dumps(d)


def response_from_json(data: str) -> CommandResponse:
    """Deserialize a response from JSON."""
    d = json.loads(data)
    return CommandResponse(
        command_id=d.get("command_id", ""),
        success=d.get("success", False),
        data=d.get("data"),
        error=d.get("error")
    )
