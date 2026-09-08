"""
Configuration management for Ligora backend.

Every external service URL, timeout, and analysis parameter is a user setting
(environment variable). No chemical or biological knowledge lives here.
"""
import os
from pathlib import Path
from typing import Optional, Dict, Any


class Config:
    """Application configuration (all values from environment or defaults)."""

    def __init__(self):
        # Paths
        self.base_dir: Path = Path(os.environ.get(
            "LIGORA_BASE_DIR",
            str(Path.home() / ".ligora")
        ))
        self.workspace_dir: Path = self.base_dir / "workspaces"
        self.cache_dir: Path = self.base_dir / "cache"
        self.logs_dir: Path = self.base_dir / "logs"
        self.engines_dir: Path = self.base_dir / "engines"

        # Network endpoints (RCSB Data API lives under /rest/v1)
        self.rcsb_data_url: str = os.environ.get(
            "LIGORA_RCSB_DATA_URL",
            "https://data.rcsb.org/rest/v1",
        )
        self.rcsb_files_url: str = os.environ.get(
            "LIGORA_RCSB_FILES_URL",
            "https://files.rcsb.org",
        )
        self.pubchem_base_url: str = os.environ.get(
            "LIGORA_PUBCHEM_URL",
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        )
        self.chembl_base_url: str = os.environ.get(
            "LIGORA_CHEMBL_URL",
            "https://www.ebi.ac.uk/chembl/api/data",
        )
        self.unichem_base_url: str = os.environ.get(
            "LIGORA_UNICHEM_URL",
            "https://www.ebi.ac.uk/unichem/rest",
        )
        # PDBBind has no open public API. When the user configures an
        # endpoint it is queried; otherwise affinity lookups report
        # "not available" instead of synthesizing values.
        self.pdbbind_url: Optional[str] = os.environ.get("LIGORA_PDBBIND_URL")
        # BindingDB REST web services (real affinity records, Ki/Kd/IC50).
        self.bindingdb_base_url: str = os.environ.get(
            "LIGORA_BINDINGDB_URL",
            "https://bindingdb.org",
        )

        # API settings
        self.request_timeout: int = int(os.environ.get(
            "LIGORA_REQUEST_TIMEOUT", "30"
        ))
        self.max_retries: int = int(os.environ.get(
            "LIGORA_MAX_RETRIES", "2"
        ))
        self.cache_ttl_seconds: int = int(os.environ.get(
            "LIGORA_CACHE_TTL", "86400"  # 24 hours
        ))

        # External tool executables (resolved from PATH when not set)
        self.vina_executable: Optional[str] = os.environ.get("VINA_PATH")
        self.gnina_executable: Optional[str] = os.environ.get("GNINA_PATH")
        self.obabel_executable: Optional[str] = os.environ.get("OBABEL_PATH")
        self.plip_executable: Optional[str] = os.environ.get("PLIP_PATH")

        # Tool timeouts
        self.plip_timeout: int = int(os.environ.get("LIGORA_PLIP_TIMEOUT", "300"))
        self.vina_timeout: int = int(os.environ.get("LIGORA_VINA_TIMEOUT", "3600"))

        # MD engine (GROMACS; user-provided, never simulated)
        self.gromacs_executable: Optional[str] = os.environ.get(
            "LIGORA_GMX_PATH")
        self.md_timeout: int = int(os.environ.get(
            "LIGORA_MD_TIMEOUT", "7200"))
        # Default force field for pdb2gmx when the user does not choose one.
        # None means GROMACS selects from its own shipped force-field list
        # (its engine-owned default); the chosen field is recorded in the
        # result. A configuration value, not chemistry knowledge.
        self.md_force_field: Optional[str] = os.environ.get(
            "LIGORA_MD_FORCEFIELD")
        self.md_water_model: str = os.environ.get(
            "LIGORA_MD_WATER", "none")

        # Docking settings
        self.default_docking_exhaustiveness: int = int(os.environ.get(
            "LIGORA_DOCKING_EXHAUSTIVENESS", "8"
        ))
        self.default_docking_num_modes: int = int(os.environ.get(
            "LIGORA_DOCKING_NUM_MODES", "9"
        ))
        self.default_docking_energy_range: float = float(os.environ.get(
            "LIGORA_DOCKING_ENERGY_RANGE", "3.0"
        ))

        # Geometry cleanup (Open Babel force field / minimization settings)
        self.geometry_cleanup_force_field: str = os.environ.get(
            "LIGORA_GEOM_FF", "mmff94"
        )
        self.geometry_cleanup_steps: int = int(os.environ.get(
            "LIGORA_GEOM_STEPS", "200"
        ))
        self.geometry_cleanup_tolerance: float = float(os.environ.get(
            "LIGORA_GEOM_TOLERANCE", "1e-6"
        ))

        # Water network analysis parameters (user settings, not chemistry rules)
        self.water_network_cutoff: float = float(os.environ.get(
            "LIGORA_WATER_NETWORK_CUTOFF", "3.5"
        ))
        self.water_contact_cutoff: float = float(os.environ.get(
            "LIGORA_WATER_CONTACT_CUTOFF", "4.0"
        ))

        # Fingerprint settings
        self.fingerprint_radius: int = int(os.environ.get(
            "LIGORA_FP_RADIUS", "2"
        ))
        self.fingerprint_nbits: int = int(os.environ.get(
            "LIGORA_FP_NBITS", "2048"
        ))

        # Ensure directories exist
        self._ensure_dirs()

    def _ensure_dirs(self):
        """Create necessary directories."""
        for d in [self.workspace_dir, self.cache_dir,
                  self.logs_dir, self.engines_dir]:
            d.mkdir(parents=True, exist_ok=True)

    def get_engine_path(self, engine_name: str) -> Optional[str]:
        """Get the configured path for an engine executable, if any."""
        return {
            "vina": self.vina_executable,
            "gnina": self.gnina_executable,
            "obabel": self.obabel_executable,
            "plip": self.plip_executable,
        }.get(engine_name)

    def set_engine_path(self, engine_name: str, path: str):
        """Set the path to an engine executable."""
        if engine_name == "vina":
            self.vina_executable = path
        elif engine_name == "gnina":
            self.gnina_executable = path
        elif engine_name == "obabel":
            self.obabel_executable = path
        elif engine_name == "plip":
            self.plip_executable = path

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return {
            "base_dir": str(self.base_dir),
            "workspace_dir": str(self.workspace_dir),
            "cache_dir": str(self.cache_dir),
            "rcsb_data_url": self.rcsb_data_url,
            "rcsb_files_url": self.rcsb_files_url,
            "pubchem_base_url": self.pubchem_base_url,
            "chembl_base_url": self.chembl_base_url,
            "unichem_base_url": self.unichem_base_url,
        }


# Global config instance
_config: Optional[Config] = None


def get_config() -> Config:
    """Get the global configuration instance."""
    global _config
    if _config is None:
        _config = Config()
    return _config


def reset_config():
    """Reset the global configuration (useful for testing)."""
    global _config
    _config = None
