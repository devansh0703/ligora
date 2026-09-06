"""
Configuration management for Ligora backend.
"""
import os
from pathlib import Path
from typing import Optional, Dict, List


class Config:
    """Application configuration."""

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

        # Network
        self.rcsb_base_url: str = os.environ.get(
            "LIGORA_RSCB_URL",
            "https://data.rcsb.org"
        )
        self.rcsb_search_url: str = os.environ.get(
            "LIGORA_RSCB_SEARCH_URL",
            "https://search.rcsb.org"
        )
        self.pubchem_base_url: str = os.environ.get(
            "LIGORA_PUBCHEM_URL",
            "https://pubchem.ncbi.nlm.nih.gov"
        )
        self.chembl_base_url: str = os.environ.get(
            "LIGORA_CHEMBL_URL",
            "https://www.ebi.ac.uk/chembl/api/data"
        )
        self.pdbbind_url: str = os.environ.get(
            "LIGORA_PDBBIND_URL",
            "https://www.pdbbind-plus.org.cn"
        )

        # API settings
        self.request_timeout: int = int(os.environ.get(
            "LIGORA_REQUEST_TIMEOUT", "30"
        ))
        self.max_retries: int = int(os.environ.get(
            "LIGORA_MAX_RETRIES", "3"
        ))
        self.cache_ttl_seconds: int = int(os.environ.get(
            "LIGORA_CACHE_TTL", "86400"  # 24 hours
        ))

        # Engine paths (can be set via env or discovered)
        self.vina_executable: Optional[str] = os.environ.get("VINA_PATH")
        self.gnina_executable: Optional[str] = os.environ.get("GNINA_PATH")
        self.openmm_python_path: Optional[str] = os.environ.get("OPENMM_PYTHON")

        # Job settings
        self.default_docking_exhaustiveness: int = int(os.environ.get(
            "LIGORA_DOCKING_EXHAUSTIVENESS", "8"
        ))
        self.default_docking_num_modes: int = int(os.environ.get(
            "LIGORA_DOCKING_NUM_MODES", "9"
        ))
        self.default_docking_energy_range: float = float(os.environ.get(
            "LIGORA_DOCKING_ENERGY_RANGE", "3.0"
        ))

        # Geometry cleanup settings
        self.geometry_cleanup_max_iterations: int = int(os.environ.get(
            "LIGORA_GEOM_CLEANUP_ITERATIONS", "50"
        ))
        self.geometry_cleanup_tolerance: float = float(os.environ.get(
            "LIGORA_GEOM_CLEANUP_TOLERANCE", "1e-6"
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
        """Get the path to an engine executable."""
        if engine_name == "vina":
            return self.vina_executable
        elif engine_name == "gnina":
            return self.gnina_executable
        return None

    def set_engine_path(self, engine_name: str, path: str):
        """Set the path to an engine executable."""
        if engine_name == "vina":
            self.vina_executable = path
        elif engine_name == "gnina":
            self.gnina_executable = path

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return {
            "base_dir": str(self.base_dir),
            "workspace_dir": str(self.workspace_dir),
            "cache_dir": str(self.cache_dir),
            "rcsb_base_url": self.rcsb_base_url,
            "pubchem_base_url": self.pubchem_base_url,
            "chembl_base_url": self.chembl_base_url,
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
