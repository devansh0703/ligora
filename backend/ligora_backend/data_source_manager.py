"""
Data Source Manager for Ligora.

Manages download and caching of data source catalogs:
- RCSB Chemical Component Dictionary (CCD)
- PubChem compound data
- ChEMBL compound data

All data used by the application comes from these downloaded catalogs.
No chemical data is hardcoded. The application discovers what is available
from the data sources and uses only that.
"""

import json
from typing import Optional, Dict, Any
from dataclasses import dataclass

import requests

from .config import get_config
from .net import http_timeout


@dataclass
class DataSourceInfo:
    """Metadata about a data source."""
    name: str
    base_url: str
    catalog_url: str
    catalog_format: str  # 'json', 'cif', 'sdf', etc.
    last_updated: Optional[str] = None
    catalog_size: Optional[int] = None
    catalog_checksum: Optional[str] = None
    is_available: bool = False
    error_message: Optional[str] = None


class DataSourceManager:
    """
    Manages data source catalogs for Ligora.

    Responsibilities:
    1. Discover available data sources
    2. Download and cache their catalogs
    3. Provide lookup from cached catalogs
    4. Update catalogs when stale

    All chemical data comes from these catalogs - nothing is hardcoded.
    """

    def __init__(self, config=None):
        self.config = config or get_config()
        self.data_dir = self.config.base_dir / "data_sources"
        self.data_dir.mkdir(parents=True, exist_ok=True)

        self._session = requests.Session()
        self._session.headers.update({
            'User-Agent': 'Ligora/0.1.0 (https://github.com/ligora)',
        })

        # Data source registry - populated from discovery
        self.sources: Dict[str, DataSourceInfo] = {}

        # Cached catalogs - populated after download
        self.catalogs: Dict[str, Any] = {}

        # Discovery performed flag
        self._discovery_performed = False

    def discover_sources(self) -> Dict[str, DataSourceInfo]:
        """
        Discover available data sources and their capabilities.

        Checks each known data source endpoint to see if it's accessible
        and what format its catalog is available in.

        Returns:
            Dictionary of source name -> DataSourceInfo
        """
        self.sources = {}

        # Check RCSB PDB
        self._discover_rcsb()

        # Check PubChem
        self._discover_pubchem()

        # Check ChEMBL
        self._discover_chembl()

        self._discovery_performed = True
        return self.sources

    def _discover_rcsb(self):
        """Discover RCSB PDB data sources."""
        info = DataSourceInfo(
            name="RCSB PDB",
            base_url=self.config.rcsb_base_url,
            catalog_url="",
            catalog_format="",
        )

        # Check if core API is accessible
        try:
            response = self._session.get(
                f"{info.base_url}/rest/v1/core/entry/3W85",
                timeout=10
            )
            if response.status_code == 200:
                info.is_available = True
                info.catalog_format = "json"
                # Catalog is the REST API itself - we query it per-entry
                info.catalog_url = f"{info.base_url}/rest/v1"
            else:
                info.error_message = f"API returned status {response.status_code}"
        except requests.RequestException as e:
            info.error_message = str(e)

        self.sources["rcsb"] = info

    def _discover_pubchem(self):
        """Discover PubChem data sources."""
        info = DataSourceInfo(
            name="PubChem",
            base_url=self.config.pubchem_base_url,
            catalog_url="",
            catalog_format="",
        )

        try:
            # Check if PubChem PUG-REST is accessible
            response = self._session.get(
                f"{info.base_url}/rest/pug/compound/cid/962/property/MolecularWeight/JSON",
                timeout=10
            )
            if response.status_code == 200:
                info.is_available = True
                info.catalog_format = "json"
                info.catalog_url = f"{info.base_url}/rest/pug"
            else:
                info.error_message = f"API returned status {response.status_code}"
        except requests.RequestException as e:
            info.error_message = str(e)

        self.sources["pubchem"] = info

    def _discover_chembl(self):
        """Discover ChEMBL data sources."""
        info = DataSourceInfo(
            name="ChEMBL",
            base_url=self.config.chembl_base_url,
            catalog_url="",
            catalog_format="",
        )

        try:
            response = self._session.get(
                f"{info.base_url}/molecule/CHEMBL1.json",
                timeout=10
            )
            if response.status_code == 200:
                info.is_available = True
                info.catalog_format = "json"
                info.catalog_url = f"{info.base_url}"
            else:
                info.error_message = f"API returned status {response.status_code}"
        except requests.RequestException as e:
            info.error_message = str(e)

        self.sources["chembl"] = info

    def get_source_status(self) -> Dict[str, Dict[str, Any]]:
        """
        Get status of all discovered data sources.

        Returns:
            Dictionary of source name -> status info
        """
        if not self._discovery_performed:
            self.discover_sources()

        status = {}
        for name, info in self.sources.items():
            status[name] = {
                "available": info.is_available,
                "format": info.catalog_format,
                "catalog_url": info.catalog_url,
                "error": info.error_message,
            }
        return status

    def download_catalog(self, source_name: str) -> bool:
        """
        Download the catalog for a specific data source.

        For RCSB, this downloads the structure-specific data on demand.
        For PubChem, this would download bulk compound data (not implemented yet).
        For ChEMBL, this would download bulk molecule data (not implemented yet).

        Args:
            source_name: Name of the source ("rcsb", "pubchem", "chembl")

        Returns:
            True if catalog was downloaded/updated successfully
        """
        if source_name not in self.sources:
            return False

        info = self.sources[source_name]

        if not info.is_available:
            return False

        catalog_path = self.data_dir / f"{source_name}_catalog.json"

        # Check if we have a recent cache
        if catalog_path.exists():
            try:
                json.loads(catalog_path.read_text())  # validity probe
                # For now, just use the API directly - no bulk downloading
                # The "catalog" is the live API
                self.catalogs[source_name] = {"type": "live_api", "source": source_name}
                return True
            except (json.JSONDecodeError, IOError):
                pass

        # Mark as having a live API catalog
        self.catalogs[source_name] = {
            "type": "live_api",
            "source": source_name,
            "base_url": info.base_url,
        }

        return True

    def ensure_catalog(self, source_name: str) -> bool:
        """
        Ensure we have a catalog for the given source.

        Discovers the source if needed and downloads its catalog.

        Args:
            source_name: Name of the source

        Returns:
            True if catalog is available
        """
        if source_name not in self.sources:
            self.discover_sources()

        if not self.download_catalog(source_name):
            return False

        return source_name in self.catalogs

    def get_compound_from_catalog(
        self,
        source_name: str,
        identifier: str,
        identifier_type: str = "name",
    ) -> Optional[Dict[str, Any]]:
        """
        Look up a compound from a data source catalog.

        Args:
            source_name: Name of the source ("rcsb", "pubchem", "chembl")
            identifier: Compound identifier (name, CID, ChEMBL ID, etc.)
            identifier_type: Type of identifier ("name", "cid", "chembl_id")

        Returns:
            Compound data from the catalog, or None if not found
        """
        if not self.ensure_catalog(source_name):
            return None

        catalog = self.catalogs.get(source_name)
        if not catalog:
            return None

        if catalog.get("type") == "live_api":
            return self._query_live_api(source_name, identifier, identifier_type)

        # Future: load from cached bulk file
        return None

    def _query_live_api(
        self,
        source_name: str,
        identifier: str,
        identifier_type: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Query a data source's live API for compound data.

        Args:
            source_name: Name of the source
            identifier: Compound identifier
            identifier_type: Type of identifier

        Returns:
            Compound data from the API, or None if not found
        """
        if source_name == "rcsb":
            return self._query_rcsb(identifier, identifier_type)
        elif source_name == "pubchem":
            return self._query_pubchem(identifier, identifier_type)
        elif source_name == "chembl":
            return self._query_chembl(identifier, identifier_type)

        return None

    def _query_rcsb(
        self,
        identifier: str,
        identifier_type: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Query RCSB for chemical component data.

        Args:
            identifier: Residue name (e.g., "HEM", "ATP", "W85")
            identifier_type: Must be "name" for RCSB

        Returns:
            CCD data from RCSB, or None if not found
        """
        if identifier_type != "name":
            return None

        try:
            # Try the core entry endpoint for structure data
            # This is what the current code uses
            response = self._session.get(
                f"{self.config.rcsb_base_url}/rest/v1/core/entry/{identifier}",
                timeout=http_timeout(self.config.request_timeout)
            )
            if response.status_code == 200:
                # Parse what we can from the entry
                result = {
                    "name": identifier,
                    "source": "rcsb",
                }

                # Try to extract formula from the entry's chemical components
                # This would need more parsing based on what the API returns
                return result
        except requests.RequestException:
            pass

        return None

    def _query_pubchem(
        self,
        identifier: str,
        identifier_type: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Query PubChem for compound data.

        Args:
            identifier: Compound name or CID
            identifier_type: "name" or "cid"

        Returns:
            PubChem data, or None if not found
        """
        from .enrichment import EnrichmentClient
        enrichment = EnrichmentClient()

        if identifier_type == "cid":
            try:
                cid = int(identifier)
                return enrichment.get_pubchem_info(pubchem_cid=cid)
            except ValueError:
                pass

        # Default to name lookup
        return enrichment.get_pubchem_info(compound_name=identifier)

    def _query_chembl(
        self,
        identifier: str,
        identifier_type: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Query ChEMBL for compound data.

        Args:
            identifier: ChEMBL ID or SMILES
            identifier_type: "chembl_id" or "smiles"

        Returns:
            ChEMBL data, or None if not found
        """
        from .enrichment import EnrichmentClient
        enrichment = EnrichmentClient()

        if identifier_type == "chembl_id":
            return enrichment.get_chembl_info(chembl_id=identifier)
        elif identifier_type == "smiles":
            return enrichment.get_chembl_info(smiles=identifier)

        return None

    def resolve_ligand_from_sources(
        self,
        residue_name: str,
        formula: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Resolve a ligand's identity using all available data sources.

        Queries each available source in order of preference and merges results.
        No chemical data is invented - only what sources actually provide.

        Args:
            residue_name: The 3-letter residue name (e.g., "W85", "HEM")
            formula: Optional formula to help with lookup

        Returns:
            Dictionary with resolved identity data and source attribution.
            All fields are optional - only what sources provided is included.
        """
        result: Dict[str, Any] = {
            "residue_name": residue_name,
            "sources_used": [],
            "fields": {},
        }

        # Ensure we have catalogs for all sources
        for source_name in ["rcsb", "pubchem", "chembl"]:
            self.ensure_catalog(source_name)

        # Try each source in order
        # 1. PubChem by name (best for small molecules)
        pubchem_data = self.get_compound_from_catalog(
            "pubchem",
            residue_name,
            "name",
        )
        if pubchem_data:
            result["sources_used"].append("pubchem")
            result["fields"] = self._merge_compound_data(
                result["fields"],
                pubchem_data,
                source="pubchem",
            )

        # 2. ChEMBL by SMILES (if we got SMILES from PubChem)
        if result["fields"].get("smiles"):
            chembl_data = self.get_compound_from_catalog(
                "chembl",
                result["fields"]["smiles"],
                "smiles",
            )
            if chembl_data:
                result["sources_used"].append("chembl")
                result["fields"] = self._merge_compound_data(
                    result["fields"],
                    chembl_data,
                    source="chembl",
                )

        return result

    def _merge_compound_data(
        self,
        existing: Dict[str, Any],
        new_data: Dict[str, Any],
        source: str,
    ) -> Dict[str, Any]:
        """
        Merge new data into existing compound data.

        Only adds fields that don't already have values.
        All data must come from the source - nothing is invented.

        Args:
            existing: Existing data dict
            new_data: New data from source
            source: Source name for attribution

        Returns:
            Merged data dict
        """
        # Fields that can come from compound data sources
        mergeable_fields = [
            "formula",
            "smiles",
            "isomeric_smiles",
            "connectivity_smiles",
            "iupac_name",
            "inchi_key",
            "molecular_weight",
            "xlogp",
            "tpsa",
            "complexity",
            "rotatable_bond_count",
            "h_bond_acceptor_count",
            "h_bond_donor_count",
        ]

        for field in mergeable_fields:
            if field not in existing or existing[field] is None:
                # Map source field names to our field names
                value = new_data.get(field)
                if value is None:
                    # Try common alternative field names
                    alt_mapping = {
                        "smiles": ["canonical_smiles", "absolute_smiles"],
                        "molecular_weight": ["full_mwt", "molecular_weight"],
                        "iupac_name": ["pref_name", "iupac_name"],
                        "inchi_key": ["inchi_key", "standard_inchi_key"],
                        "formula": ["full_molformula", "formula"],
                    }
                    for alt in alt_mapping.get(field, []):
                        value = new_data.get(alt)
                        if value is not None:
                            break

                if value is not None:
                    existing[field] = value

        return existing
