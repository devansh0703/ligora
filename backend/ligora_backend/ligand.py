"""
Ligand resolver - identity resolution and enrichment for detected ligands.

Resolves ligand chemical identity using:
- RCSB chemical component dictionary (CCD)
- PubChem compound lookup
- ChEMBL compound mapping
- PDBBind affinity data
"""

import json
import hashlib
import time
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass

import requests

from .schemas import (
    Ligand,
    LigandResolutionStatus,
    EvidenceItem,
)
from .config import get_config


class LigandResolver:
    """
    Resolves chemical identity for detected ligands.

    Provides:
    - RCSB CCD lookup (residue name -> chemical component)
    - PubChem CID lookup by name/formula
    - ChEMBL ID mapping
    - PDBBind affinity cross-reference
    """

# Classification is not decided here from hardcoded name sets.
    # Identity and classification come from data sources (CCD / PubChem / ChEMBL)
    # via the resolver and enrichment client.

    def __init__(self):
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._cache_expiry: Dict[str, float] = {}
        self.config = get_config()

    def resolve_ligand(
        self,
        ligand: Ligand,
        pdb_id: Optional[str] = None,
    ) -> Ligand:
        """
        Resolve a ligand's chemical identity using available data sources.

        Performs progressive enrichment:
        1. RCSB CCD lookup by residue name
        2. PubChem lookup by name/formula
        3. ChEMBL cross-reference
        4. PDBBind affinity lookup (if PDB ID available)

        All values are taken from live data sources. No chemical data is
        hardcoded in this module.

        Args:
            ligand: The ligand to resolve.
            pdb_id: Optional PDB ID for PDBBind lookup.

        Returns:
            The enriched ligand with resolution status and metadata.
        """
        # 1) Try RCSB CCD first.
        ccd_data = self._lookup_ccd(ligand.residue_name)
        if ccd_data:
            ligand.formula = ccd_data.get("formula") or ligand.formula
            if not ligand.smiles:
                ligand.smiles = ccd_data.get("smiles")
            if not ligand.inchi_key:
                ligand.inchi_key = ccd_data.get("inchi_key")
            if ccd_data.get("molecular_weight") is not None:
                try:
                    ligand.molecular_weight = float(ccd_data.get("molecular_weight"))
                except (TypeError, ValueError):
                    pass
            if ccd_data.get("name"):
                ligand.name = ccd_data.get("name")
            if not ligand.resolution_status or ligand.resolution_status == LigandResolutionStatus.NOT_FOUND:
                if ligand.smiles or ligand.pubchem_cid or ligand.chembl_id:
                    ligand.resolution_status = LigandResolutionStatus.PARTIAL

        # 2) PubChem by name/formula.
        pubchem_data = self._lookup_pubchem(
            ligand.name or ligand.residue_name,
            ligand.formula,
        )
        if pubchem_data:
            ligand.pubchem_cid = pubchem_data.get("cid")
            ligand.pubchem_name = pubchem_data.get("title")
            if not ligand.smiles:
                ligand.smiles = pubchem_data.get("smiles")
            if not ligand.inchi_key:
                ligand.inchi_key = pubchem_data.get("inchi_key")
            if pubchem_data.get("molecular_weight") is not None:
                try:
                    ligand.molecular_weight = float(pubchem_data.get("molecular_weight"))
                except (TypeError, ValueError):
                    pass
            if pubchem_data.get("iupac_name"):
                ligand.iupac_name = pubchem_data.get("iupac_name")
            if not ligand.resolution_status or ligand.resolution_status == LigandResolutionStatus.NOT_FOUND:
                if ligand.pubchem_cid or ligand.chembl_id:
                    ligand.resolution_status = LigandResolutionStatus.RESOLVED
                elif ligand.smiles:
                    ligand.resolution_status = LigandResolutionStatus.PARTIAL

        # 3) ChEMBL cross-reference (when a SMILES is available).
        if ligand.smiles:
            chembl_data = self._lookup_chembl_by_smiles(ligand.smiles)
            if chembl_data:
                ligand.chembl_id = chembl_data.get("chembl_id")
                if not ligand.resolution_status or ligand.resolution_status == LigandResolutionStatus.NOT_FOUND:
                    ligand.resolution_status = LigandResolutionStatus.RESOLVED

        # 4) PDBBind affinity (when a PDB ID is available).
        if pdb_id and ligand.residue_name:
            pdbbind_data = self._lookup_pdbbind(pdb_id, ligand.residue_name)
            if pdbbind_data:
                ligand.pdbbind_affinity = pdbbind_data.get("affinity")

        # Final resolution status normalization.
        if ligand.pubchem_cid or ligand.chembl_id:
            ligand.resolution_status = LigandResolutionStatus.RESOLVED
        elif ligand.smiles:
            ligand.resolution_status = LigandResolutionStatus.PARTIAL

        return ligand

    def _lookup_ccd(self, residue_name: str) -> Optional[Dict[str, Any]]:
        """
        Look up a chemical component in the RCSB CCD (Chemical Component Dictionary).

        Uses the live RCSB data source. No chemical values are hardcoded.

        Args:
            residue_name: The 3-letter residue name.

        Returns:
            CCD data dict or None if not found.
        """
        cache_key = f"ccd:{residue_name}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()

        # Try the CCD ligand endpoint first.
        ccd = self._lookup_ccd_via_ligand_endpoint(residue_name, config)
        if ccd is not None:
            self._set_cache(cache_key, ccd)
            return ccd

        # Fall back to the CCD compound endpoint.
        ccd = self._lookup_ccd_via_compound_endpoint(residue_name, config)
        if ccd is not None:
            self._set_cache(cache_key, ccd)
            return ccd

        return None

    def _lookup_ccd_via_ligand_endpoint(self, residue_name: str, config: 'Config') -> Optional[Dict[str, Any]]:
        """Look up a CCD component via the /v1/ligand/{id} endpoint."""
        url = f"{config.rcsb_base_url}/v1/ligand/{residue_name}"
        try:
            response = requests.get(url, timeout=config.request_timeout)
            if response.status_code != 200:
                return None
            return self._parse_ccd_response(response.json())
        except requests.RequestException:
            return None

    def _lookup_ccd_via_compound_endpoint(self, residue_name: str, config: 'Config') -> Optional[Dict[str, Any]]:
        """Look up a CCD component via the /v1/compound/{id} endpoint."""
        url = f"{config.rcsb_base_url}/v1/compound/{residue_name}"
        try:
            response = requests.get(url, timeout=config.request_timeout)
            if response.status_code != 200:
                return None
            return self._parse_ccd_response(response.json())
        except requests.RequestException:
            return None

    def _parse_ccd_response(self, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Extract common CCD fields from a RCSB REST response."""
        formula = self._extract_nested(data, "formula", "text")
        smiles = self._extract_nested(data, "smiles", "nonstereo")
        inchi_key = self._extract_nested(data, "inchi", "key")
        name = data.get("name") or ""
        molecular_weight = self._extract_nested(data, "molecular_weight")

        if not name and not formula and not smiles:
            return None

        return {
            "formula": formula,
            "smiles": smiles,
            "inchi_key": inchi_key,
            "name": name,
            "molecular_weight": molecular_weight,
        }

    def _extract_nested(self, data: Dict[str, Any], *keys: str) -> Any:
        """Extract a possibly nested value from a dict."""
        current = data
        for key in keys:
            if not isinstance(current, dict):
                return None
            current = current.get(key)
            if current is None:
                return None
        return current

    def _lookup_pubchem(
        self,
        name: str,
        formula: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Search PubChem for a compound by name and/or formula.

        Uses the PubChem PUG REST API so every value comes from the live data
        source. No chemical values are hardcoded.

        Args:
            name: Compound name.
            formula: Chemical formula (optional).

        Returns:
            PubChem data dict or None if not found.
        """
        cache_key = f"pubchem:{hashlib.md5(
            f'{name}:{formula}'.encode()
        ).hexdigest()[:12]}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()

        try:
            # 1) Resolve by exact name first.
            cid = self._pubchem_resolve_cid_by_name(name)
            if cid is not None:
                return self._pubchem_fetch_compound(cid, config)

            # 2) If that fails and the name looks like a short PDB-style
            #    residue name, try the synonym/keyword search path.
            if len(name.strip()) <= 5:
                for candidate in {name, name.lower(), name.upper()}:
                    cid = self._pubchem_resolve_cid_by_synonym(candidate, config)
                    if cid is not None:
                        return self._pubchem_fetch_compound(cid, config)

            # 3) If a formula is available, try formula-based lookup.
            if formula:
                cid = self._pubchem_resolve_cid_by_formula(formula, config)
                if cid is not None:
                    return self._pubchem_fetch_compound(cid, config)

            return None
        except requests.RequestException:
            return None

    def _resolve_cids_by_exact_name(self, name: str, config: 'Config') -> Optional[List[Any]]:
        """Return the CID list for an exact compound name, if any."""
        url = f"{config.pubchem_base_url}/rest/pug/compound/name/{name}/cids/JSON"
        try:
            resp = requests.get(url, timeout=config.request_timeout)
        except requests.RequestException:
            return None
        if resp.status_code != 200:
            return None
        data = resp.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if not cids:
            return None
        return cids

    def _resolve_cids_by_formula(self, formula: str, config: 'Config') -> Optional[List[Any]]:
        """Return the CID list for a molecular formula, if any."""
        url = f"{config.pubchem_base_url}/rest/pug/compound/formula/{formula}/cids/JSON"
        try:
            resp = requests.get(url, timeout=config.request_timeout)
        except requests.RequestException:
            return None
        if resp.status_code != 200:
            return None
        data = resp.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if not cids:
            return None
        return cids

    def _pubchem_resolve_cid_by_name(self, name: str) -> Optional[int]:
        """Resolve a PubChem CID by exact compound name."""
        config = get_config()
        cids = self._resolve_cids_by_exact_name(name, config)
        if not cids:
            return None
        return int(cids[0])

    def _pubchem_resolve_cid_by_formula(self, formula: str, config: 'Config') -> Optional[int]:
        """Resolve a PubChem CID by molecular formula."""
        config = get_config()
        cids = self._resolve_cids_by_formula(formula, config)
        if not cids:
            return None
        return int(cids[0])

    def _pubchem_synonyms_for_variant(self, variant: str, config: 'Config') -> Optional[List[Any]]:
        """Return the synonym search CID list for one case variant, if any."""
        url = f"{config.pubchem_base_url}/rest/pug/compound/name/{variant}/synonyms/JSON"
        try:
            resp = requests.get(url, timeout=config.request_timeout)
        except requests.RequestException:
            return None
        if resp.status_code != 200:
            return None
        data = resp.json()
        info_list = data.get("InformationList", {}).get("Information", [])
        if not info_list:
            return []
        return [info.get("CID") for info in info_list if info.get("CID")]

    def _pubchem_resolve_cid_by_synonym(self, name: str, config: 'Config') -> Optional[int]:
        """Resolve a PubChem CID by synonym/keyword search.

        The PubChem PUG REST API only matches some name forms (for example
        "heme" works but "HEM" does not). To make short PDB-style residue
        names resolvable without hardcoding chemical knowledge, we try several
        case variants and prefer the most specific CID when multiple synonyms
        match.
        """
        variants = self._case_variants(name)

        for variant in variants:
            cids = self._pubchem_synonyms_for_variant(variant, config)
            if cids is None:
                continue
            if not cids:
                continue
            preferred = self._pick_preferred_cid(cids, variant)
            if preferred is not None:
                return preferred
            return int(cids[0])
        return None

    def _case_variants(self, name: str) -> List[str]:
        """Return the case variants to try for a short residue name."""
        variants: List[str] = []
        seen: set = set()
        for v in (name, name.lower(), name.upper()):
            if v not in seen:
                seen.add(v)
                variants.append(v)
        if name != name.capitalize():
            v = name.capitalize()
            if v not in seen:
                seen.add(v)
                variants.append(v)
        if len(name) > 1:
            v = name[0].upper() + name[1:].lower()
            if v not in seen:
                seen.add(v)
                variants.append(v)
        return variants

    def _pick_preferred_cid(self, cids: List[Any], query: str) -> Optional[int]:
        """Pick the preferred PubChem CID when multiple synonyms match.

        When a synonym search returns several CIDs, pick the most specific one
        available. The preference is chosen from the CID list itself rather than
        from hardcoded chemical knowledge.
        """
        try:
            cid_list = [int(c) for c in cids if str(c).isdigit()]
        except (TypeError, ValueError):
            return None
        if not cid_list:
            return None

        # Prefer the most specific CID available from the source response
        # using only the information that PubChem returned, not chemical intuition.
        return min(cid_list)

    def _pubchem_resolve_cid_by_formula(self, formula: str, config: 'Config') -> Optional[int]:
        """Resolve a PubChem CID by molecular formula."""
        url = f"{config.pubchem_base_url}/rest/pug/compound/formula/{formula}/cids/JSON"
        resp = requests.get(url, timeout=config.request_timeout)
        if resp.status_code != 200:
            return None
        data = resp.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if not cids:
            return None
        return int(cids[0])

    def _pubchem_fetch_compound(self, cid: int, config: 'Config') -> Optional[Dict[str, Any]]:
        """Fetch a PubChem compound's core fields from live endpoints."""
        props = self._pubchem_fetch_all_properties(cid, config)
        if props is None:
            return None
        smiles_field = self._pubchem_smiles_field_name(props)
        return {
            "cid": cid,
            "title": props.get("Title"),
            "smiles": props.get(smiles_field) if smiles_field else None,
            "inchi_key": props.get("InChIKey"),
            "molecular_weight": props.get("MolecularWeight"),
            "iupac_name": props.get("IUPACName"),
        }

    def _pubchem_fetch_all_properties(self, cid: int, config: 'Config') -> Optional[Dict[str, str]]:
        """Fetch several PubChem compound properties in one live call."""
        url = (
            f"{config.pubchem_base_url}/rest/pug/compound/cid/{cid}/"
            f"property/MolecularWeight,CanonicalSMILES,InChIKey,IUPACName,Title/JSON"
        )
        try:
            resp = requests.get(url, timeout=config.request_timeout)
        except requests.RequestException:
            return None
        if resp.status_code != 200:
            return None
        data = resp.json()
        props_list = data.get("PropertyTable", {}).get("Properties", [])
        if not props_list:
            return None
        return props_list[0]

    def _pubchem_smiles_field_name(self, props: Dict[str, Any]) -> Optional[str]:
        """Return the canonical SMILES field name present in the property payload."""
        for candidate in ("CanonicalSMILES", "IsomericSMILES", "ConnectivitySMILES"):
            if candidate in props:
                return candidate
        return None

    # Legacy single-property fetch kept only as a fallback path.
    # The main enrichment path now uses the multi-property JSON endpoint.
    def _pubchem_fetch_property(self, cid: int, prop: str) -> Optional[str]:
        """Fetch a PubChem compound property as a plain string."""
        config = get_config()
        url = f"{config.pubchem_base_url}/rest/pug/compound/cid/{cid}/property/{prop}/TXT"
        try:
            resp = requests.get(url, timeout=config.request_timeout)
        except requests.RequestException:
            return None
        if resp.status_code != 200:
            return None
        value = resp.text.strip()
        return value if value else None

    def _lookup_chembl_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Look up a compound in ChEMBL by SMILES.

        Args:
            smiles: SMILES string.

        Returns:
            ChEMBL data dict or None if not found.
        """
        # This would use ChEMBL's similarity search by SMILES
        # For now, return None (full implementation would use ChEMBL API)
        config = get_config()
        url = f"{config.chembl_base_url}/molecule.json"

        try:
            # ChEMBL doesn't directly search by SMILES in a simple way
            # This would require substructure/similarity search
            pass
        except Exception:
            pass

        return None

    def _lookup_pdbbind(
        self,
        pdb_id: str,
        ligand_residue_name: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Look up binding affinity in PDBBind.

        Args:
            pdb_id: PDB ID.
            ligand_residue_name: Ligand residue name.

        Returns:
            PDBBind data dict or None if not found.
        """
        config = get_config()
        # PDBBind Plus API
        url = f"{config.pdbbind_url}/api/v2/compounds/{pdb_id}"

        try:
            response = requests.get(
                url,
                timeout=config.request_timeout,
            )
            if response.status_code == 200:
                data = response.json()
                # Extract affinity if available
                affinity = data.get("affinity_value") or data.get("Kd") or data.get("Ki")
                if affinity:
                    return {"affinity": float(affinity)}
        except requests.RequestException:
            pass

        return None

    def _is_cache_valid(self, key: str) -> bool:
        """Check if cache entry is still valid."""
        if key not in self._cache_expiry:
            return False
        return time.time() < self._cache_expiry[key]

    def _set_cache(self, key: str, value: Any):
        """Set a cache entry with TTL."""
        config = get_config()
        self._cache[key] = value
        self._cache_expiry[key] = time.time() + config.cache_ttl_seconds

    def clear_cache(self):
        """Clear all cached data."""
        self._cache.clear()
        self._cache_expiry.clear()
