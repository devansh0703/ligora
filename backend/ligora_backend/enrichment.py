"""
Enrichment client for live data lookup from external sources.

Integrates with:
- RCSB PDB (structure metadata, annotations, related structures)
- PubChem (compound identity, properties, bioactivity)
- ChEMBL (bioactivity data, targets, binding data)
- PDBBind (binding affinity data)
"""

import json
import hashlib
import time
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field

import requests

from .schemas import EvidenceItem, LigandResolutionStatus
from .config import get_config

Config = get_config


class EnrichmentClient:
    """
    Client for fetching enrichment data from external APIs.

    Provides:
    - RCSB CCD lookup for chemical component identity and classification
    - RCSB structure annotations
    - PubChem compound lookup
    - ChEMBL bioactivity data
    - PDBBind affinity lookup

    All chemical values come from those sources. No chemical data is
    hardcoded or synthesized locally.
    """

    def __init__(self):
        """Initialize the enrichment client."""
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._cache_expiry: Dict[str, float] = {}
        self._session = requests.Session()
        self._session.headers.update({
            'User-Agent': 'Ligora/0.1.0 (https://github.com/ligora)',
        })

    def get_compound_from_ccd(self, residue_name: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve chemical component data from the RCSB CCD.

        Tries the compound endpoint first, then the ligand endpoint.
        """
        compound = self.get_compound_from_ccd_v1_compound(residue_name)
        if compound is None:
            compound = self.get_compound_from_ccd_v1_ligand(residue_name)
        return compound

    def get_compound_from_ccd_v1_compound(self, residue_name: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve chemical component data from the RCSB /v1/compound/{id} endpoint.
        """
        cache_key = f"ccd:compound:{residue_name}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()
        url = f"{config.rcsb_base_url}/v1/compound/{residue_name}"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code != 200:
                return None
            data = response.json()
            return self._parse_ccd_compound_v1_compound(data)
        except requests.RequestException:
            return None

    def get_compound_from_ccd_v1_ligand(self, residue_name: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve chemical component data from the RCSB /v1/ligand/{id} endpoint.
        """
        cache_key = f"ccd:ligand:{residue_name}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()
        url = f"{config.rcsb_base_url}/v1/ligand/{residue_name}"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code != 200:
                return None
            data = response.json()
            return self._parse_ccd_compound_v1_ligand(data)
        except requests.RequestException:
            return None

    def _parse_ccd_compound_v1_compound(self, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Retrieve chemical component data from the RCSB CCD.

        Args:
            residue_name: The 3-letter CCD identifier (for example HEM, ATP, HOH).

        Returns:
            CCD data useful for ligand identity and classification, or None.
        """
        cache_key = f"ccd:{residue_name}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()
        url = f"{config.rcsb_base_url}/v1/compound/{residue_name}"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code != 200:
                return None
            data = response.json()
            return self._parse_ccd_compound(data)
        except requests.RequestException:
            return None

    def get_compound_from_ccd_v1_compound(self, residue_name: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve chemical component data from the RCSB /v1/compound/{id} endpoint.

        Args:
            residue_name: The 3-letter CCD identifier (for example HEM, ATP, HOH).

        Returns:
            CCD data useful for ligand identity and classification, or None.
        """
        cache_key = f"ccd:compound:{residue_name}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()
        url = f"{config.rcsb_base_url}/v1/compound/{residue_name}"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code != 200:
                return None
            data = response.json()
            return self._parse_ccd_compound_v1_compound(data)
        except requests.RequestException:
            return None

    def _parse_ccd_compound_v1_compound(self, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Extract authoritative CCD fields from the /v1/compound response shape."""
        compound = data.get('compound') or data
        if not isinstance(compound, dict):
            compound = data

        name = compound.get('name') or ''
        formula = compound.get('formula') or compound.get('chemical_formula') or compound.get('formula_string') or ''
        inchi_key = compound.get('inchi_key') or ''
        pdbx_type = compound.get('pdbx_type') or compound.get('type') or compound.get('chemical_type') or ''
        iupac_name = compound.get('iupac_name') or ''
        molecular_weight = compound.get('molecular_weight')
        smiles = compound.get('smiles') or ''

        if not smiles:
            rep = compound.get('representative') or {}
            if not smiles:
                smiles = rep.get('smiles') or ''
            if not inchi_key:
                inchi_key = rep.get('inchi_key') or ''

        result = {
            'name': name,
            'formula': formula,
            'inchi_key': inchi_key,
            'pdbx_type': pdbx_type,
            'iupac_name': iupac_name,
            'molecular_weight': molecular_weight,
            'smiles': smiles,
        }
        if not name and not formula and not pdbx_type:
            return None
        return result

    def _parse_ccd_compound_v1_ligand(self, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Extract authoritative CCD fields from the /v1/ligand response shape."""
        lig = data.get('ligand') or data
        if not isinstance(lig, dict):
            lig = data

        name = lig.get('name') or ''
        formula = lig.get('formula') or ''
        inchi_key = lig.get('inchi_key') or ''
        pdbx_type = lig.get('pdbx_type') or ''
        iupac_name = lig.get('iupac_name') or ''
        molecular_weight = lig.get('molecular_weight')
        smiles = lig.get('smiles') or ''

        if not smiles:
            rep = lig.get('representative') or {}
            if not smiles:
                smiles = rep.get('smiles') or ''
            if not inchi_key:
                inchi_key = rep.get('inchi_key') or ''

        result = {
            'name': name,
            'formula': formula,
            'inchi_key': inchi_key,
            'pdbx_type': pdbx_type,
            'iupac_name': iupac_name,
            'molecular_weight': molecular_weight,
            'smiles': smiles,
        }
        if not name and not formula and not pdbx_type:
            return None
        return result

    def _parse_ccd_compound_v1_ligand(self, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Extract authoritative CCD fields from the /v1/ligand response shape."""
        lig = data.get('ligand') or data
        if not isinstance(lig, dict):
            lig = data

        name = lig.get('name') or ''
        formula = lig.get('formula') or ''
        inchi_key = lig.get('inchi_key') or ''
        pdbx_type = lig.get('pdbx_type') or ''
        iupac_name = lig.get('iupac_name') or ''
        molecular_weight = lig.get('molecular_weight')
        smiles = lig.get('smiles') or ''

        if not smiles:
            rep = lig.get('representative') or {}
            if not smiles:
                smiles = rep.get('smiles') or ''
            if not inchi_key:
                inchi_key = rep.get('inchi_key') or ''

        result = {
            'name': name,
            'formula': formula,
            'inchi_key': inchi_key,
            'pdbx_type': pdbx_type,
            'iupac_name': iupac_name,
            'molecular_weight': molecular_weight,
            'smiles': smiles,
        }
        if not name and not formula and not pdbx_type:
            return None
        return result

    def _parse_ccd_compound(self, data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Extract authoritative CCD fields from a RCSB compound response.

        CCD is the source of record for chemical component identity in the PDB.
        This parser tolerates several response shapes returned by the live API
        and only keeps fields that come directly from that source.
        """
        compound = data.get('compound') or data
        if not isinstance(compound, dict):
            compound = data

        name = compound.get('name') or ''
        formula = compound.get('formula') or ''
        inchi_key = compound.get('inchi_key') or ''
        pdbx_type = compound.get('pdbx_type') or ''
        iupac_name = compound.get('iupac_name') or ''
        molecular_weight = compound.get('molecular_weight')
        smiles = compound.get('smiles') or ''

        # CCD compound responses sometimes embed the representation block
        # under different keys; try a couple of common shapes.
        if not smiles:
            rep = compound.get('representative') or {}
            if not smiles:
                smiles = rep.get('smiles') or ''
            if not inchi_key:
                inchi_key = rep.get('inchi_key') or ''

        # Alternate shapes observed in live responses
        if not formula:
            formula = compound.get('chemical_formula') or compound.get('formula_string') or ''
        if not pdbx_type:
            pdbx_type = compound.get('type') or compound.get('chemical_type') or ''

        result = {
            'name': name,
            'formula': formula,
            'inchi_key': inchi_key,
            'pdbx_type': pdbx_type,
            'iupac_name': iupac_name,
            'molecular_weight': molecular_weight,
            'smiles': smiles,
        }
        if not name and not formula and not pdbx_type:
            return None
        return result

    def get_structure_metadata(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """
        Fetch metadata for a PDB structure from RCSB.

        Args:
            pdb_id: PDB ID (e.g., "1ABC").

        Returns:
            Structure metadata or None if unavailable.
        """
        cache_key = f"metadata:{pdb_id}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()
        url = f"{config.rcsb_base_url}/v1/core/entry/{pdb_id}"

        try:
            response = self._session.get(
                url,
                timeout=config.request_timeout,
            )
            if response.status_code == 200:
                data = response.json()
                result = {
                    'pdb_id': pdb_id,
                    'title': data.get('title', ''),
                    'resolution': data.get('structure_determination_method',
                                          {}).get('resolution'),
                    'experiment_type': data.get('structure_determination_method',
                                                {}).get('method'),
                    'deposition_date': data.get('deposited'),
                    'release_date': data.get('released'),
                    ' macromolecule_type': data.get('macromolecule_type'),
                    'processing': data.get('assembly', {}).get('name'),
                    'organism': data.get('source_organism', {}).get('scientific_name'),
                    'mutations': data.get('macromolecule_details', {}).get('mutation'),
                }
                self._set_cache(cache_key, result)
                return result
        except requests.RequestException:
            pass

        return None

    def get_related_structures(
        self,
        pdb_id: str,
        by_ligand: bool = True,
        by_sequence: bool = False,
        by_3d: bool = False,
    ) -> List[Dict[str, Any]]:
        """
        Find related structures in the PDB.

        Args:
            pdb_id: PDB ID to find relatives for.
            by_ligand: Find structures with same ligand.
            by_sequence: Find structures with similar sequence.
            by_3d: Find structures with similar 3D structure.

        Returns:
            List of related structure info.
        """
        # Use RCSB Search API for related structures
        config = get_config()
        results = []

        if by_ligand:
            # Find structures with same ligands
            url = f"{config.rcsb_search_url}/search/v2/query"

            query = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": pdb_id,
                    }
                },
                "request": {
                    "return_type": "entry",
                    "request_mode": "interleaved",
                    "pagination": {"start": 0, "rows": 20},
                }
            }

            try:
                response = self._session.post(
                    url,
                    json=query,
                    timeout=config.request_timeout,
                )
                if response.status_code == 200:
                    data = response.json()
                    result_set = data.get('result_set', [])
                    for entry in result_set:
                        pdb_id_entry = entry.get('identifier')
                        if pdb_id_entry != pdb_id:
                            results.append({
                                'pdb_id': pdb_id_entry,
                                'similarity': entry.get('score', 0),
                                'type': 'same_ligand',
                            })
            except requests.RequestException:
                pass

        # For sequence/3D similarity, would use different query types
        # This is a simplified implementation

        return results[:20]

    def get_ligand_3d_info(self, ligand: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Get 3D structure information for a ligand.

        Args:
            ligand: Ligand data with CCD ID.

        Returns:
            3D structure info or None.
        """
        ccd_id = ligand.get('ccd_id') or ligand.get('residue_name')
        if not ccd_id:
            return None

        config = get_config()
        url = f"{config.rcsb_base_url}/v1/ligand/{ccd_id}/structure"

        try:
            response = self._session.get(
                url,
                timeout=config.request_timeout,
            )
            if response.status_code == 200:
                data = response.json()
                return {
                    'ccd_id': ccd_id,
                    'has_3d': data.get('has_3d_structure', False),
                    '3d_source': data.get('3d_source'),
                    '3d_url': data.get('3d_url'),
                }
        except requests.RequestException:
            pass

        return None

    def get_pubchem_info(
        self,
        compound_name: Optional[str] = None,
        smiles: Optional[str] = None,
        formula: Optional[str] = None,
        pubchem_cid: Optional[int] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Get compound information from PubChem.

        Args:
            compound_name: Compound name.
            smiles: SMILES string.
            formula: Chemical formula.
            pubchem_cid: PubChem CID (if known).

        Returns:
            PubChem compound info or None.
        """
        cache_key_parts = [
            str(pubchem_cid) if pubchem_cid else '',
            str(compound_name) if compound_name else '',
            str(smiles) if smiles else '',
            str(formula) if formula else '',
        ]
        cache_key = f"pubchem:{hashlib.md5(
            ':'.join(cache_key_parts).encode()
        ).hexdigest()[:16]}"
        if self._is_cache_valid(cache_key):
            return self._cache.get(cache_key)

        config = get_config()
        base_url = f"{config.pubchem_base_url}/rest"

        try:
            # If we have CID, fetch directly
            if pubchem_cid:
                url = f"{base_url}/compound/cid/{pubchem_cid}/property/" \
                      "IUPACName,MolecularWeight,CanonicalSMILES," \
                      "IsomericSMILES,InChIKey,Formula,XLogP,TSPA"
                response = self._session.get(url, timeout=config.request_timeout)
                if response.status_code == 200:
                    props = response.json().get('PropertyTable', {}).get('Properties', [])
                    if props:
                        result = {
                            'cid': pubchem_cid,
                            'iupac_name': props[0].get('IUPACName'),
                            'title': props[0].get('Title'),
                            'smiles': props[0].get('CanonicalSMILES'),
                            'isomeric_smiles': props[0].get('IsomericSMILES'),
                            'inchi_key': props[0].get('InChIKey'),
                            'formula': props[0].get('Formula'),
                            'molecular_weight': props[0].get('MolecularWeight'),
                            'xlogp': props[0].get('XLogP'),
                            'tpsa': props[0].get('TSPA'),
                            'source': 'pubchem',
                        }
                        self._set_cache(cache_key, result)
                        return result

            # Otherwise search by name
            if compound_name:
                search_url = f"{base_url}/compound/name/{compound_name.replace(' ', '%20')}"
                response = self._session.get(search_url, timeout=config.request_timeout)
                if response.status_code == 200:
                    cids = response.json().get('IdentifierList', {}).get('CID', [])
                    if cids:
                        return self.get_pubchem_info(pubchem_cid=cids[0])

            # Search by SMILES
            if smiles:
                # Use SMILES search in PubChem
                search_url = f"{base_url}/compound/smiles/{smiles.replace(' ', '%20')}"
                response = self._session.get(search_url, timeout=config.request_timeout)
                if response.status_code == 200:
                    cids = response.json().get('IdentifierList', {}).get('CID', [])
                    if cids:
                        return self.get_pubchem_info(pubchem_cid=cids[0])

        except requests.RequestException:
            pass

        return None

    def get_chembl_info(
        self,
        smiles: Optional[str] = None,
        pubchem_cid: Optional[int] = None,
        chembl_id: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Get compound information from ChEMBL.

        Args:
            smiles: SMILES string.
            pubchem_cid: PubChem CID for cross-reference.
            chembl_id: ChEMBL compound ID.

        Returns:
            ChEMBL compound info or None.
        """
        if chembl_id:
            return self._get_chembl_by_id(chembl_id)

        if smiles:
            return self._get_chembl_by_smiles(smiles)

        if pubchem_cid:
            return self._get_chembl_by_pubchem(pubchem_cid)

        return None

    def _get_chembl_by_id(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """Fetch ChEMBL compound by ID."""
        config = get_config()
        url = f"{config.chembl_base_url}/molecule/{chembl_id}.json"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code == 200:
                data = response.json()
                molecule = data.get('molecule', {})
                return {
                    'chembl_id': chembl_id,
                    'pref_name': molecule.get('pref_name'),
                    'smiles': molecule.get('molecule_properties', {}).get('canonical_smiles'),
                    'molecular_weight': molecule.get('molecule_properties', {}).get('full_mwt'),
                    'alogp': molecule.get('molecule_properties', {}).get('alogp'),
                    'hba': molecule.get('molecule_properties', {}).get('hba'),
                    'hbd': molecule.get('molecule_properties', {}).get('hbd'),
                    'source': 'chembl',
                }
        except requests.RequestException:
            pass

        return None

    def _get_chembl_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """Search ChEMBL by SMILES using similarity search."""
        config = get_config()
        url = f"{config.chembl_base_url}/similarity.json?smiles={smiles}&similarity_score=70"
        
        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code == 200:
                data = response.json()
                molecules = data.get('molecules', [])
                if molecules:
                    mol = molecules[0]
                    mol_props = mol.get('molecule_properties', {})
                    return {
                        'chembl_id': mol.get('molecule_chembl_id'),
                        'pref_name': mol.get('pref_name'),
                        'smiles': mol_props.get('canonical_smiles'),
                        'molecular_weight': mol_props.get('full_mwt'),
                        'alogp': mol_props.get('alogp'),
                        'hba': mol_props.get('hba'),
                        'hbd': mol_props.get('hbd'),
                        'source': 'chembl',
                    }
        except requests.RequestException:
            pass

        return None

    def _get_chembl_by_pubchem(self, pubchem_cid: int) -> Optional[Dict[str, Any]]:
        """Get ChEMBL ID from PubChem CID cross-reference."""
        config = get_config()
        url = f"{config.pubchem_base_url}/rest/compound/cid/{pubchem_cid}/xrefs"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code == 200:
                data = response.json()
                # Look for ChEMBL cross-references
                for xref_list in data.get('InformationList', {}).get('Information', [{}])[0].get('XrefList', []):
                    if xref_list.get('Name') == 'ChEMBL':
                        chembl_id = xref_list.get('Qualifier') or xref_list.get('ID')
                        if chembl_id:
                            return self._get_chembl_by_id(chembl_id)
        except requests.RequestException:
            pass

        return None

    def get_chembl_activity(
        self,
        chembl_id: str,
        target_type: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Get bioactivity data for a ChEMBL compound.

        Args:
            chembl_id: ChEMBL compound ID.
            target_type: Optional target type filter.

        Returns:
            List of activity records.
        """
        config = get_config()
        activities = []

        # Get activities for this molecule
        url = f"{config.chembl_base_url}/activity.json?molecule_chembl_id={chembl_id}&limit=100"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code == 200:
                data = response.json()
                activities_list = data.get('activities', [])
                for activity in activities_list[:50]:  # Limit to first 50
                    if target_type and target_type not in str(activity):
                        continue
                    activities.append({
                        'chembl_id': chembl_id,
                        'activity_id': activity.get('activity_id'),
                        'target_chembl_id': activity.get('target_chembl_id'),
                        'target_type': activity.get('target_type'),
                        'standard_type': activity.get('standard_type'),
                        'standard_value': activity.get('standard_value'),
                        'standard_units': activity.get('standard_units'),
                        'relation': activity.get('relation'),
                        'document_id': activity.get('document_id'),
                        'source': 'chembl',
                    })
        except requests.RequestException:
            pass

        return activities

    def get_pdbbind_info(
        self,
        pdb_id: str,
        ligand_residue_name: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Get binding affinity information from PDBBind.

        Args:
            pdb_id: PDB ID.
            ligand_residue_name: Ligand residue name.

        Returns:
            PDBBind entry info or None.
        """
        config = get_config()
        url = f"{config.pdbbind_url}/api/v2/compounds/{pdb_id}"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code == 200:
                data = response.json()
                # Extract affinity
                affinity = data.get('affinity_value') or data.get('Kd') or data.get('Ki')
                if affinity:
                    return {
                        'pdb_id': pdb_id,
                        'ligand': ligand_residue_name,
                        'affinity': float(affinity),
                        'affinity_type': data.get('affinity_type'),
                        'source': 'pdbbind',
                    }
        except requests.RequestException:
            pass

        return None

    def get_pubchem_bioactivity(
        self,
        pubchem_cid: int,
    ) -> List[Dict[str, Any]]:
        """
        Get bioactivity data from PubChem for a compound.

        Args:
            pubchem_cid: PubChem CID.

        Returns:
            List of bioactivity records.
        """
        config = get_config()
        activities = []

        url = f"{config.pubchem_base_url}/rest/pcassay/aid/{pubchem_cid}/cids/{pubchem_cid}/bioactivity"

        try:
            response = self._session.get(url, timeout=config.request_timeout)
            if response.status_code == 200:
                data = response.json()
                # Process bioactivity data
                pass  # Simplified
        except requests.RequestException:
            pass

        return activities

    def enrich_ligand(
        self,
        ligand: Dict[str, Any],
        pdb_id: Optional[str] = None,
        ccd_data: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Enrich a ligand with data from multiple sources.

        Chemical identity, formula, molecular weight, SMILES, and
        classification must come from live data sources (CCD, PubChem,
        ChEMBL, PDBBind), not from local heuristics.

        Args:
            ligand: Ligand data.
            pdb_id: Optional PDB ID.
            ccd_data: Optional RCSB CCD data for the ligand residue name.

        Returns:
            Enriched ligand data with evidence chain.
        """
        enriched = ligand.copy()
        evidence: List[Dict[str, Any]] = []

        # Attach authoritative CCD data first, because it is the source
        # of record for chemical component identity in the PDB.
        if ccd_data:
            ccd_fields = []
            if ccd_data.get('formula'):
                enriched['formula'] = ccd_data.get('formula')
                ccd_fields.append('formula')
            if ccd_data.get('molecular_weight') is not None:
                try:
                    enriched['molecular_weight'] = float(ccd_data.get('molecular_weight'))
                except (TypeError, ValueError):
                    pass
                ccd_fields.append('molecular_weight')
            if ccd_data.get('smiles'):
                enriched['smiles'] = ccd_data.get('smiles')
                ccd_fields.append('smiles')
            if ccd_data.get('inchi_key'):
                enriched['inchi_key'] = ccd_data.get('inchi_key')
                ccd_fields.append('inchi_key')
            if ccd_data.get('name'):
                enriched['name'] = ccd_data.get('name')
                ccd_fields.append('name')
            if ccd_data.get('pdbx_type'):
                enriched['classification_hint'] = ccd_data.get('pdbx_type')
                ccd_fields.append('pdbx_type')
            if ccd_data.get('iupac_name'):
                enriched['iupac_name'] = ccd_data.get('iupac_name')
                ccd_fields.append('iupac_name')
            if ccd_fields:
                evidence.append({
                    'source': 'rcsb_ccd',
                    'fields': ccd_fields,
                    'url': f"{config.rcsb_base_url}/v1/compound/{ligand.get('residue_name')}",
                })

        # Try PubChem by name/formula when no authoritative SMILES is available.
        if not enriched.get('smiles') and (
            ligand.get('name') or ligand.get('residue_name')
        ):
            pubchem_info = self.get_pubchem_info(
                compound_name=ligand.get('name') or ligand.get('residue_name'),
                formula=ligand.get('formula'),
            )
            if pubchem_info:
                enriched['pubchem_cid'] = pubchem_info.get('cid')
                enriched['pubchem_name'] = pubchem_info.get('title')
                if not enriched.get('smiles'):
                    enriched['smiles'] = pubchem_info.get('smiles')
                if not enriched.get('inchi_key'):
                    enriched['inchi_key'] = pubchem_info.get('inchi_key')
                if pubchem_info.get('molecular_weight') is not None and not enriched.get('molecular_weight'):
                    enriched['molecular_weight'] = pubchem_info.get('molecular_weight')
                if pubchem_info.get('xlogp') is not None:
                    enriched['xlogp'] = pubchem_info.get('xlogp')
                if pubchem_info.get('tpsa') is not None:
                    enriched['tpsa'] = pubchem_info.get('tpsa')
                if pubchem_info.get('iupac_name'):
                    enriched['iupac_name'] = pubchem_info.get('iupac_name')
                evidence.append({
                    'source': 'pubchem',
                    'fields': ['cid', 'title', 'smiles', 'inchi_key', 'molecular_weight'],
                    'url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_info.get('cid')}",
                })

        # Try ChEMBL for bioactivity context when a SMILES is available.
        if enriched.get('smiles'):
            chembl_info = self.get_chembl_info(smiles=enriched['smiles'])
            if chembl_info:
                enriched['chembl_id'] = chembl_info.get('chembl_id')
                if not enriched.get('iupac_name'):
                    enriched['iupac_name'] = chembl_info.get('pref_name')
                evidence.append({
                    'source': 'chembl',
                    'fields': ['chembl_id', 'pref_name', 'smiles'],
                    'url': f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_info.get('chembl_id')}",
                })

        # Try PDBBind for affinity context when a PDB ID is available.
        if pdb_id:
            pdbbind_info = self.get_pdbbind_info(
                pdb_id,
                ligand.get('residue_name') or ligand.get('name'),
            )
            if pdbbind_info:
                enriched['pdbbind_affinity'] = pdbbind_info.get('affinity')
                enriched['pdbbind_affinity_type'] = pdbbind_info.get('affinity_type')
                evidence.append({
                    'source': 'pdbbind',
                    'fields': ['affinity'],
                    'url': f"http://www.pdbbind.org.cn/?pdbid={pdb_id}",
                })

        enriched['evidence'] = evidence
        enriched['resolution_status'] = self._determine_resolution_status(enriched)

        return enriched

    def _determine_resolution_status(self, enriched: Dict[str, Any]) -> str:
        """Determine ligand resolution status."""
        if enriched.get('pubchem_cid') or enriched.get('chembl_id'):
            return LigandResolutionStatus.RESOLVED.value
        if enriched.get('smiles'):
            return LigandResolutionStatus.PARTIAL.value
        return LigandResolutionStatus.NOT_FOUND.value

    def get_all_evidence(
        self,
        ligand: Dict[str, Any],
        pdb_id: Optional[str] = None,
    ) -> List[EvidenceItem]:
        """
        Get all available evidence items for a ligand.

        Args:
            ligand: Ligand data.
            pdb_id: Optional PDB ID.

        Returns:
            List of evidence items.
        """
        evidence_items = []

        # PubChem
        if ligand.get('pubchem_cid'):
            evidence_items.append(EvidenceItem(
                source='pubchem',
                field='compound_identity',
                value={'cid': ligand['pubchem_cid']},
                url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{ligand['pubchem_cid']}",
            ))

        # ChEMBL
        if ligand.get('chembl_id'):
            evidence_items.append(EvidenceItem(
                source='chembl',
                field='bioactivity',
                value={'chembl_id': ligand['chembl_id']},
                url=f"https://www.ebi.ac.uk/chembl/compound_report_card/{ligand['chembl_id']}",
            ))

        # PDBBind
        if ligand.get('pdbbind_affinity'):
            evidence_items.append(EvidenceItem(
                source='pdbbind',
                field='binding_affinity',
                value={'affinity': ligand['pdbbind_affinity'],
                       'units': ligand.get('pdbbind_affinity_type', 'nM')},
                url=f"http://www.pdbbind.org.cn/?pdbid={pdb_id or ''}",
            ))

        return evidence_items

    def _is_cache_valid(self, key: str) -> bool:
        """Check if cache entry is valid."""
        if key not in self._cache_expiry:
            return False
        return time.time() < self._cache_expiry[key]

    def _set_cache(self, key: str, value: Any):
        """Set cache entry with TTL."""
        config = get_config()
        self._cache[key] = value
        self._cache_expiry[key] = time.time() + config.cache_ttl_seconds

    def clear_cache(self):
        """Clear all cached data."""
        self._cache.clear()
        self._cache_expiry.clear()

    def health_check(self) -> Dict[str, str]:
        """
        Check availability of all data sources.

        Returns:
            Dictionary of source status.
        """
        sources = {
            'rscb': 'unknown',
            'pubchem': 'unknown',
            'chembl': 'unknown',
            'pdbbind': 'unknown',
        }

        config = get_config()

        # Check RCSB
        try:
            response = self._session.get(
                f"{config.rcsb_base_url}/v1/health",
                timeout=5,
            )
            sources['rscb'] = 'ok' if response.status_code == 200 else 'error'
        except requests.RequestException:
            sources['rscb'] = 'unavailable'

        # Check PubChem
        try:
            response = self._session.get(
                f"{config.pubchem_base_url}/rest/health",
                timeout=5,
            )
            sources['pubchem'] = 'ok' if response.status_code == 200 else 'error'
        except requests.RequestException:
            sources['pubchem'] = 'unavailable'

        # Check ChEMBL
        try:
            response = self._session.get(
                f"{config.chembl_base_url}/health",
                timeout=5,
            )
            sources['chembl'] = 'ok' if response.status_code == 200 else 'error'
        except requests.RequestException:
            sources['chembl'] = 'unavailable'

        # Check PDBBind
        try:
            response = self._session.get(
                f"{config.pdbbind_url}/api/v2/health",
                timeout=5,
            )
            sources['pdbbind'] = 'ok' if response.status_code == 200 else 'error'
        except requests.RequestException:
            sources['pdbbind'] = 'unavailable'

        return sources
