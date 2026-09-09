"""
Enrichment client for live data lookup from external open sources.

Integrates:
- RCSB Data API (structure metadata, CCD chemical components)
- PubChem PUG-REST (compound identity, properties, bioactivity)
- ChEMBL REST (compound context, bioactivity)
- UniChem (identifier cross-mapping by InChIKey)

Every value comes from a live source response; nothing is synthesized
locally. Failures return None / empty and the UI reports the gap.
"""

import hashlib
import time
from typing import Optional, List, Dict, Any
from urllib.parse import quote

import requests

from .schemas import EvidenceItem
from .config import get_config
from .net import http_timeout


class EnrichmentClient:
    """Client for fetching enrichment data from external APIs."""

    def __init__(self):
        self._cache: Dict[str, Any] = {}
        self._cache_expiry: Dict[str, float] = {}
        self._session = requests.Session()
        self._session.headers.update({
            'User-Agent': 'Ligora/0.1.0 (open-source molecular workstation)',
        })

    # ------------------------------------------------------------------
    # Cache helpers
    # ------------------------------------------------------------------

    def _is_cache_valid(self, key: str) -> bool:
        if key not in self._cache_expiry:
            return False
        return time.time() < self._cache_expiry[key]

    def _cached(self, key: str, fetch):
        if self._is_cache_valid(key):
            return self._cache.get(key)
        value = fetch()
        config = get_config()
        self._cache[key] = value
        self._cache_expiry[key] = time.time() + config.cache_ttl_seconds
        return value

    def clear_cache(self):
        self._cache.clear()
        self._cache_expiry.clear()

    # ------------------------------------------------------------------
    # RCSB structure metadata
    # ------------------------------------------------------------------

    def get_structure_metadata(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """Entry metadata from the RCSB Data API."""
        pdb_id = pdb_id.strip().upper()
        cache_key = f"metadata:{pdb_id}"

        def fetch():
            config = get_config()
            url = f"{config.rcsb_data_url}/core/entry/{pdb_id}"
            try:
                response = self._session.get(
                    url, timeout=http_timeout(config.request_timeout))
                if response.status_code != 200:
                    return None
                data = response.json()

                rcsb_info = data.get('rcsb_entry_info', {})
                exptl = data.get('exptl', [{}])
                refine = data.get('refine', [{}])

                result = {
                    'pdb_id': pdb_id,
                    'title': (data.get('struct', {}) or {}).get('title'),
                    'experiment_type': (exptl[0] if exptl else {}).get(
                        'method'),
                    'deposition_date': (data.get('rcsb_accession_info', {})
                                        or {}).get('deposit_date'),
                    'release_date': (data.get('rcsb_accession_info', {})
                                     or {}).get('initial_release_date'),
                    'resolution': rcsb_info.get('resolution_combined'),
                    'experimental_method': rcsb_info.get(
                        'experimental_method'),
                    'entry_info': rcsb_info,
                }
                if not result['resolution'] and refine:
                    result['resolution'] = (refine[0] or {}).get(
                        'ls_d_res_high')
                return result
            except (requests.RequestException, ValueError):
                return None

        return self._cached(cache_key, fetch)

    # ------------------------------------------------------------------
    # CCD chemical component (classification: pdbx_type)
    # ------------------------------------------------------------------

    def get_compound_from_ccd(self, residue_name: str
                              ) -> Optional[Dict[str, Any]]:
        """Chemical component data from the RCSB CCD."""
        if not residue_name:
            return None
        residue_name = residue_name.upper()
        cache_key = f"ccd:{residue_name}"

        def fetch():
            config = get_config()
            url = (f"{config.rcsb_data_url}/core/chemcomp/"
                   f"{quote(residue_name)}")
            try:
                response = self._session.get(
                    url, timeout=http_timeout(config.request_timeout))
                if response.status_code != 200:
                    return None
                cc = (response.json().get('chem_comp') or {})
                if not cc:
                    return None
                return {
                    'id': cc.get('id'),
                    'name': cc.get('name'),
                    'formula': cc.get('formula'),
                    'formula_weight': cc.get('formula_weight'),
                    'type': cc.get('type'),
                    'pdbx_type': cc.get('pdbx_type'),
                }
            except (requests.RequestException, ValueError):
                return None

        return self._cached(cache_key, fetch)

    # ------------------------------------------------------------------
    # PubChem
    # ------------------------------------------------------------------

    PUBCHEM_PROPERTY_FIELDS = [
        'MolecularWeight', 'ConnectivitySMILES', 'IsomericSMILES',
        'IUPACName', 'InChIKey', 'XLogP', 'TPSA', 'HBondAcceptorCount',
        'HBondDonorCount', 'RotatableBondCount', 'ExactMass', 'Complexity',
    ]

    def _fetch_pubchem_properties(self, cid: int) -> Optional[Dict[str, Any]]:
        """Fetch PubChem properties in one PUG-REST call."""
        config = get_config()
        fields = ','.join(self.PUBCHEM_PROPERTY_FIELDS)
        url = (f"{config.pubchem_base_url}/compound/cid/{cid}/"
               f"property/{fields}/JSON")
        try:
            response = self._session.get(
                url, timeout=http_timeout(config.request_timeout))
            if response.status_code != 200:
                return None
            props = response.json().get('PropertyTable', {}).get(
                'Properties', [])
            return props[0] if props else None
        except (requests.RequestException, ValueError):
            return None

    def get_pubchem_info(self, compound_name: Optional[str] = None,
                         smiles: Optional[str] = None,
                         pubchem_cid: Optional[int] = None,
                         ) -> Optional[Dict[str, Any]]:
        """Compound info from PubChem by CID, name, or SMILES."""
        cache_key = "pubchem:" + hashlib.md5(
            f"{pubchem_cid}:{compound_name}:{smiles}".encode()
        ).hexdigest()[:16]

        def fetch():
            config = get_config()
            cid = pubchem_cid

            if cid is None and compound_name:
                url = (f"{config.pubchem_base_url}/compound/name/"
                       f"{quote(compound_name)}/cids/JSON")
                try:
                    resp = self._session.get(
                        url, timeout=http_timeout(config.request_timeout))
                    if resp.status_code == 200:
                        cids = resp.json().get('IdentifierList', {}).get(
                            'CID', [])
                        cid = cids[0] if cids else None
                except (requests.RequestException, ValueError):
                    return None

            # SMILES search via PubChem's fast identity endpoint
            if cid is None and smiles:
                url = (f"{config.pubchem_base_url}/compound/smiles/"
                       f"{quote(smiles)}/cids/JSON")
                try:
                    resp = self._session.get(
                        url, timeout=http_timeout(config.request_timeout))
                    if resp.status_code == 200:
                        cids = resp.json().get('IdentifierList', {}).get(
                            'CID', [])
                        cid = cids[0] if cids else None
                except (requests.RequestException, ValueError):
                    return None

            if cid is None:
                return None

            props = self._fetch_pubchem_properties(cid)
            if not props:
                return None
            props['cid'] = cid
            props['title'] = props.get('Title')
            props['smiles'] = (props.get('IsomericSMILES') or
                               props.get('ConnectivitySMILES'))
            return props

        return self._cached(cache_key, fetch)

    # ------------------------------------------------------------------
    # ChEMBL
    # ------------------------------------------------------------------

    def get_chembl_info(self, chembl_id: Optional[str] = None,
                        smiles: Optional[str] = None,
                        inchi_key: Optional[str] = None,
                        ) -> Optional[Dict[str, Any]]:
        """Compound info from ChEMBL by ID, or via UniChem InChIKey map."""
        if chembl_id:
            return self._get_chembl_by_id(chembl_id)

        if inchi_key:
            mapped = self._get_chembl_id_from_unichem(inchi_key)
            if mapped:
                return self._get_chembl_by_id(mapped)

        if smiles:
            return self._get_chembl_by_smiles(smiles)
        return None

    def _get_chembl_id_from_unichem(self, inchi_key: str) -> Optional[str]:
        """Map an InChIKey to a ChEMBL ID via the UniChem service.

        A transient timeout is retried once (same policy as the status
        checks) before being reported as unavailable — a slow EBI moment
        must not be misreported as 'no mapping exists'.
        """
        config = get_config()
        url = f"{config.unichem_base_url}/inchikey/{inchi_key}"
        for attempt in range(2):
            try:
                resp = self._session.get(
                    url, timeout=http_timeout(config.request_timeout))
                if resp.status_code != 200:
                    return None
                mappings = resp.json()
                if isinstance(mappings, list):
                    for m in mappings:
                        if str(m.get('src_id')) == '1':  # 1 = ChEMBL
                            return m.get('src_compound_id')
                return None
            except (requests.RequestException, ValueError):
                if attempt == 1:
                    return None
        return None

    def _get_chembl_by_id(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        config = get_config()
        url = f"{config.chembl_base_url}/molecule/{chembl_id}.json"
        try:
            response = self._session.get(
                url, timeout=http_timeout(config.request_timeout))
            if response.status_code != 200:
                return None
            molecule = (response.json().get('molecule') or {})
            props = (molecule.get('molecule_properties') or {})
            structs = (molecule.get('molecule_structures') or {})
            return {
                'chembl_id': chembl_id,
                'pref_name': molecule.get('pref_name'),
                'smiles': props.get('canonical_smiles') or structs.get(
                    'canonical_smiles'),
                'molecular_weight': props.get('full_mwt'),
                'formula': props.get('full_molformula'),
                'alogp': props.get('alogp'),
                'hba': props.get('hba'),
                'hbd': props.get('hbd'),
                'source': 'chembl',
            }
        except (requests.RequestException, ValueError):
            return None

    def _get_chembl_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """ChEMBL similarity search by SMILES."""
        config = get_config()
        url = (f"{config.chembl_base_url}/similarity/"
               f"{quote(smiles)}/70.json?limit=1")
        try:
            response = self._session.get(
                url, timeout=http_timeout(config.request_timeout))
            if response.status_code != 200:
                return None
            molecules = response.json().get('molecules', [])
            if not molecules:
                return None
            mol = molecules[0]
            props = (mol.get('molecule_properties') or {})
            return {
                'chembl_id': mol.get('molecule_chembl_id'),
                'pref_name': mol.get('pref_name'),
                'smiles': props.get('canonical_smiles'),
                'molecular_weight': props.get('full_mwt'),
                'source': 'chembl_similarity',
            }
        except (requests.RequestException, ValueError):
            return None

    def get_chembl_activity(self, chembl_id: str,
                            limit: int = 50) -> List[Dict[str, Any]]:
        """Bioactivity records for a ChEMBL compound."""
        config = get_config()
        url = (f"{config.chembl_base_url}/activity.json?"
               f"molecule_chembl_id={chembl_id}&limit={limit}")
        activities: List[Dict[str, Any]] = []
        try:
            response = self._session.get(
                url, timeout=http_timeout(config.request_timeout))
            if response.status_code != 200:
                return activities
            for activity in response.json().get('activities', []):
                activities.append({
                    'activity_id': activity.get('activity_id'),
                    'target_chembl_id': activity.get('target_chembl_id'),
                    'target_type': activity.get('target_type'),
                    'standard_type': activity.get('standard_type'),
                    'standard_value': activity.get('standard_value'),
                    'standard_units': activity.get('standard_units'),
                    'relation': activity.get('relation'),
                    'source': 'chembl',
                })
        except requests.RequestException:
            pass
        return activities

    # ------------------------------------------------------------------
    # Evidence assembly
    # ------------------------------------------------------------------

    def get_all_evidence(self, ligand: Dict[str, Any],
                         pdb_id: Optional[str] = None
                         ) -> List[EvidenceItem]:
        """Evidence items for what sources actually provided."""
        evidence_items = []

        if ligand.get('pubchem_cid'):
            evidence_items.append(EvidenceItem(
                source='pubchem',
                field='compound_identity',
                value={'cid': ligand['pubchem_cid']},
                url=(f"https://pubchem.ncbi.nlm.nih.gov/compound/"
                     f"{ligand['pubchem_cid']}"),
            ))

        if ligand.get('chembl_id'):
            evidence_items.append(EvidenceItem(
                source='chembl',
                field='compound_context',
                value={'chembl_id': ligand['chembl_id']},
                url=(f"https://www.ebi.ac.uk/chembl/compound_report_card/"
                     f"{ligand['chembl_id']}"),
            ))

        if ligand.get('pdbbind_affinity'):
            evidence_items.append(EvidenceItem(
                source='pdbbind',
                field='binding_affinity',
                value={'affinity': ligand['pdbbind_affinity']},
                url=f"http://www.pdbbind.org.cn/?pdbid={pdb_id or ''}",
            ))

        if ligand.get('classification_hint'):
            ligand_id = (ligand.get('residue_name', '')
                         or ligand.get('name', ''))
            evidence_items.append(EvidenceItem(
                source='rcsb_ccd',
                field='component_classification',
                value={'pdbx_type': ligand['classification_hint']},
                url=("https://www.rcsb.org/ligands/"
                     f"{quote(ligand_id)}"),
            ))

        return evidence_items

    # ------------------------------------------------------------------
    # Health checks (real requests)
    # ------------------------------------------------------------------

    def health_check(self) -> Dict[str, str]:
        """Check availability of all data sources with real requests."""
        sources = {}
        config = get_config()

        checks = {
            'rcsb': (f"{config.rcsb_data_url}/core/entry/1CRN",
                     lambda r: r.status_code == 200),
            'pubchem': (
                f"{config.pubchem_base_url}/compound/cid/702/property/"
                f"Title/JSON",
                lambda r: r.status_code == 200),
            'chembl': (f"{config.chembl_base_url}/status.json",
                       lambda r: r.status_code == 200),
            'unichem': (
                f"{config.unichem_base_url}/inchikey/"
                f"LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                lambda r: r.status_code == 200),
        }
        if config.pdbbind_url:
            checks['pdbbind'] = (
                config.pdbbind_url, lambda r: r.status_code < 500)

        for name, (url, ok) in checks.items():
            # One retry: several of these services (notably UniChem) are
            # intermittently slow; a single dropped read is a false negative.
            # A longer-than-usual per-attempt timeout avoids marking a busy
            #-but-alive service unavailable right after the app has issued
            # many other live calls. Every attempt still hits the live
            # service - nothing is cached or assumed.
            sources[name] = 'unavailable'
            for _attempt in range(2):
                try:
                    response = self._session.get(url, timeout=15)
                    if ok(response):
                        sources[name] = 'ok'
                        break
                    sources[name] = 'error'
                    break
                except requests.RequestException:
                    continue
        return sources
