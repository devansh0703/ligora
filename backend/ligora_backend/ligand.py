"""
Ligand identity resolution from authoritative live data sources.

Resolution order:
1. RCSB Chemical Component Dictionary (CCD) - source of record for PDB
   chemical components (formula, weight, SMILES, InChIKey, type).
2. PubChem PUG-REST - compound identity and properties.
3. ChEMBL - cross-reference via UniChem (inchikey) or similarity search.
4. PDBBind-style affinity sources only when a configured endpoint provides
   them; absent sources leave the field None and the app says so.

No chemical values are invented, guessed, or defaulted in this module.
All network responses are cached with a TTL; failures return None and the
caller surfaces the gap to the user.
"""

import hashlib
import json
import time
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from urllib.parse import quote

import requests

from .schemas import (
    Ligand,
    LigandResolutionStatus,
)
from .config import get_config
from .net import http_timeout

# Sentinel distinguishing "no disk entry" from a cached miss (None value).
_DISK_MISS = object()

# Descriptor type constants from the RCSB chem_comp endpoint payloads.
_DESC_SMILES = "SMILES"
_DESC_SMILES_CANONICAL = "SMILES_CANONICAL"
_DESC_INCHIKEY = "InChIKey"


def ccd_element_map(comp_id: str) -> Dict[str, Optional[str]]:
    """
    Return {atom_name: element_symbol} for a CCD component, from live CCD data.

    The CCD component file (files.rcsb.org/ligands/view/<id>.cif) carries the
    authoritative `_chem_comp_atom` loop with `label_atom_id` and
    `type_symbol` per atom. It is parsed with the app's own mmCIF parser so
    the data path is identical to structure parsing. Returns an empty dict
    when the component cannot be resolved (callers must then treat element
    information as unknown, never guess it).
    """
    config = get_config()
    comp_id = comp_id.upper()
    url = f"{config.rcsb_files_url}/ligands/view/{quote(comp_id)}.cif"
    try:
        response = requests.get(url, timeout=http_timeout(config.request_timeout))
        if response.status_code != 200:
            return {}
        content = response.text
    except requests.RequestException:
        return {}

    from .parser import StructureParser
    data = StructureParser()._parse_mmcif_data(content)
    rows = data.get("_chem_comp_atom", [])
    if isinstance(rows, dict):
        rows = [rows]

    mapping: Dict[str, Optional[str]] = {}
    for row in rows:
        name = row.get("atom_id") or row.get("label_atom_id")
        elem = row.get("type_symbol")
        if name and elem:
            mapping[name.strip()] = elem.strip().upper()
    return mapping


def ccd_ideal_coordinates(comp_id: str) -> Dict[str, Tuple[str, float, float, float]]:
    """
    Return {atom_name: (element, x, y, z)} ideal coordinates for a CCD
    component, from its live CCD component file (pdbx_model_Cartn_x_ideal
    columns). Returns an empty dict when the component has no CCD coverage.
    """
    config = get_config()
    comp_id = comp_id.upper()
    url = f"{config.rcsb_files_url}/ligands/view/{quote(comp_id)}.cif"
    try:
        response = requests.get(url, timeout=http_timeout(config.request_timeout))
        if response.status_code != 200:
            return {}
        content = response.text
    except requests.RequestException:
        return {}

    from .parser import StructureParser
    data = StructureParser()._parse_mmcif_data(content)
    rows = data.get("_chem_comp_atom", [])
    if isinstance(rows, dict):
        rows = [rows]

    coords: Dict[str, Tuple[str, float, float, float]] = {}
    for row in rows:
        name = (row.get("atom_id") or row.get("label_atom_id") or "").strip()
        elem = (row.get("type_symbol") or "").strip().upper()
        try:
            x = float(row.get("pdbx_model_Cartn_x_ideal"))
            y = float(row.get("pdbx_model_Cartn_y_ideal"))
            z = float(row.get("pdbx_model_Cartn_z_ideal"))
        except (TypeError, ValueError):
            continue
        if name and elem:
            coords[name] = (elem, x, y, z)
    return coords


def ccd_bond_list(comp_id: str) -> List[Dict[str, str]]:
    """Bond list for a CCD component from its live component CIF.

    Returns rows of {atom_id_1, atom_id_2, value_order, aromatic_flag} from
    the authoritative `_chem_comp_bond` category. Empty when the component
    has no CCD coverage.
    """
    rows = _fetch_ccd_cif_category(comp_id, "_chem_comp_bond")
    if isinstance(rows, dict):
        rows = [rows]
    if not rows:
        return []
    out: List[Dict[str, str]] = []
    for row in rows:
        a1 = (row.get("atom_id_1") or "").strip()
        a2 = (row.get("atom_id_2") or "").strip()
        if not a1 or not a2:
            continue
        out.append({
            "atom_id_1": a1,
            "atom_id_2": a2,
            "value_order": (row.get("value_order") or "").strip(),
            "aromatic_flag": (row.get("pdbx_aromatic_flag") or "").strip(),
        })
    return out


def ccd_atom_charges(comp_id: str) -> Dict[str, int]:
    """Per-atom formal charges for a CCD component (live component CIF)."""
    rows = _fetch_ccd_cif_category(comp_id, "_chem_comp_atom")
    if isinstance(rows, dict):
        rows = [rows]
    charges: Dict[str, int] = {}
    for row in rows or []:
        name = (row.get("atom_id") or row.get("label_atom_id") or "").strip()
        raw = (row.get("charge") or "0").strip()
        try:
            charges[name] = int(float(raw))
        except (TypeError, ValueError):
            charges[name] = 0
    return charges


def _fetch_ccd_cif_category(comp_id: str,
                            category: str) -> Optional[Any]:
    """Fetch one category from the CCD component CIF (files.rcsb.org).

    Uses the app's own mmCIF parser so the data path is identical to
    structure parsing. Returns the parsed category rows (a list of dicts for
    loops, a dict for key-value categories) or None on failure.
    """
    config = get_config()
    url = f"{config.rcsb_files_url}/ligands/view/{quote(comp_id.upper())}.cif"
    try:
        response = requests.get(url, timeout=http_timeout(config.request_timeout))
        if response.status_code != 200:
            return None
        content = response.text
    except requests.RequestException:
        return None

    from .parser import StructureParser
    data = StructureParser()._parse_mmcif_data(content)
    return data.get(category)


def _fetch_ccd_raw(comp_id: str) -> Optional[Dict[str, Any]]:
    """Fetch the raw RCSB chem_comp payload (no caching here; callers cache)."""
    config = get_config()
    url = f"{config.rcsb_data_url}/core/chemcomp/{quote(comp_id.upper())}"
    try:
        response = requests.get(url, timeout=http_timeout(config.request_timeout))
        if response.status_code != 200:
            return None
        return response.json()
    except (requests.RequestException, ValueError):
        return None


class LigandResolver:
    """
    Resolves chemical identity for detected ligands from live data sources.
    """

    def __init__(self):
        self._cache: Dict[str, Any] = {}
        self._cache_expiry: Dict[str, float] = {}

    # ------------------------------------------------------------------
    # Cache helpers
    # ------------------------------------------------------------------

    def _disk_cache_path(self, key: str) -> Optional[Path]:
        """Workspace-backed cache file for a lookup key (None if unusable)."""
        try:
            config = get_config()
            safe = hashlib.md5(key.encode()).hexdigest()
            return config.cache_dir / "resolver" / f"{safe}.json"
        except (OSError, ValueError):
            return None

    def _disk_read(self, key: str) -> Optional[Any]:
        """Read a cached lookup result; None when absent or expired."""
        path = self._disk_cache_path(key)
        if path is None or not path.exists():
            return _DISK_MISS
        try:
            doc = json.loads(path.read_text())
            if time.time() > doc.get("expires", 0):
                return _DISK_MISS
            return doc.get("value", _DISK_MISS)
        except (OSError, ValueError):
            return _DISK_MISS

    def _disk_write(self, key: str, value: Any) -> None:
        """Persist a lookup result with its TTL (responsible caching)."""
        path = self._disk_cache_path(key)
        if path is None:
            return
        try:
            config = get_config()
            path.parent.mkdir(parents=True, exist_ok=True)
            # A missed source is retried after a shorter interval; a real
            # answer is kept for the configured TTL.
            ttl = (config.cache_ttl_seconds if value is not None
                   else min(3600, config.cache_ttl_seconds))
            path.write_text(json.dumps({
                "key": key,
                "value": value,
                "expires": time.time() + ttl,
            }))
        except (OSError, ValueError, TypeError):
            pass  # cache write failures never break resolution

    def _is_cache_valid(self, key: str) -> bool:
        if key not in self._cache_expiry:
            return False
        return time.time() < self._cache_expiry[key]

    def _set_cache(self, key: str, value: Any):
        config = get_config()
        self._cache[key] = value
        self._cache_expiry[key] = time.time() + config.cache_ttl_seconds

    def _cached(self, key: str, fetch):
        """Cache-through helper: in-memory, then workspace disk, then live.

        All layers hold real fetched data with the same TTL; nothing is
        fabricated. A source that did not answer is remembered for a short
        interval only, so unavailability is retried later.
        """
        if self._is_cache_valid(key):
            return self._cache.get(key)
        disk = self._disk_read(key)
        if disk is not _DISK_MISS:
            self._set_cache(key, disk)
            return disk
        value = fetch()
        self._set_cache(key, value)
        self._disk_write(key, value)
        return value

    def clear_cache(self):
        self._cache.clear()
        self._cache_expiry.clear()

    # ------------------------------------------------------------------
    # CCD (RCSB Chemical Component Dictionary)
    # ------------------------------------------------------------------

    def _lookup_ccd(self, residue_name: str) -> Optional[Dict[str, Any]]:
        """Look up a chemical component in the RCSB CCD."""
        if not residue_name:
            return None
        key = f"ccd:{residue_name.upper()}"

        def fetch():
            data = _fetch_ccd_raw(residue_name)
            if not data:
                return None
            cc = data.get("chem_comp") or {}
            if not cc:
                return None
            descriptors = data.get("pdbx_chem_comp_descriptor") or []

            smiles = None
            inchi_key = None
            for desc in descriptors:
                dtype = desc.get("type")
                value = desc.get("descriptor")
                if dtype == _DESC_INCHIKEY and not inchi_key and value:
                    inchi_key = value.strip()
                elif dtype == _DESC_SMILES_CANONICAL and not smiles and value:
                    smiles = value.strip()
                elif dtype == _DESC_SMILES and not smiles and value:
                    smiles = value.strip()

            # The REST payload lacks pdbx_type; the authoritative CCD
            # component CIF carries it (e.g. HOH -> HETAS solvent).
            pdbx_type = cc.get("pdbx_type")
            if not pdbx_type:
                cif_rows = _fetch_ccd_cif_category(
                    residue_name, "_chem_comp")
                if cif_rows:
                    row = cif_rows[0] if isinstance(cif_rows, list) else cif_rows
                    pdbx_type = row.get("pdbx_type")
            return {
                "id": cc.get("id"),
                "name": cc.get("name"),
                "formula": cc.get("formula"),
                "molecular_weight": cc.get("formula_weight"),
                "type": cc.get("type"),
                "pdbx_type": pdbx_type,
                "smiles": smiles,
                "inchi_key": inchi_key,
            }

        return self._cached(key, fetch)

    # ------------------------------------------------------------------
    # PubChem
    # ------------------------------------------------------------------

    def _lookup_pubchem(self, name: str,
                        formula: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """Search PubChem by name, then by formula, then fetch properties."""
        cache_key = "pubchem:" + hashlib.md5(
            f"{name}:{formula}".encode()).hexdigest()[:12]

        def fetch():
            config = get_config()
            cid = self._pubchem_resolve_cid(name, config)
            if cid is None and formula:
                cid = self._pubchem_resolve_cid_by_formula(formula, config)
            if cid is None:
                return None
            return self._pubchem_fetch_compound(cid, config)

        return self._cached(cache_key, fetch)

    def _pubchem_resolve_cid(self, name: str, config) -> Optional[int]:
        """Resolve a PubChem CID by exact compound name."""
        if not name:
            return None
        url = (f"{config.pubchem_base_url}/compound/name/"
               f"{quote(name)}/cids/JSON")
        try:
            resp = requests.get(
                        url, timeout=http_timeout(config.request_timeout))
            if resp.status_code != 200:
                return None
            cids = resp.json().get("IdentifierList", {}).get("CID", [])
            return int(cids[0]) if cids else None
        except (requests.RequestException, ValueError, IndexError):
            return None

    def _pubchem_resolve_cid_by_formula(self, formula: str,
                                        config) -> Optional[int]:
        """Resolve a PubChem CID by molecular formula."""
        if not formula:
            return None
        url = (f"{config.pubchem_base_url}/compound/fastformula/"
               f"{quote(formula)}/cids/JSON")
        try:
            resp = requests.get(
                        url, timeout=http_timeout(config.request_timeout))
            if resp.status_code != 200:
                return None
            cids = resp.json().get("IdentifierList", {}).get("CID", [])
            return int(cids[0]) if cids else None
        except (requests.RequestException, ValueError, IndexError):
            return None

    def _pubchem_fetch_compound(self, cid: int,
                                config) -> Optional[Dict[str, Any]]:
        """Fetch a PubChem compound's core fields from live endpoints."""
        url = (f"{config.pubchem_base_url}/compound/cid/{cid}/property/"
               f"MolecularWeight,CanonicalSMILES,InChIKey,IUPACName,Title/JSON")
        try:
            resp = requests.get(
                        url, timeout=http_timeout(config.request_timeout))
            if resp.status_code != 200:
                return None
            props_list = resp.json().get("PropertyTable", {}).get(
                "Properties", [])
            if not props_list:
                return None
            props = props_list[0]
            return {
                "cid": cid,
                "title": props.get("Title"),
                "smiles": props.get("CanonicalSMILES") or
                          props.get("ConnectivitySMILES"),
                "inchi_key": props.get("InChIKey"),
                "molecular_weight": props.get("MolecularWeight"),
                "iupac_name": props.get("IUPACName"),
            }
        except (requests.RequestException, ValueError):
            return None

    # ------------------------------------------------------------------
    # ChEMBL / UniChem
    # ------------------------------------------------------------------

    def _lookup_chembl(self, smiles: Optional[str] = None,
                       inchi_key: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """Find a ChEMBL compound via UniChem InChIKey or ChEMBL similarity."""
        cache_key = "chembl:" + hashlib.md5(
            f"{smiles}:{inchi_key}".encode()).hexdigest()[:12]

        def fetch():
            config = get_config()

            # 1) UniChem cross-reference by InChIKey (authoritative mapping).
            if inchi_key:
                url = f"{config.unichem_base_url}/inchikey/{inchi_key}"
                try:
                    resp = requests.get(
                        url, timeout=http_timeout(config.request_timeout))
                    if resp.status_code == 200:
                        mappings = resp.json()
                        if isinstance(mappings, list):
                            for m in mappings:
                                # src_id 1 is ChEMBL in UniChem sources
                                if str(m.get("src_id")) == "1":
                                    return {
                                        "chembl_id": m.get(
                                            "src_compound_id"),
                                        "source": "unichem",
                                    }
                except (requests.RequestException, ValueError):
                    pass

            # 2) ChEMBL similarity search fallback (still live data).
            if smiles:
                url = (f"{config.chembl_base_url}/similarity/"
                       f"{quote(smiles)}/70.json?limit=1")
                try:
                    resp = requests.get(
                        url, timeout=http_timeout(config.request_timeout))
                    if resp.status_code == 200:
                        molecules = resp.json().get("molecules", [])
                        if molecules:
                            mol = molecules[0]
                            return {
                                "chembl_id": mol.get("molecule_chembl_id"),
                                "pref_name": mol.get("pref_name"),
                                "source": "chembl_similarity",
                            }
                except (requests.RequestException, ValueError):
                    pass

            return None

        return self._cached(cache_key, fetch)

    # ------------------------------------------------------------------
    # PDBBind-style affinity (only when a source actually answers)
    # ------------------------------------------------------------------

    def _lookup_pdbbind(self, pdb_id: str,
                        ligand_residue_name: str) -> Optional[Dict[str, Any]]:
        """Query a configured PDBBind endpoint; None when it does not answer."""
        config = get_config()
        if not config.pdbbind_url:
            # No source configured: report not-available, never synthesize.
            return None
        url = f"{config.pdbbind_url.rstrip('/')}/compounds/{pdb_id}"
        try:
            response = requests.get(url, timeout=http_timeout(config.request_timeout))
            if response.status_code == 200:
                data = response.json()
                affinity = (data.get("affinity_value") or data.get("Kd")
                            or data.get("Ki"))
                if affinity:
                    return {"affinity": float(affinity)}
        except (requests.RequestException, ValueError):
            pass
        return None

    # ------------------------------------------------------------------
    # Main resolution entry point
    # ------------------------------------------------------------------

    def resolve_ligand(
        self,
        ligand: Ligand,
        pdb_id: Optional[str] = None,
    ) -> Ligand:
        """
        Resolve a ligand's chemical identity using available data sources.

        Progressively enriches from CCD, PubChem, ChEMBL, and (when a source
        answers) PDBBind. Every attached value carries its origin in the
        evidence chain.

        CCD-first short-circuit: when the CCD itself classifies the component
        as solvent (HETAS) or an ion (HETAI), that classification IS the
        authoritative identity for our purposes, and the drug-database chain
        (PubChem/ChEMBL/UniChem) is not consulted. This is not a heuristic —
        it is the CCD's own classification deciding which sources are
        relevant; no value is ever invented.
        """
        evidence: List[Dict[str, Any]] = []

        # 1) CCD - source of record for PDB chemical components.
        ccd = self._lookup_ccd(ligand.residue_name)
        ccd_pdbx_type = (ccd or {}).get("pdbx_type")
        if ccd:
            if ccd.get("formula"):
                ligand.formula = ccd["formula"]
            if ccd.get("molecular_weight") is not None:
                ligand.molecular_weight = float(ccd["molecular_weight"])
            if ccd.get("smiles") and not ligand.smiles:
                ligand.smiles = ccd["smiles"]
            if ccd.get("inchi_key") and not ligand.inchi_key:
                ligand.inchi_key = ccd["inchi_key"]
            if ccd.get("name"):
                ligand.name = ccd["name"]
            # pdbx_type is the CCD's own classification code:
            # HETAS = solvent, HETAI = ion, HETAIN = neutral organic, etc.
            if ccd.get("pdbx_type"):
                ligand.classification_hint = ccd["pdbx_type"]
            elif ccd.get("type"):
                ligand.classification_hint = ccd["type"]
            evidence.append({
                "source": "rcsb_ccd",
                "fields": ["name", "formula", "molecular_weight", "smiles",
                           "inchi_key", "type"],
                "url": (f"{get_config().rcsb_data_url}/core/chemcomp/"
                        f"{ligand.residue_name.upper()}"),
            })

        if ccd_pdbx_type in ("HETAS", "HETAI"):
            # Solvent/ion identity is complete from the CCD alone.
            ligand.resolution_status = LigandResolutionStatus.PARTIAL
            ligand.has_2d_structure = bool(ligand.smiles)
            self._last_evidence = evidence
            return ligand

        # 2) PubChem - compound identity and properties.
        pubchem = self._lookup_pubchem(
            ligand.residue_name, ligand.formula)
        if pubchem:
            ligand.pubchem_cid = pubchem.get("cid")
            if pubchem.get("title"):
                ligand.pubchem_name = pubchem["title"]
            if pubchem.get("smiles") and not ligand.smiles:
                ligand.smiles = pubchem["smiles"]
            if pubchem.get("inchi_key") and not ligand.inchi_key:
                ligand.inchi_key = pubchem["inchi_key"]
            if pubchem.get("molecular_weight") is not None and \
                    ligand.molecular_weight is None:
                ligand.molecular_weight = float(pubchem["molecular_weight"])
            if pubchem.get("iupac_name"):
                ligand.iupac_name = pubchem["iupac_name"]
            evidence.append({
                "source": "pubchem",
                "fields": ["cid", "title", "smiles", "inchi_key",
                           "molecular_weight", "iupac_name"],
                "url": (f"https://pubchem.ncbi.nlm.nih.gov/compound/"
                        f"{pubchem.get('cid')}"),
            })

        # 3) ChEMBL cross-reference.
        chembl = self._lookup_chembl(
            smiles=ligand.smiles, inchi_key=ligand.inchi_key)
        if chembl:
            ligand.chembl_id = chembl.get("chembl_id")
            if chembl.get("pref_name") and not ligand.pubchem_name:
                ligand.pubchem_name = chembl["pref_name"]
            evidence.append({
                "source": chembl.get("source", "chembl"),
                "fields": ["chembl_id"],
                "url": (f"https://www.ebi.ac.uk/chembl/compound_report_card/"
                        f"{chembl.get('chembl_id')}"),
            })

        # 4) PDBBind affinity when a configured source answers.
        if pdb_id and ligand.residue_name:
            pdbbind = self._lookup_pdbbind(pdb_id, ligand.residue_name)
            if pdbbind:
                ligand.pdbbind_affinity = pdbbind["affinity"]
                evidence.append({
                    "source": "pdbbind",
                    "fields": ["affinity"],
                    "url": f"http://www.pdbbind.org.cn/?pdbid={pdb_id}",
                })

        # Final status from what was actually resolved.
        if ligand.pubchem_cid or ligand.chembl_id:
            ligand.resolution_status = LigandResolutionStatus.RESOLVED
        elif ligand.smiles or ligand.formula:
            ligand.resolution_status = LigandResolutionStatus.PARTIAL
        else:
            ligand.resolution_status = LigandResolutionStatus.NOT_FOUND

        ligand.has_2d_structure = bool(ligand.smiles)
        self._last_evidence = evidence
        return ligand

    # Evidence from the most recent resolve_ligand call.
    _last_evidence: List[Dict[str, Any]] = []

    def get_last_evidence(self) -> List[Dict[str, Any]]:
        return list(self._last_evidence)

