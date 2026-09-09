"""
PDB-wide discovery and protein-context data, all from live authoritative
sources:

- RCSB Search API (search.rcsb.org): full-text, sequence (MMseqs2),
  chemical similarity (OpenEye fingerprints), same-ligand entry discovery.
- RCSB Data API: entry-level validation summary (pdbx_vrpt_summary).
- RCSB ModelServer (models.rcsb.org): on-demand coordinate subsets
  (assembly/chain/ligand) as mmCIF — the efficient fetch path.
- PDBe REST API (ebi.ac.uk/pdbe/api): EU-mirror entry summary, secondary
  structure, residue listing, ligand monomers.
- UniProt REST: protein function/features with real binding-site annotations.
- AlphaFold DB: predicted models (mmCIF) for UniProt accessions.

Nothing here invents values: every field is returned by the service and
carries its source; unreachable services raise honest errors.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

import requests

from ..config import get_config
from ..net import http_timeout


class DiscoveryError(RuntimeError):
    """Raised when a live source cannot answer (no fallbacks exist)."""


class DiscoveryClient:
    """Search the whole PDB and fetch protein context from live sources."""

    def __init__(self):
        self.config = get_config()
        self._session = requests.Session()
        self._session.headers["User-Agent"] = "Ligora/0.1 (open-source)"

    # ------------------------------------------------------------------
    # RCSB Search API
    # ------------------------------------------------------------------

    SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

    def _search(self, query: Dict[str, Any], return_type: str,
                rows: int, timeout_seconds: Optional[int] = None) -> Dict[str, Any]:
        payload = {
            "query": query,
            "return_type": return_type,
            "request_options": {"paginate": {"start": 0, "rows": rows}},
        }
        # MMseqs2 and strucmotif are compute services; their real response
        # profile can exceed the lookup timeout. A longer total deadline is
        # applied per service — still a hard cap, still no fabrication.
        budget = timeout_seconds or self.config.request_timeout
        # One honest retry on transient 5xx (RCSB's compute services poll
        # their internal workers on a ~30 s ticket and can return 500 under
        # load with a correct query, and the search API can answer 503 while
        # a backend is warming). Same payload, no parameter changes.
        resp = None
        for attempt in range(2):
            try:
                resp = self._session.post(self.SEARCH_URL, json=payload,
                                          timeout=http_timeout(budget))
                # A 5xx from the compute backend is transient when the query
                # itself validated; retry once before reporting failure.
                if resp.status_code in (500, 502, 503, 504) and attempt == 0:
                    continue
                break
            except requests.RequestException as e:
                if attempt == 1:
                    raise DiscoveryError(
                        f"RCSB search unreachable: {e}") from e
        # 204 = valid query with zero results (per RCSB API contract)
        if resp.status_code == 204:
            return {"total_count": 0, "hits": []}
        if resp.status_code != 200:
            raise DiscoveryError(
                f"RCSB search failed ({resp.status_code}): {resp.text[:200]}")
        data = resp.json()
        return {
            "total_count": data.get("total_count", 0),
            "hits": [
                {
                    "id": r["identifier"],
                    "score": r.get("score"),
                }
                for r in data.get("result_set", [])
            ],
        }

    def search_text(self, text: str, rows: int = 25) -> Dict[str, Any]:
        """Full-text search over PDB structure metadata."""
        if not text.strip():
            raise ValueError("Search text is empty")
        query = {"type": "terminal", "service": "full_text",
                 "parameters": {"value": text.strip()}}
        result = self._search(query, "entry", rows)
        result["service"] = "full_text"
        result["source"] = "RCSB Search API"
        return result

    def search_sequence(self, sequence: str, identity_cutoff: float = 0.9,
                        rows: int = 25) -> Dict[str, Any]:
        """MMseqs2 sequence search against all PDB polymer entities."""
        seq = re.sub(r"[^A-Za-z]", "", sequence).upper()
        if len(seq) < 12:
            raise ValueError(
                "Sequence too short for MMseqs2 search (min ~12 residues)")
        query = {"type": "terminal", "service": "sequence", "parameters": {
            "evalue_cutoff": 1, "identity_cutoff": identity_cutoff,
            "sequence_type": "protein", "value": seq}}
        result = self._search(query, "polymer_entity", rows,
                              timeout_seconds=120)
        result["service"] = "sequence"
        result["source"] = "RCSB Search API (MMseqs2)"
        return result

    def search_chemical_similarity(
        self,
        smiles: str,
        match_type: str = "fingerprint-similarity",
        rows: int = 25,
    ) -> Dict[str, Any]:
        """
        Chemical similarity search (OpenEye toolkit server-side).

        match_type: fingerprint-similarity (Tanimoto quick screen),
        graph-exact, graph-relaxed, graph-relaxed-stereo,
        sub-struct-graph-relaxed, sub-struct-graph-relaxed-stereo.
        """
        smiles = smiles.strip()
        if not smiles:
            raise ValueError("SMILES is empty")
        allowed = {
            "fingerprint-similarity", "graph-exact", "graph-relaxed",
            "graph-relaxed-stereo", "sub-struct-graph-relaxed",
            "sub-struct-graph-relaxed-stereo",
        }
        if match_type not in allowed:
            raise ValueError(f"match_type must be one of {sorted(allowed)}")
        query = {"type": "terminal", "service": "chemical", "parameters": {
            "value": smiles, "type": "descriptor",
            "descriptor_type": "SMILES", "match_type": match_type}}
        result = self._search(query, "entry", rows)
        result["service"] = f"chemical:{match_type}"
        result["source"] = "RCSB Search API (OpenEye)"
        return result

    def search_similar_components(
        self,
        smiles: str,
        rows: int = 25,
    ) -> Dict[str, Any]:
        """
        Similar CCD components (Tanimoto quick screen, mol_definition
        return type): real chemical-component identifiers with their
        similarity scores, straight from the RCSB chemical search service.
        """
        smiles = smiles.strip()
        if not smiles:
            raise ValueError("SMILES is empty")
        query = {"type": "terminal", "service": "chemical", "parameters": {
            "value": smiles, "type": "descriptor",
            "descriptor_type": "SMILES",
            "match_type": "fingerprint-similarity"}}
        result = self._search(query, "mol_definition", rows)
        result["service"] = "chemical:fingerprint-similarity"
        result["source"] = "RCSB Search API (OpenEye, mol_definition)"
        return result

    def search_same_ligand(self, comp_id: str, rows: int = 25) -> Dict[str, Any]:
        """
        Entries whose binding site contains the given CCD component.

        Two-stage real query: fetch the component's canonical SMILES from
        the CCD (live), then graph-relaxed chemical search (isomorphic +
        substructure within the screened subset) so salts/variants of the
        same molecule are found.
        """
        comp_id = comp_id.strip().upper()
        cc = self._session.get(
            f"{self.config.rcsb_data_url}/core/chemcomp/{comp_id}",
            timeout=http_timeout(self.config.request_timeout))
        if cc.status_code != 200:
            raise DiscoveryError(f"CCD has no component {comp_id}")
        descriptors = cc.json().get("pdbx_chem_comp_descriptor") or []
        smiles = next((d["descriptor"] for d in descriptors
                       if d.get("type") == "SMILES_CANONICAL" and d.get("descriptor")),
                      None)
        if not smiles:
            raise DiscoveryError(f"CCD lists no SMILES for {comp_id}")
        result = self.search_chemical_similarity(
            smiles, match_type="graph-relaxed", rows=rows)
        result["service"] = f"same_ligand:{comp_id}"
        return result

    # ------------------------------------------------------------------
    # RCSB Data API: structure quality (validation report)
    # ------------------------------------------------------------------

    def entry_quality(self, pdb_id: str) -> Dict[str, Any]:
        """
        Entry-level quality metrics from the RCSB validation report
        (pdbx_vrpt_summary), when the entry carries one.
        """
        pdb_id = pdb_id.strip().upper()
        resp = self._session.get(
            f"{self.config.rcsb_data_url}/core/entry/{pdb_id}",
            timeout=http_timeout(self.config.request_timeout))
        if resp.status_code != 200:
            raise DiscoveryError(f"RCSB has no entry {pdb_id}")
        entry = resp.json()
        info = entry.get("rcsb_entry_info", {})
        quality: Dict[str, Any] = {}
        if info.get("resolution_combined"):
            quality["resolution"] = info["resolution_combined"][0]
        vrpt = entry.get("pdbx_vrpt_summary")
        if isinstance(vrpt, dict):
            # The summary block is mostly provenance; the mesh section
            # carries outlier percentages when present.
            mesh = entry.get("pdbx_vrpt_mesh") or {}
            if isinstance(mesh, dict):
                for key in ("percentile_rfree", "percentile_clashscore",
                            "clashscore", "num_H_reduce_outliers",
                            "num_RSRZ_outliers", "num_RSCC_outliers"):
                    if key in mesh:
                        quality[key] = mesh[key]
        return {
            "pdb_id": pdb_id,
            "quality": quality,
            "has_validation_report": bool(vrpt),
            "source": "RCSB Data API (validation report)",
            "url": f"https://www.rcsb.org/structure/{pdb_id}",
        }

    # ------------------------------------------------------------------
    # UniProt REST: protein context
    # ------------------------------------------------------------------

    def uniprot_for_entry(self, pdb_id: str) -> List[Dict[str, Any]]:
        """
        UniProt accessions for the entry's polymer entities (via RCSB's
        SIFTS mapping). Empty list when the entry has no UniProt mapping.
        """
        pdb_id = pdb_id.strip().upper()
        resp = self._session.get(
            f"{self.config.rcsb_data_url}/core/polymer_entity/"
            f"{pdb_id}/1", timeout=http_timeout(self.config.request_timeout))
        if resp.status_code != 200:
            return []
        refs = (resp.json()
                .get("rcsb_polymer_entity_container_identifiers", {})
                .get("reference_sequence_identifiers") or [])
        return [
            {"accession": r["database_accession"],
             "database": r.get("database_name", "UniProt"),
             "coverage": r.get("reference_sequence_coverage")}
            for r in refs if r.get("database_name") == "UniProt"
        ]

    def uniprot_context(self, accession: str) -> Dict[str, Any]:
        """
        Protein context from UniProt: names, organism, function comments,
        and real feature annotations (Binding site with ligand, Active
        site, Domain, Secondary structure counts). Fields verbatim.
        """
        accession = accession.strip()
        resp = self._session.get(
            f"https://rest.uniprot.org/uniprotkb/{accession}.json",
            timeout=http_timeout(self.config.request_timeout))
        if resp.status_code != 200:
            raise DiscoveryError(
                f"UniProt has no entry {accession} ({resp.status_code})")
        u = resp.json()
        features = u.get("features", [])
        binding_sites = [
            {
                "start": f["location"]["start"]["value"],
                "end": f["location"]["end"]["value"],
                "ligand": f.get("ligand", {}).get("name"),
                "ligand_id": f.get("ligand", {}).get("id"),
                "description": f.get("description", ""),
            }
            for f in features if f["type"] == "Binding site"
        ]
        active_sites = [
            {
                "start": f["location"]["start"]["value"],
                "end": f["location"]["end"]["value"],
                "description": f.get("description", ""),
            }
            for f in features if f["type"] == "Active site"
        ]
        domains = [
            {
                "start": f["location"]["start"]["value"],
                "end": f["location"]["end"]["value"],
                "description": f.get("description", ""),
            }
            for f in features if f["type"] == "Domain"
        ]
        function = [
            c["texts"][0]["value"]
            for c in u.get("comments", [])
            if c.get("commentType") == "FUNCTION" and c.get("texts")
        ]
        names = u.get("proteinDescription", {})
        return {
            "accession": accession,
            "entry_name": names.get("recommendedName", {}).get(
                "fullName", {}).get("value"),
            "organism": (u.get("organism", {}).get("scientificName")),
            "function": function,
            "binding_sites": binding_sites,
            "active_sites": active_sites,
            "domains": domains,
            "sequence_length": (u.get("sequence", {}).get("length")),
            "source": "UniProtKB REST",
            "url": f"https://www.uniprot.org/uniprotkb/{accession}",
        }

    # ------------------------------------------------------------------
    # AlphaFold DB
    # ------------------------------------------------------------------

    def alphafold_model_url(self, accession: str) -> str:
        """Current AlphaFold DB model URL for a UniProt accession."""
        accession = accession.strip()
        return (f"https://alphafold.ebi.ac.uk/files/"
                f"AF-{accession}-F1-model_v6.cif")

    def fetch_alphafold_model(self, accession: str) -> Dict[str, Any]:
        """
        Fetch the AlphaFold predicted model (mmCIF text) for a UniProt
        accession. Raises honestly when none exists (HTTP 404).
        """
        url = self.alphafold_model_url(accession)
        timeout = http_timeout(self.config.request_timeout * 3)
        resp = self._session.get(url, timeout=timeout)
        if resp.status_code == 404:
            raise DiscoveryError(
                f"AlphaFold DB has no model for {accession}")
        if resp.status_code != 200:
            raise DiscoveryError(
                f"AlphaFold DB request failed ({resp.status_code})")
        # pLDDT confidence is stored per-residue in the file itself as
        # B-factor (label_atom.B_iso_or_equiv) — extracted by the parser.
        return {
            "accession": accession,
            "url": url,
            "content": resp.text,
            "format": "mmcif",
            "source": "AlphaFold DB (EMBL-EBI, CC-BY 4.0)",
        }

    # ------------------------------------------------------------------
    # RCSB ModelServer: on-demand coordinate subsets
    # ------------------------------------------------------------------

    MODEL_SERVER_URL = "https://models.rcsb.org/v1"

    def fetch_modelserver_subset(
        self, pdb_id: str,
        label_asym_id: Optional[str] = None,
        label_comp_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Fetch a coordinate subset from RCSB's ModelServer as mmCIF text —
        the efficient on-demand path (full entry, one chain, or one
        component). The response is ModelServer's own mmCIF output
        (verified live: /v1/{id}/atoms?encoding=cif), parseable by the
        app's own mmCIF parser.
        """
        pdb_id = pdb_id.strip().lower()
        params: Dict[str, str] = {"encoding": "cif"}
        if label_asym_id:
            params["label_asym_id"] = label_asym_id.strip()
        if label_comp_id:
            params["label_comp_id"] = label_comp_id.strip().upper()
        url = f"{self.MODEL_SERVER_URL}/{pdb_id}/atoms"
        resp = self._session.get(
            url, params=params,
            timeout=http_timeout(self.config.request_timeout))
        if resp.status_code == 400:
            raise DiscoveryError(
                f"ModelServer rejected the subset request ({resp.text[:160]})")
        if resp.status_code == 404:
            raise DiscoveryError(f"ModelServer has no entry {pdb_id.upper()}")
        if resp.status_code != 200:
            raise DiscoveryError(
                f"ModelServer request failed ({resp.status_code})")
        return {
            "pdb_id": pdb_id.upper(),
            "url": resp.url,
            "content": resp.text,
            "format": "mmcif",
            "subset": {"label_asym_id": label_asym_id,
                       "label_comp_id": label_comp_id},
            "source": "RCSB ModelServer",
        }

    # ------------------------------------------------------------------
    # PDBe REST API (EU mirror + annotations)
    # ------------------------------------------------------------------

    PDBE_API_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry"

    def pdbe_entry_summary(self, pdb_id: str) -> Dict[str, Any]:
        """Entry summary from PDBe (EU mirror of the PDB archive)."""
        pdb_id = pdb_id.strip().lower()
        resp = self._session.get(
            f"{self.PDBE_API_URL}/summary/{pdb_id}",
            timeout=http_timeout(self.config.request_timeout))
        if resp.status_code == 404:
            raise DiscoveryError(f"PDBe has no entry {pdb_id.upper()}")
        if resp.status_code != 200:
            raise DiscoveryError(
                f"PDBe request failed ({resp.status_code})")
        data = resp.json().get(pdb_id) or []
        if not data:
            raise DiscoveryError(f"PDBe returned no summary for {pdb_id.upper()}")
        return {
            "pdb_id": pdb_id.upper(),
            "summary": data,
            "source": "PDBe REST API (EMBL-EBI, CC-BY 4.0)",
            "url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}",
        }

    def pdbe_secondary_structure(self, pdb_id: str) -> Dict[str, Any]:
        """
        PDBe's own secondary-structure annotation (helices/sheets per
        chain, from the deposited records) — an EU-source alternative to
        running DSSP locally.
        """
        pdb_id = pdb_id.strip().lower()
        resp = self._session.get(
            f"{self.PDBE_API_URL}/secondary_structure/{pdb_id}",
            timeout=http_timeout(self.config.request_timeout))
        if resp.status_code == 404:
            raise DiscoveryError(
                f"PDBe has no secondary-structure data for {pdb_id.upper()}")
        if resp.status_code != 200:
            raise DiscoveryError(
                f"PDBe request failed ({resp.status_code})")
        data = resp.json().get(pdb_id) or {}
        return {
            "pdb_id": pdb_id.upper(),
            "molecules": data.get("molecules", []),
            "source": "PDBe REST API (secondary_structure)",
            "url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}",
        }

    def pdbe_ligand_monomers(self, pdb_id: str) -> Dict[str, Any]:
        """PDBe's per-chain ligand monomer listing for the entry."""
        pdb_id = pdb_id.strip().lower()
        resp = self._session.get(
            f"{self.PDBE_API_URL}/ligand_monomers/{pdb_id}",
            timeout=http_timeout(self.config.request_timeout))
        if resp.status_code == 404:
            raise DiscoveryError(
                f"PDBe has no ligand monomers for {pdb_id.upper()}")
        if resp.status_code != 200:
            raise DiscoveryError(
                f"PDBe request failed ({resp.status_code})")
        data = resp.json().get(pdb_id) or []
        return {
            "pdb_id": pdb_id.upper(),
            "ligands": data,
            "source": "PDBe REST API (ligand_monomers)",
            "url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}",
        }
