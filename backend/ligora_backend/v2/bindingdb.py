"""
BindingDB REST client — real affinity data (Ki/Kd/IC50/EC50) from the
current BindingDB web services (https://bindingdb.org/rwd/bind/
BindingDBRESTfulAPI.jsp, verified live).

Endpoints used (GET, JSON via &response=application/json):
- /rest/getLigandsByPDBs?pdb={ids}&cutoff={nM}&identity={%}
      → binding measurements for protein–ligand pairs whose protein has
        at least the requested sequence identity to the PDB entry's protein.
- /rest/getLigandsByUniprot?uniprot={acc};{cutoff}
      → binding measurements for a UniProt target.
- /rest/getTargetByCompound?smiles={SMILES}&cutoff={sim}
      → compounds similar to the query SMILES with their targets/affinities.

Every returned field (monomer id, SMILES, affinity type, affinity value in
nM, PMID, DOI) is BindingDB's own output — nothing is recomputed or
invented. When BindingDB has no record for an entry it currently answers
HTTP 500 (a server-side error page); that is surfaced honestly as "no
binding data available", never converted into a fabricated empty success
or a placeholder number.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import requests

from ..config import get_config
from ..net import http_timeout


class BindingDBError(RuntimeError):
    """Raised when the live source cannot answer (no fallbacks exist)."""


class BindingDBClient:
    """Query BindingDB's REST web services for real affinity records."""

    def __init__(self):
        self.config = get_config()
        self._session = requests.Session()
        self._session.headers["User-Agent"] = "Ligora/0.1 (open-source)"

    @property
    def _base(self) -> str:
        # Config default is https://bindingdb.org; endpoint paths below are
        # the service's own documented REST paths.
        return (self.config.bindingdb_base_url or
                "https://bindingdb.org").rstrip("/")

    def _get_json(self, path: str, params: Dict[str, str]) -> Dict[str, Any]:
        url = f"{self._base}{path}"
        try:
            resp = self._session.get(
                url, params=params,
                timeout=http_timeout(self.config.request_timeout))
        except requests.RequestException as e:
            raise BindingDBError(f"BindingDB unreachable: {e}") from e

        if resp.status_code != 200:
            # BindingDB's own contract quirk: an entry with no record can
            # come back as HTTP 500 (its service errors server-side).
            # Report it honestly with the distinction the service itself
            # makes, and never invent an empty-but-successful answer.
            detail = resp.text[:160].replace("\n", " ")
            if resp.status_code == 500 and "polymerid IN (  )" in resp.text:
                raise BindingDBError(
                    "BindingDB has no record for this entry "
                    "(service returned an error for an unknown PDB ID)")
            raise BindingDBError(
                f"BindingDB request failed ({resp.status_code}): {detail}")

        try:
            return resp.json()
        except ValueError as e:
            # A 200 with a zero-byte body is a real service behavior for
            # some accessions (observed live on getLigandsByUniprot for
            # P00004). It is not an empty-but-valid answer.
            if not resp.text.strip():
                raise BindingDBError(
                    "BindingDB returned an empty response body "
                    "for this query") from e
            raise BindingDBError(
                f"BindingDB returned non-JSON content: "
                f"{resp.text[:120]!r}") from e

    # ------------------------------------------------------------------
    # Response normalization (BindingDB's own field names, verbatim)
    # ------------------------------------------------------------------

    @staticmethod
    def _response_payload(data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Locate the response envelope. BindingDB's JSON wraps records in a
        service-specific key ('getLindsByPDBsResponse',
        'getLindsByUniprotResponse', ...) whose exact spelling is the
        service's own — so the envelope is discovered, not guessed.
        """
        for value in data.values():
            if isinstance(value, dict) and (
                    "affinities" in value or "bdb.affinities" in value):
                return value
        return None

    @classmethod
    def _normalize_affinities(
            cls, data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Extract affinity records. Two real response shapes exist:
        - getLigandsByPDBs: bare keys (monomerid, smile, affinity_type,
          affinity, pmid, doi)
        - getLigandsByUniprot / getTargetByCompound: 'bdb.'-prefixed keys
          (bdb.monomerid, bdb.smile, bdb.affinity_type, bdb.affinity)
        """
        payload = cls._response_payload(data)
        if payload is None:
            raw = []
        elif "affinities" in payload:
            raw = payload.get("affinities") or []
        else:
            raw = payload.get("bdb.affinities") or []

        records: List[Dict[str, Any]] = []
        for item in raw:
            if not isinstance(item, dict):
                continue
            records.append({
                "monomer_id": item.get("monomerid",
                                       item.get("bdb.monomerid")),
                "smiles": item.get("smile", item.get("bdb.smile")),
                "affinity_type": item.get("affinity_type",
                                          item.get("bdb.affinity_type")),
                "affinity_nm": item.get("affinity",
                                        item.get("bdb.affinity")),
                "pmid": item.get("pmid", item.get("bdb.pmid")),
                "doi": item.get("doi", item.get("bdb.doi")),
                "query": item.get("query"),
            })
        return records

    # ------------------------------------------------------------------
    # Public queries
    # ------------------------------------------------------------------

    def ligands_by_pdb(self, pdb_id: str, affinity_cutoff_nm: int = 10000,
                       identity_cutoff: Optional[int] = None) -> Dict[str, Any]:
        """
        Binding measurements for protein–ligand pairs tied to a PDB ID.

        BindingDB matches the entry's protein against its polymer index at
        a sequence-identity threshold. At a strict threshold some real
        entries fall outside the index (observed live: 1HRC errors at 90%
        but returns genuine cytochrome-c binders at 40%), so when the
        service reports no polymer match the query is re-issued at the
        service's own looser documented thresholds — same endpoint, same
        source, and the identity actually used is reported.
        """
        pdb = (pdb_id or "").strip().upper()
        if len(pdb) != 4 or not pdb.isalnum():
            raise BindingDBError(f"Invalid PDB ID: {pdb_id!r}")

        thresholds = ([int(identity_cutoff)] if identity_cutoff is not None
                      else [90, 70, 40])
        last_error: Optional[BindingDBError] = None
        for identity in thresholds:
            try:
                return self._ligands_by_pdb_at(
                    pdb, int(affinity_cutoff_nm), identity)
            except BindingDBError as e:
                last_error = e
                if "no record for this entry" not in str(e):
                    raise
        raise last_error  # type: ignore[misc]

    def _ligands_by_pdb_at(self, pdb: str, affinity_cutoff_nm: int,
                           identity: int) -> Dict[str, Any]:
        data = self._get_json("/rest/getLigandsByPDBs", {
            "pdb": pdb,
            "cutoff": str(int(affinity_cutoff_nm)),
            "identity": str(int(identity)),
            "response": "application/json",
        })
        records = self._normalize_affinities(data)
        return {
            "source": "BindingDB REST (getLigandsByPDBs)",
            "url": f"{self._base}/rest/getLigandsByPDBs",
            "pdb_id": pdb,
            "affinity_cutoff_nm": int(affinity_cutoff_nm),
            "identity_cutoff": int(identity),
            "record_count": len(records),
            "records": records,
        }

    def ligands_by_uniprot(self, accession: str,
                           affinity_cutoff_nm: int = 10000
                           ) -> Dict[str, Any]:
        """Binding measurements for a UniProt target."""
        acc = (accession or "").strip()
        if not acc:
            raise BindingDBError("A UniProt accession is required")
        data = self._get_json("/rest/getLigandsByUniprot", {
            "uniprot": f"{acc};{int(affinity_cutoff_nm)}",
            "response": "application/json",
        })
        records = self._normalize_affinities(data)
        payload = self._response_payload(data) or {}
        return {
            "source": "BindingDB REST (getLigandsByUniprot)",
            "url": f"{self._base}/rest/getLigandsByUniprot",
            "uniprot": acc,
            # BindingDB's own service-side hit count and accession set.
            "hit_count": payload.get("bdb.hit"),
            "primary_accession": payload.get("bdb.primary"),
            "affinity_cutoff_nm": int(affinity_cutoff_nm),
            "record_count": len(records),
            "records": records,
        }

    def targets_by_compound(self, smiles: str,
                            similarity_cutoff: float = 0.85
                            ) -> Dict[str, Any]:
        """Targets/affinities for compounds similar to the query SMILES."""
        smi = (smiles or "").strip()
        if not smi:
            raise BindingDBError(
                "A SMILES string is required (the ligand has none resolved)")
        data = self._get_json("/rest/getTargetByCompound", {
            "smiles": smi,
            "cutoff": str(float(similarity_cutoff)),
            "response": "application/json",
        })
        records = self._normalize_affinities(data)
        return {
            "source": "BindingDB REST (getTargetByCompound)",
            "url": f"{self._base}/rest/getTargetByCompound",
            "smiles": smi,
            "similarity_cutoff": float(similarity_cutoff),
            "record_count": len(records),
            "records": records,
        }

    # ------------------------------------------------------------------
    # Convenience: entry-level lookup with UniProt fallback
    # ------------------------------------------------------------------

    def affinity_for_entry(
        self, pdb_id: str,
        uniprot_accession: Optional[str] = None,
        affinity_cutoff_nm: int = 10000,
    ) -> Dict[str, Any]:
        """
        Binding data for a PDB entry. Queries by PDB ID first; when the
        entry is not in BindingDB's PDB coverage and a UniProt accession
        is supplied, queries by UniProt instead. When neither answers,
        the BindingDBError from the last attempt is raised — honestly.
        """
        errors: List[str] = []
        try:
            return self.ligands_by_pdb(
                pdb_id, affinity_cutoff_nm=affinity_cutoff_nm)
        except BindingDBError as e:
            errors.append(f"by PDB ID: {e}")
        if uniprot_accession:
            try:
                return self.ligands_by_uniprot(
                    uniprot_accession, affinity_cutoff_nm=affinity_cutoff_nm)
            except BindingDBError as e:
                errors.append(f"by UniProt: {e}")
        raise BindingDBError(
            "BindingDB has no binding data for this entry ("
            + "; ".join(errors) + ")")
