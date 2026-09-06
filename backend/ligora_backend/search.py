"""
BM25 search over live-fetched documents.

The index is a cache of real documents: RCSB entry metadata, CCD/PubChem
chemical records, and the user's own workspace artifacts. Indexing fetches
from the providers listed in data/DATA_SOURCES.md; query-time scoring is
pure Okapi BM25 math over those documents (k1=1.5, b=0.75) — no network at
query time and no synthesized hits. An empty index yields zero hits, which
the UI presents as "no indexed documents yet".
"""

from __future__ import annotations

import json
import math
import re
import time
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from .config import get_config

_TOKEN_RE = re.compile(r"[a-z0-9]+")

# Okapi BM25 parameters (standard, not tunable per query).
_K1 = 1.5
_B = 0.75

# Field boosts: identifiers match stronger than free text.
_FIELD_WEIGHTS = {
    "id": 3.0,
    "residue_name": 3.0,
    "title": 1.6,
    "name": 1.6,
    "journal_abbrev": 1.1,
    "formula": 1.3,
    "smiles": 1.0,
    "classification": 1.2,
    "text": 1.0,
}


def _tokenize(text: str) -> List[str]:
    return _TOKEN_RE.findall(str(text).lower())


class _Doc:
    __slots__ = ("doc_id", "scope", "fields", "fetched_at", "source", "url",
                 "tfs", "length")

    def __init__(self, doc_id: str, scope: str, fields: Dict[str, str],
                 fetched_at: str, source: str, url: str = ""):
        self.doc_id = doc_id
        self.scope = scope
        self.fields = fields
        self.fetched_at = fetched_at
        self.source = source
        self.url = url
        weighted: List[str] = []
        for field, value in fields.items():
            weight = _FIELD_WEIGHTS.get(field, 1.0)
            tokens = _tokenize(value)
            weighted.extend(tokens * max(1, int(round(weight))))
        self.tfs = Counter(weighted)
        self.length = len(weighted)


class SearchIndex:
    """Okapi BM25 index over documents fetched from live sources."""

    def __init__(self, cache_dir: Optional[Path] = None):
        self.config = get_config()
        self._session = requests.Session()
        self._session.headers["User-Agent"] = (
            "Ligora/0.1 (open-source molecular analysis workstation)")
        self._cache_dir = Path(
            cache_dir or Path(self.config.base_dir) / "search_index")
        self._docs: Dict[str, _Doc] = {}
        self._df: Counter = Counter()
        self._load()

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _index_file(self) -> Path:
        return self._cache_dir / "index.json"

    def _load(self):
        try:
            raw = json.loads(self._index_file().read_text())
        except (OSError, json.JSONDecodeError):
            return
        for entry in raw.get("documents", []):
            doc = _Doc(
                entry["doc_id"], entry["scope"], entry["fields"],
                entry["fetched_at"], entry["source"], entry.get("url", ""))
            self._docs[doc.doc_id] = doc
        self._recompute_df()

    def _save(self):
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "documents": [
                {
                    "doc_id": d.doc_id, "scope": d.scope, "fields": d.fields,
                    "fetched_at": d.fetched_at, "source": d.source,
                    "url": d.url,
                }
                for d in self._docs.values()
            ]
        }
        self._index_file().write_text(json.dumps(payload, indent=1))

    def _recompute_df(self):
        self._df = Counter()
        for doc in self._docs.values():
            for term in doc.tfs:
                self._df[term] += 1

    # ------------------------------------------------------------------
    # Fetching (live sources only)
    # ------------------------------------------------------------------

    @staticmethod
    def _now() -> str:
        return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    def add_structure(self, pdb_id: str) -> Optional[_Doc]:
        """Fetch real entry metadata from the RCSB Data API and index it."""
        pdb_id = pdb_id.strip().upper()
        # rcsb_data_url already ends with /rest/v1.
        url = f"{self.config.rcsb_data_url}/core/entry/{pdb_id}"
        try:
            resp = self._session.get(url, timeout=self.config.request_timeout)
            if resp.status_code != 200:
                return None
            entry = resp.json()
        except requests.RequestException:
            return None
        citation = entry.get("rcsb_primary_citation") or {}
        info = entry.get("rcsb_entry_info") or {}
        fields = {
            "id": pdb_id,
            "title": entry.get("struct", {}).get("title", ""),
            "journal_abbrev": citation.get("journal_abbrev", "") or "",
            "text": citation.get("title", "") or "",
            "classification": info.get("structure_determination_methodology", "")
            or entry.get("exptl", {}).get("method", "") or "",
        }
        doc = _Doc(f"struct:{pdb_id}", "structures", fields, self._now(),
                   "RCSB PDB Data API", url)
        self._docs[doc.doc_id] = doc
        self._recompute_df()
        self._save()
        return doc

    def add_chemical(self, comp_id: str) -> Optional[_Doc]:
        """Fetch a real CCD component (classification, formula, synonyms)."""
        comp_id = comp_id.strip().upper()
        # rcsb_data_url already ends with /rest/v1.
        url = f"{self.config.rcsb_data_url}/core/chemcomp/{comp_id}"
        try:
            resp = self._session.get(url, timeout=self.config.request_timeout)
            if resp.status_code != 200:
                return None
            payload = resp.json()
            comp = payload.get("chem_comp", {})
        except requests.RequestException:
            return None
        if not comp:
            return None
        synonyms = " ".join(
            (s.get("name") or "")
            for s in payload.get("rcsb_chem_comp_synonyms", []))
        fields = {
            "id": comp_id,
            "residue_name": comp_id,
            "name": comp.get("name", "") or "",
            "formula": comp.get("formula", "") or "",
            "classification": comp.get("type", "") or "",
            "text": synonyms or (comp.get("name", "") or ""),
        }
        doc = _Doc(f"chem:{comp_id}", "chem", fields, self._now(),
                   "RCSB Chemical Component Dictionary", url)
        self._docs[doc.doc_id] = doc
        self._recompute_df()
        self._save()
        return doc

    def add_pubchem_compound(self, cid: int) -> Optional[_Doc]:
        """Fetch real PubChem properties for a CID and index them."""
        props = "Title,MolecularFormula,ConnectivitySMILES,IUPACName"
        url = (f"{self.config.pubchem_base_url}/compound/cid/{cid}/"
               f"property/{props}/JSON")
        try:
            resp = self._session.get(url, timeout=self.config.request_timeout)
            if resp.status_code != 200:
                return None
            prop = (resp.json().get("PropertyTable", {})
                    .get("Properties") or [{}])[0]
        except requests.RequestException:
            return None
        if not prop:
            return None
        fields = {
            "id": f"CID{prop.get('CID', cid)}",
            "name": prop.get("Title", "") or "",
            "formula": prop.get("MolecularFormula", "") or "",
            "smiles": prop.get("ConnectivitySMILES", "") or "",
            "text": prop.get("IUPACName", "") or "",
        }
        doc = _Doc(f"pubchem:{prop.get('CID', cid)}", "chem", fields,
                   self._now(), "PubChem PUG REST", url)
        self._docs[doc.doc_id] = doc
        self._recompute_df()
        self._save()
        return doc

    def add_workspace_artifact(self, path: Path) -> Optional[_Doc]:
        """Index the user's own saved analysis artifact (real local data)."""
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            return None
        fields = {
            "id": str(data.get("structure_id", path.stem)),
            "title": str(data.get("structure_title", "")),
            "name": str(data.get("ligand_name", "")),
            "text": " ".join(
                str(c.get("contact_type", ""))
                for c in data.get("contacts", [])[:50]),
        }
        doc = _Doc(f"artifact:{path.name}", "workspace", fields,
                   self._now(), "local workspace", str(path))
        self._docs[doc.doc_id] = doc
        self._recompute_df()
        self._save()
        return doc

    def refresh(self, pdb_ids: Optional[List[str]] = None,
                chem_names: Optional[List[str]] = None) -> Dict[str, int]:
        """Re-fetch the given entities from live sources (and index them)."""
        added = 0
        for pdb_id in pdb_ids or []:
            if self.add_structure(pdb_id):
                added += 1
        for name in chem_names or []:
            if self.add_chemical(name):
                added += 1
        # Workspace artifacts under the active base dir.
        for summary in Path(self.config.base_dir).glob(
                "workspaces/*/analysis_summary.json"):
            if self.add_workspace_artifact(summary):
                added += 1
        return {"indexed": added, "total": len(self._docs)}

    # ------------------------------------------------------------------
    # BM25 scoring
    # ------------------------------------------------------------------

    def _bm25(self, query_terms: List[str], scope: Optional[str],
              limit: int) -> List[Dict[str, Any]]:
        n_docs = len(self._docs)
        if n_docs == 0:
            return []
        avgdl = sum(d.length for d in self._docs.values()) / n_docs
        scores: List[tuple] = []
        for doc in self._docs.values():
            if scope and doc.scope != scope:
                continue
            score = 0.0
            for term in query_terms:
                tf = doc.tfs.get(term, 0)
                if not tf:
                    continue
                df = self._df.get(term, 0)
                idf = math.log(1.0 + (n_docs - df + 0.5) / (df + 0.5))
                denom = tf + _K1 * (1 - _B + _B * doc.length / (avgdl or 1.0))
                score += idf * (tf * (_K1 + 1)) / denom
            if score > 0:
                scores.append((score, doc))
        scores.sort(key=lambda pair: -pair[0])
        hits = []
        for score, doc in scores[:limit]:
            hits.append({
                "doc_id": doc.doc_id,
                "scope": doc.scope,
                "score": round(score, 4),
                "fields": doc.fields,
                "source": doc.source,
                "url": doc.url,
                "fetched_at": doc.fetched_at,
            })
        return hits

    def search(self, query: str, scope: Optional[str] = None,
               limit: int = 20) -> Dict[str, Any]:
        """Rank indexed documents for the query with Okapi BM25."""
        terms = _tokenize(query)
        if not terms:
            return {"query": query, "total": 0, "hits": []}
        hits = self._bm25(terms, scope, limit)
        return {
            "query": query,
            "total": len(hits),
            "documents_indexed": len(self._docs),
            "hits": hits,
        }
