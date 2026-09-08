"""
Covalent-ligand handling, grounded entirely in the structure file's own
records — never in distance-based guessing:

- mmCIF: the `_struct_conn` category as deposited, rows whose
  `conn_type_id` is `covale` or `disulf` (wwPDB mmCIF dictionary). The
  record carries both partners (label_asym/comp/seq/atom ids) and the
  deposited `pdbx_dist_value`.
- PDB: `LINK` records (fixed-column legacy format). Partner distance is
  computed from the structure's real atom coordinates (same math the
  measurement tool uses) because the legacy record's own distance column
  is present only in newer files.

When a selected ligand is covalently attached to the polymer, Ligora
reports it on the ligand card and in the analysis: covalent attachment is
depositor-stated chemistry, not an interaction type, so PLIP's
non-covalent contact table stays untouched — it is annotated instead.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..config import get_config
from ..net import http_timeout
from ..schemas import Structure

# conn_type_id values (mmCIF _struct_conn dictionary) that state a covalent
# connection. Anything else (hydrog/saltbr/...) is PLIP's/deppositor's
# non-covalent annotation space and is not handled here.
_COVALENT_CONN_TYPES = {"covale", "disulf"}


def _resolve_file(structure: Structure) -> Optional[Path]:
    """Locate the structure's own file on disk (never fabricated)."""
    if structure.file_path:
        p = Path(structure.file_path)
        if p.exists():
            return p
    # Sessions opened by PDB ID cache the fetched file in the workspace;
    # when that is gone but the entry is a real PDB ID, fetch it once.
    pdb_id = (structure.id or "").strip().upper()
    if structure.source in ("rcsb", "pdb", "csm") and len(pdb_id) == 4:
        import requests
        config = get_config()
        # Fetch the same format the structure was parsed from, so the
        # records read match the record set the analysis was built on.
        extensions = ((".pdb", ".cif") if structure.file_format == "pdb"
                      else (".cif", ".pdb"))
        for ext in extensions:
            url = f"{config.rcsb_files_url}/download/{pdb_id}{ext}"
            try:
                resp = requests.get(
                    url, timeout=http_timeout(config.request_timeout))
            except requests.RequestException:
                continue
            if resp.status_code == 200 and resp.text.strip():
                import tempfile
                tmp = Path(tempfile.gettempdir()) / f"ligora_{pdb_id}{ext}"
                tmp.write_text(resp.text, encoding="utf-8")
                return tmp
    return None


def _read_struct_conn(text: str) -> List[Dict[str, Any]]:
    """Parse `_struct_conn` rows from mmCIF using the app's own CIF reader."""
    from ..parser import StructureParser
    parser = StructureParser()
    data = parser._parse_mmcif_data(text)
    rows = parser._get_loop_rows(data, "_struct_conn")
    links: List[Dict[str, Any]] = []
    for row in rows:
        conn_type = (row.get("conn_type_id") or "").strip().lower()
        if conn_type not in _COVALENT_CONN_TYPES:
            continue

        links.append({
            "conn_type": conn_type,
            "partner1": {
                "chain": row.get("ptnr1_label_asym_id"),
                "comp": row.get("ptnr1_label_comp_id"),
                "seq": row.get("ptnr1_label_seq_id"),
                "atom": row.get("ptnr1_label_atom_id"),
            },
            "partner2": {
                "chain": row.get("ptnr2_label_asym_id"),
                "comp": row.get("ptnr2_label_comp_id"),
                "seq": row.get("ptnr2_label_seq_id"),
                "atom": row.get("ptnr2_label_atom_id"),
            },
            "distance": _as_float(row.get("pdbx_dist_value")),
            "record": "_struct_conn",
        })
    return links


def _as_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


_LINK_RE = re.compile(r"^LINK[R]? ")

# PDB LINK fixed columns (wwPDB v3.3 format):
#   partner 1: atom 13-16, altLoc 17, resName 18-20, chain 22, resSeq 23-26
#   partner 2: atom 43-46, altLoc 47, resName 48-50, chain 52, resSeq 53-56


def _pdb_field(line: str, start: int, end: int) -> str:
    return line[start:end].strip() if len(line) >= start else ""


def _read_pdb_links(text: str) -> List[Dict[str, Any]]:
    links: List[Dict[str, Any]] = []
    for line in text.splitlines():
        if not _LINK_RE.match(line):
            continue
        links.append({
            "conn_type": "link",
            "partner1": {
                "atom": _pdb_field(line, 12, 16),
                "comp": _pdb_field(line, 17, 20),
                "chain": _pdb_field(line, 21, 22),
                "seq": _pdb_field(line, 22, 26),
            },
            "partner2": {
                "atom": _pdb_field(line, 42, 46),
                "comp": _pdb_field(line, 47, 50),
                "chain": _pdb_field(line, 51, 52),
                "seq": _pdb_field(line, 52, 56),
            },
            "distance": None,
            "record": line[:6].strip(),
        })
    return links


def _atom_lookup(structure: Structure) -> tuple:
    """
    Atom indexes for distance computation:
    - primary: (chain, resseq, atomname) -> Atom (polymer + any chain atom)
    - fallback: (comp, atomname) -> Atom from the ligand objects (non-
      polymer components live in `structure.ligands`, not in chains).
    """
    lookup: Dict[tuple, Any] = {}
    by_comp: Dict[tuple, Any] = {}
    for chain in structure.chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                lookup[(chain.id, residue.id, atom.name)] = atom
                by_comp.setdefault(
                    (residue.name.strip().upper(), atom.name), atom)
    for ligand in structure.ligands:
        for atom in ligand.atoms:
            by_comp.setdefault(
                (ligand.residue_name.strip().upper(), atom.name), atom)
    return lookup, by_comp


def _distance(a: Any, b: Any) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def detect_covalent_links(structure: Structure) -> Dict[str, Any]:
    """
    Report every covalent connection the structure file itself declares
    between a non-polymer ligand component and any other component.

    Returns {"available", "links", "covalent_ligands", "source"}.
    Distances come from the deposited record (mmCIF pdbx_dist_value) or
    from real atom coordinates (PDB LINK); never from a threshold guess.
    """
    path = _resolve_file(structure)
    if not path:
        return {
            "available": False,
            "error": ("The structure's source file is not available, so "
                      "its covalent-link records cannot be read. Nothing "
                      "is inferred from geometry."),
            "links": [],
            "covalent_ligands": [],
        }

    text = path.read_text(encoding="utf-8", errors="replace")
    # Format: trust the parser's own record first, then the file's
    # extension; only sniff the content when both are silent (the PDB
    # format's REMARKs can legitimately mention mmCIF category names).
    if structure.file_format == "mmcif":
        is_cif = True
    elif structure.file_format == "pdb":
        is_cif = False
    else:
        is_cif = path.suffix in (".cif", ".cif.gz") or \
            "_struct_conn." in text[:200000]
    if is_cif:
        raw_links = _read_struct_conn(text)
    else:
        raw_links = _read_pdb_links(text)

    ligand_comps = {
        (lig.residue_name or "").strip().upper()
        for lig in structure.ligands
    }

    lookup, by_comp = _atom_lookup(structure)
    links: List[Dict[str, Any]] = []
    covalent_ligands: Dict[str, List[Dict[str, Any]]] = {}

    for link in raw_links:
        p1, p2 = link["partner1"], link["partner2"]
        comp1 = (p1.get("comp") or "").strip().upper()
        comp2 = (p2.get("comp") or "").strip().upper()
        if comp1 not in ligand_comps and comp2 not in ligand_comps:
            continue

        # Distance: deposited value when the record carries one; otherwise
        # real coordinates of the two partner atoms.
        distance = link.get("distance")
        if distance is None:
            def _find(p: Dict[str, Any]) -> Any:
                atom = p.get("atom") or ""
                try:
                    seq = int(str(p.get("seq") or "").strip())
                except (TypeError, ValueError):
                    seq = None
                hit = None
                if seq is not None:
                    hit = lookup.get((p.get("chain") or "", seq, atom))
                if hit is None:
                    hit = by_comp.get(((p.get("comp") or "").strip()
                                       .upper(), atom))
                return hit
            a, b = _find(p1), _find(p2)
            if a is not None and b is not None:
                distance = round(_distance(a, b), 3)

        record = {
            "conn_type": link["conn_type"],
            "record": link["record"],
            "partner1": p1,
            "partner2": p2,
            "distance": distance,
        }
        links.append(record)
        ligand_comp = comp1 if comp1 in ligand_comps else comp2
        covalent_ligands.setdefault(ligand_comp, []).append(record)

    return {
        "available": True,
        "source": ("mmCIF _struct_conn (as deposited)"
                   if is_cif else "PDB LINK records (as deposited)"),
        "file": str(path),
        "links": links,
        "covalent_ligands": covalent_ligands,
    }


def covalent_flag_for_ligand(
    structure: Structure, residue_name: str,
    detection: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """
    The covalent attachments declared for one ligand component, for the
    ligand card. Detection runs once and can be passed in.
    """
    det = detection or detect_covalent_links(structure)
    if not det.get("available"):
        return []
    return det["covalent_ligands"].get(
        (residue_name or "").strip().upper(), [])
