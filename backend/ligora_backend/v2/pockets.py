"""
Pocket detection via the real fpocket binary.

fpocket (Le Guilloux et al. 2009, J. Chem. Inf. Model; maintained by
Discngine) is invoked as an external engine — never reimplemented. All
pocket geometry, descriptors and scores come from fpocket's own output
files; nothing is approximated here.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..config import get_config
from ..schemas import Ligand, Structure


def _write_protein_pdb(structure: Structure, path: Path) -> int:
    """
    Write the structure's polymer atoms as a plain PDB for fpocket.

    Same strict column layout the PLIP writer uses (ATOM records only,
    element right-justified in 77-78). Returns the atom count.
    """
    lines: List[str] = []
    serial = 1
    for chain in structure.chains:
        if not chain.is_polymer:
            continue
        for residue in chain.residues:
            for atom in residue.atoms:
                element = (atom.element or "").upper()[:2].rjust(2)
                lines.append(
                    f"{'ATOM':<6}{serial:>5}  "
                    f"{atom.name:<4}"
                    f"{residue.name:<3} "
                    f"{chain.id:<1}"
                    f"{residue.id:>4}    "
                    f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
                    f"{atom.occupancy:>6.2f}{atom.b_factor:>6.2f}"
                    f"{'':>10}"
                    f"{element}"
                )
                serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return serial - 1


_DESC_KEYS = [
    "Score", "Druggability Score", "Number of Alpha Spheres",
    "Total SASA", "Polar SASA", "Apolar SASA", "Volume",
    "Mean local hydrophobic density", "Mean alpha sphere radius",
    "Mean alp. sph. solvent access", "Apolar alpha sphere proportion",
    "Hydrophobicity score", "Volume score", "Polarity score",
    "Charge score", "Proportion of polar atoms", "Alpha sphere density",
    "Cent. of mass - Alpha Sphere max dist", "Flexibility",
]


def _parse_fpocket_info(info_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Parse fpocket's <stem>_info.txt into {pocket_id: {descriptor: value}}.

    The file is fpocket's own machine-generated report; values are kept
    verbatim (numeric when numeric) — no derived quantities added.
    """
    result: Dict[str, Dict[str, Any]] = {}
    current: Optional[str] = None
    for line in info_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"^Pocket\s+(\d+)\s*:", line.strip())
        if m:
            current = m.group(1)
            result[current] = {}
            continue
        if current is None or ":" not in line:
            continue
        key, _, value = line.partition(":")
        key = key.strip().lstrip("\t").strip()
        value = value.strip()
        if not key:
            continue
        try:
            result[current][key] = float(value)
        except ValueError:
            result[current][key] = value
    return result


class PocketDetector:
    """Detect binding pockets with the real fpocket engine."""

    def __init__(self, binary_path: Optional[str] = None):
        self.config = get_config()
        self._binary_override = binary_path or os.environ.get(
            "LIGORA_FPOCKET_PATH")

    def _binary(self) -> Optional[str]:
        for candidate in (self._binary_override,):
            if candidate and Path(candidate).exists():
                return candidate
        from shutil import which
        found = which("fpocket")
        if found:
            return found
        snap = Path("/snap/bin/fpocket")
        if snap.exists():
            return str(snap)
        return None

    def is_available(self) -> bool:
        """True when the real fpocket binary is present and runnable."""
        binary = self._binary()
        if not binary:
            return False
        try:
            subprocess.run(
                [binary], capture_output=True, timeout=15, text=True,
                input="")
            return True  # the binary launched; that is the availability check
        except (OSError, subprocess.TimeoutExpired):
            return False

    def detect_pockets(
        self,
        structure: Structure,
        workdir: Path,
        ligand: Optional[Ligand] = None,
    ) -> Dict[str, Any]:
        """
        Run fpocket on the structure's polymer atoms and parse its output.

        Returns a dict with:
          - available: False when the engine is absent (honest state)
          - pockets: fpocket's own pocket records (descriptors verbatim)
          - ligand_overlap: when a ligand is given, per-pocket minimum
            distance and atom-containment counts from real coordinates
        """
        binary = self._binary()
        if not binary:
            return {
                "available": False,
                "error": "fpocket binary not found "
                         "(snap install fpocket or LIGORA_FPOCKET_PATH)",
                "pockets": [],
            }

        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)

        # A snap-confined fpocket can only see non-hidden paths under
        # $HOME (the snap 'home' interface excludes dot-directories), so
        # the run is staged there and results are copied back into the
        # session workspace for the record.
        stage_root = Path(os.environ.get(
            "LIGORA_FPOCKET_STAGE", str(Path.home() / "ligora-fpocket")))
        stage = stage_root / f"job_{os.getpid()}"
        stage.mkdir(parents=True, exist_ok=True)
        pdb_path = stage / "fpocket_input.pdb"
        n_atoms = _write_protein_pdb(structure, pdb_path)
        if n_atoms == 0:
            return {
                "available": True,
                "error": "Structure has no polymer atoms to analyze",
                "pockets": [],
            }

        try:
            proc = subprocess.run(
                [binary, "-f", str(pdb_path)],
                capture_output=True, timeout=600, text=True,
                cwd=str(stage))
        except (OSError, subprocess.TimeoutExpired) as e:
            return {
                "available": True,
                "error": f"fpocket failed to run: {e}",
                "pockets": [],
            }

        if proc.returncode != 0:
            return {
                "available": True,
                "error": (f"fpocket exited {proc.returncode}: "
                          f"{(proc.stderr or proc.stdout)[-300:]}"),
                "pockets": [],
            }

        # fpocket 4.x writes <stem>_out/ next to the input:
        #   <stem>_out/<stem>_info.txt        per-pocket descriptors
        #   <stem>_out/pockets/pocketN_atm.pdb  pocket atom sets
        out_dir = stage / (pdb_path.stem + "_out")
        pockets_dir = out_dir / "pockets"
        info_path = out_dir / (pdb_path.stem + "_info.txt")
        if not pockets_dir.is_dir():
            return {
                "available": True,
                "error": "fpocket produced no pockets directory",
                "pockets": [],
            }
        descriptors = _parse_fpocket_info(info_path) if info_path.exists() else {}

        # Preserve fpocket's own outputs in the session workspace.
        record_dir = workdir / "fpocket_out"
        try:
            import shutil
            if record_dir.exists():
                shutil.rmtree(record_dir)
            shutil.copytree(out_dir, record_dir)
        except OSError:
            pass

        pockets: List[Dict[str, Any]] = []
        for pf in sorted(pockets_dir.glob("pocket*_atm.pdb")):
            m = re.match(r"pocket(\d+)_atm", pf.stem)
            pid = m.group(1) if m else pf.stem
            coords = []
            for line in pf.read_text().splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        coords.append([float(line[30:38]),
                                       float(line[38:46]),
                                       float(line[46:54])])
                    except ValueError:
                        continue
            if not coords:
                continue
            pockets.append({
                "pocket_id": pid,
                "atom_count": len(coords),
                "center": [round(sum(c[i] for c in coords) / len(coords), 3)
                           for i in range(3)],
                "coordinates": coords,
                "descriptors": descriptors.get(pid, {}),
            })

        response: Dict[str, Any] = {
            "available": True,
            "binary": binary,
            "pocket_count": len(pockets),
            "pockets": pockets,
        }
        if ligand is not None and ligand.atoms:
            response["ligand_overlap"] = self.pocket_ligand_overlap(
                pockets, ligand)
        return response

    def pocket_ligand_overlap(
        self,
        pockets: List[Dict[str, Any]],
        ligand: Ligand,
    ) -> List[Dict[str, Any]]:
        """
        Per-pocket relationship to the co-crystallized ligand, from real
        coordinates only: minimum atom-atom distance and how many ligand
        atoms sit within 4 Å of any pocket atom. No invented scores.
        """
        import numpy as np
        lig_coords = np.array(
            [[a.x, a.y, a.z] for a in ligand.atoms], dtype=float)
        if lig_coords.ndim != 2 or lig_coords.shape[0] == 0:
            return []
        results = []
        for p in pockets:
            pc = np.array(p["coordinates"], dtype=float)
            if pc.ndim != 2 or pc.shape[0] == 0:
                continue
            d = np.sqrt(((pc[:, None, :] - lig_coords[None, :, :]) ** 2)
                        .sum(-1))
            results.append({
                "pocket_id": p["pocket_id"],
                "min_distance": round(float(d.min()), 2),
                "ligand_atoms_within_4A": int((d.min(axis=0) < 4.0).sum()),
                "ligand_atom_count": int(lig_coords.shape[0]),
                "score": p["descriptors"].get("Score"),
                "druggability": p["descriptors"].get("Druggability Score"),
            })
        results.sort(key=lambda r: r["min_distance"])
        return results
