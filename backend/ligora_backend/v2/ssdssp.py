"""
Secondary structure assignment via the real DSSP engine.

DSSP (Kabsch & Sander 1983; current implementation by NKI, 4.x) is invoked
as an external binary — never reimplemented. Assignments, accessible
surface, and hydrogen-bond statistics come from DSSP's own output.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..config import get_config
from ..schemas import Structure
from .pockets import _write_protein_pdb


class DSSPAnalyzer:
    """Run mkdssp on the structure's polymer chains and parse its output."""

    def __init__(self, binary_path: Optional[str] = None):
        self.config = get_config()
        self._binary_override = binary_path or os.environ.get(
            "LIGORA_DSSP_PATH")

    def _binary(self) -> Optional[str]:
        # An explicitly configured binary is honored as-is: when it does
        # not exist we must report that honestly, never silently swap in
        # whatever else happens to be on PATH.
        if self._binary_override:
            return self._binary_override \
                if Path(self._binary_override).exists() else None
        from shutil import which
        for name in ("mkdssp", "dssp"):
            found = which(name)
            if found:
                return found
        return None

    def is_available(self) -> bool:
        binary = self._binary()
        if not binary:
            return False
        try:
            proc = subprocess.run(
                [binary, "--version"], capture_output=True, timeout=15,
                text=True)
            return proc.returncode == 0 or bool(proc.stdout)
        except (OSError, subprocess.TimeoutExpired):
            return False

    # DSSP 8-state codes and their labels (from DSSP's own legend).
    SS_CODES = {
        "H": "alpha-helix",
        "B": "beta-bridge",
        "E": "beta-strand",
        "G": "3-10 helix",
        "I": "pi-helix",
        "P": "polyproline helix",
        "T": "turn",
        "S": "bend",
        " ": "coil",
    }

    def analyze(self, structure: Structure, workdir: Path) -> Dict[str, Any]:
        """
        Run DSSP and return per-residue assignments plus summary stats.

        The output columns are DSSP's own; nothing is recomputed.
        """
        binary = self._binary()
        if not binary:
            reason = (f"configured DSSP binary not found at "
                      f"{self._binary_override}"
                      if self._binary_override
                      else "DSSP binary not found (apt install dssp or "
                           "LIGORA_DSSP_PATH)")
            return {
                "available": False,
                "error": reason,
            }

        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        pdb_path = workdir / "dssp_input.pdb"
        n_atoms = _write_protein_pdb(structure, pdb_path)
        if n_atoms == 0:
            return {
                "available": True,
                "error": "Structure has no polymer atoms to analyze",
            }

        out_path = workdir / "dssp_output.dssp"
        try:
            # DSSP 4.x CLI: mkdssp --output-format dssp <in> <out>
            proc = subprocess.run(
                [binary, "--output-format", "dssp",
                 str(pdb_path), str(out_path)],
                capture_output=True, timeout=300, text=True,
                cwd=str(workdir))
        except (OSError, subprocess.TimeoutExpired) as e:
            return {
                "available": True,
                "error": f"DSSP failed to run: {e}",
            }

        if proc.returncode != 0:
            return {
                "available": True,
                "error": (f"DSSP exited {proc.returncode}: "
                          f"{(proc.stderr or proc.stdout)[-300:]}"),
            }
        if not out_path.exists():
            return {"available": True, "error": "DSSP produced no output"}

        text = out_path.read_text(encoding="utf-8", errors="replace")
        residues: List[Dict[str, Any]] = []
        in_header = True
        for line in text.splitlines():
            if in_header:
                if line.startswith("  #  RESIDUE"):
                    in_header = False
                continue
            # DSSP residue line fixed columns (see its own legend):
            # resseq 1-5, chain 12, aa 14, ss 17, ASA 34-38
            if len(line) < 17:
                continue
            try:
                resseq = int(line[5:10])
            except ValueError:
                continue
            chain_id = line[11:12].strip() or " "
            aa = line[13:14].strip()
            ss = line[16:17]
            asa_token = line[34:38].strip() if len(line) >= 38 else ""
            try:
                asa = float(asa_token) if asa_token else None
            except ValueError:
                asa = None
            residues.append({
                "chain": chain_id,
                "residue_id": resseq,
                "aa": aa,
                "ss_code": ss,
                "ss": self.SS_CODES.get(ss, ss or "coil"),
                "accessible_surface": asa,
            })

        counts: Dict[str, int] = {}
        for r in residues:
            counts[r["ss"]] = counts.get(r["ss"], 0) + 1

        return {
            "available": True,
            "binary": binary,
            "residue_count": len(residues),
            "residues": residues,
            "ss_counts": counts,
            "source": "DSSP (Kabsch & Sander 1983; NKI 4.x)",
        }
