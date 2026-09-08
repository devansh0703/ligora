"""
Engine adapters for external simulation engines.

Engines are integrated, never reimplemented:
- AutoDock Vina (CLI): molecular docking. Ligand/receptor PDBQT preparation
  is delegated to Open Babel; charges come from Open Babel's Gasteiger
  computation (a real model), not an invented table.
- Open Babel: geometry cleanup (force-field minimization).
- gnina (CLI): deep-learning-guided docking when installed.

When an engine is not installed the app reports it as unavailable. Engine
failures are surfaced as errors - results are never simulated.
"""

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, List, Dict, Any

from .schemas import (
    Structure,
    Ligand,
    Atom,
    DockingPose,
    EngineType,
)
from .config import get_config


def _which(configured: Optional[str], name: str) -> Optional[str]:
    """Resolve an executable from explicit config or PATH."""
    if configured:
        return configured if Path(configured).exists() else None
    return shutil.which(name)


class VinaAdapter:
    """
    AutoDock Vina docking via the official command-line program.

    Input preparation (PDBQT conversion, torsion tree, Gasteiger charges)
    is performed by Open Babel. When Vina or Open Babel is missing, docking
    reports unavailability - no fake affinities are produced.
    """

    def __init__(self):
        self.config = get_config()
        self.engine_type = EngineType.VINA

    @property
    def name(self) -> str:
        return "vina"

    # ------------------------------------------------------------------
    # Availability
    # ------------------------------------------------------------------

    def _vina_path(self) -> Optional[str]:
        return _which(self.config.vina_executable, "vina")

    def _obabel_path(self) -> Optional[str]:
        return _which(self.config.obabel_executable, "obabel")

    def is_available(self) -> bool:
        return self._vina_path() is not None

    def preparation_available(self) -> bool:
        """PDBQT preparation requires Open Babel."""
        return self._obabel_path() is not None

    # ------------------------------------------------------------------
    # Input preparation (Open Babel)
    # ------------------------------------------------------------------

    def prepare(
        self,
        structure: Structure,
        ligand: Ligand,
        output_dir: Path,
        box_center: Optional[List[float]] = None,
        box_size: Optional[List[float]] = None,
        shared_receptor_pdbqt: Optional[Path] = None,
    ) -> Dict[str, Any]:
        """
        Prepare receptor and ligand PDBQT files with Open Babel.

        When shared_receptor_pdbqt is given and already exists it is
        reused directly (docking queues over one structure prepare the
        receptor exactly once); otherwise it is created there first.

        Raises RuntimeError with a clear message when preparation is not
        possible (Open Babel missing, ligand chemistry unknown, etc.).
        """
        obabel = self._obabel_path()
        if not obabel:
            raise RuntimeError(
                "Open Babel (obabel) is required for docking input "
                "preparation but is not installed")

        output_dir.mkdir(parents=True, exist_ok=True)
        ligand_sdf = output_dir / "ligand.sdf"
        ligand_pdbqt = output_dir / "ligand.pdbqt"

        # Receptor: polymer chains to PDB, then rigid PDBQT conversion.
        # Prepared once per structure and reused across ligands: the
        # receptor conversion is by far the slowest preparation step.
        if shared_receptor_pdbqt is not None:
            receptor_pdbqt = Path(shared_receptor_pdbqt)
            receptor_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        else:
            receptor_pdbqt = output_dir / "receptor.pdbqt"
        receptor_pdb = receptor_pdbqt.with_suffix(".pdb")
        if not (receptor_pdbqt.is_file()
                and receptor_pdbqt.stat().st_size > 0):
            self._write_receptor_pdb(structure, receptor_pdb)
            self._prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)

        # Ligand: 3D SDF via RDKit bond perception (real chemistry).
        from .cheminformatics import Cheminformatics
        chem = Cheminformatics()
        sdf_text = chem.sdf_export_from_atoms(ligand.atoms)
        if sdf_text is None:
            raise RuntimeError(
                f"Ligand {ligand.residue_name} has incomplete element data; "
                f"cannot build a valid 3D molecule for docking. Complete "
                f"element data from the CCD first.")
        ligand_sdf.write_text(sdf_text, encoding="utf-8")
        self._run_obabel(obabel, [
            str(ligand_sdf), "-O", str(ligand_pdbqt),
            "--gen3d",   # ensure 3D + torsion tree
            "-h",        # keep hydrogens for correct docking chemistry
        ])

        if box_center is None or box_size is None:
            box_center, box_size = self.compute_box(ligand)

        return {
            "receptor_pdbqt": str(receptor_pdbqt),
            "ligand_pdbqt": str(ligand_pdbqt),
            "box_center": box_center,
            "box_size": box_size,
        }

    def _prepare_receptor_pdbqt(self, receptor_pdb: Path,
                                receptor_pdbqt: Path):
        """Rigid receptor PDBQT conversion via Open Babel."""
        obabel = self._obabel_path()
        if not obabel:
            raise RuntimeError(
                "Open Babel (obabel) is required for docking input "
                "preparation but is not installed")
        self._run_obabel(obabel, [
            str(receptor_pdb), "-O", str(receptor_pdbqt),
            "-xr",  # rigid receptor
        ])

    def _run_obabel(self, obabel: str, args: List[str]):
        result = subprocess.run(
            [obabel] + args, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            raise RuntimeError(
                f"Open Babel failed: {(result.stderr or result.stdout)[-400:]}")
        if "0 molecules converted" in (result.stdout + result.stderr):
            raise RuntimeError("Open Babel converted 0 molecules")

    def _write_receptor_pdb(self, structure: Structure, path: Path):
        """Write polymer chains in PDB format (columns per PDB spec).

        ATOM columns: record 1-6, serial 7-11, name 13-16, altLoc 17,
        resName 18-20, chainID 22, resSeq 23-26, x/y/z 31-54,
        occupancy 55-60, tempFactor 61-66, element 77-78. No HEADER or
        other record types: strict parsers (Vina) reject unexpected tags
        in rigid receptor files.
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

    def prepare_receptor(self, structure: Structure, receptor_pdbqt: Path):
        """
        Prepare the rigid receptor PDBQT once for reuse across many
        ligand dockings (docking queue). Idempotent: skips conversion
        when the file already exists with content.
        """
        receptor_pdbqt = Path(receptor_pdbqt)
        if receptor_pdbqt.is_file() and receptor_pdbqt.stat().st_size > 0:
            return
        receptor_pdbqt.parent.mkdir(parents=True, exist_ok=True)
        receptor_pdb = receptor_pdbqt.with_suffix(".pdb")
        self._write_receptor_pdb(structure, receptor_pdb)
        self._prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)

    def compute_box(self, ligand: Ligand) -> tuple:
        """
        Docking box centered on the ligand with 8 A padding around its
        extent. The ligand pose defines the search space (co-crystallized
        ligand mode); box parameters remain user-overridable via payload.
        """
        import numpy as np
        if not ligand.atoms:
            raise ValueError("Cannot define a docking box without ligand "
                             "atoms")
        positions = np.array([[a.x, a.y, a.z] for a in ligand.atoms])
        center = np.mean(positions, axis=0)
        extent = np.ptp(positions, axis=0)
        size = extent + 8.0
        return ([round(float(c), 3) for c in center],
                [round(float(s), 3) for s in np.maximum(size, 8.0)])

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def run(self, input_files: Dict[str, Any],
            parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Run Vina docking with the official CLI."""
        vina = self._vina_path()
        if not vina:
            raise RuntimeError(
                "AutoDock Vina is not installed. Install it (e.g. 'apt "
                "install autodock-vina' or from ccsb.scripps.edu) to run "
                "docking.")

        out_dir = Path(parameters.get("output_dir") or tempfile.mkdtemp(
            prefix="ligora_vina_"))
        out_dir.mkdir(parents=True, exist_ok=True)
        out_pdbqt = out_dir / "out.pdbqt"
        log_file = out_dir / "vina.log"

        cmd = [
            vina,
            "--receptor", input_files["receptor_pdbqt"],
            "--ligand", input_files["ligand_pdbqt"],
            "--out", str(out_pdbqt),
            "--center_x", str(input_files["box_center"][0]),
            "--center_y", str(input_files["box_center"][1]),
            "--center_z", str(input_files["box_center"][2]),
            "--size_x", str(input_files["box_size"][0]),
            "--size_y", str(input_files["box_size"][1]),
            "--size_z", str(input_files["box_size"][2]),
            "--exhaustiveness", str(parameters.get(
                "exhaustiveness", self.config.default_docking_exhaustiveness)),
            "--num_modes", str(parameters.get(
                "num_modes", self.config.default_docking_num_modes)),
            "--energy_range", str(parameters.get(
                "energy_range", self.config.default_docking_energy_range)),
            "--verbosity", "2",
        ]

        import time
        start = time.time()
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True,
                timeout=self.config.vina_timeout)
        except subprocess.TimeoutExpired:
            raise RuntimeError("Vina docking timed out")
        runtime = time.time() - start

        log = result.stdout + result.stderr
        log_file.write_text(log, encoding="utf-8")

        if result.returncode != 0:
            raise RuntimeError(f"Vina failed: {log[-600:]}")

        poses = self._parse_output_pdbqt(out_pdbqt)
        return {
            "output_pdbqt": str(out_pdbqt),
            "log": log,
            "log_file": str(log_file),
            "runtime_seconds": runtime,
            "poses": poses,
        }

    def _parse_output_pdbqt(self, pdbqt_path: Path) -> List[DockingPose]:
        """Parse Vina PDBQT output into poses (REMARK VINA RESULT lines)."""
        poses: List[DockingPose] = []
        current: Optional[DockingPose] = None
        pose_id = 0

        for line in pdbqt_path.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT:") or \
                    line.startswith("REMARK  VINA RESULT:"):
                parts = line.split()
                try:
                    affinity = float(parts[3])
                except (IndexError, ValueError):
                    affinity = float("nan")
                if current is not None:
                    poses.append(current)
                pose_id += 1
                current = DockingPose(pose_id=pose_id, affinity=affinity)
            elif line.startswith(("ATOM", "HETATM")) and current is not None:
                # PDBQT fixed columns: name 13-16, resname 18-20,
                # x/y/z 31-54, partial charge 71-76, atom type 78-79.
                name = line[12:16].strip()
                resname = line[17:20].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                element = (line[77:79].strip() or
                           re.sub(r"[^A-Za-z]", "", line[12:16])[:1]).upper()
                current.atoms.append(Atom(
                    id=len(current.atoms) + 1,
                    name=name,
                    residue_name=resname,
                    residue_id=1,
                    chain_id="L",
                    x=x, y=y, z=z,
                    element=element if element else None,
                ))
        if current is not None:
            poses.append(current)
        return poses


class GninaAdapter:
    """
    gnina docking (CNN scoring) via its CLI, following the same data flow
    as the Vina adapter. Availability requires the gnina executable.
    """

    def __init__(self):
        self.config = get_config()
        self.engine_type = EngineType.GNINA

    @property
    def name(self) -> str:
        return "gnina"

    def _gnina_path(self) -> Optional[str]:
        return _which(self.config.gnina_executable, "gnina")

    def is_available(self) -> bool:
        return self._gnina_path() is not None

    def run(self, input_files: Dict[str, Any],
            parameters: Dict[str, Any]) -> Dict[str, Any]:
        gnina = self._gnina_path()
        if not gnina:
            raise RuntimeError("gnina is not installed")

        out_dir = Path(parameters.get("output_dir") or tempfile.mkdtemp(
            prefix="ligora_gnina_"))
        out_dir.mkdir(parents=True, exist_ok=True)
        out_pdbqt = out_dir / "out.pdbqt"

        cmd = [
            gnina,
            "--receptor", input_files["receptor_pdbqt"],
            "--ligand", input_files["ligand_pdbqt"],
            "--out", str(out_pdbqt),
            "--center_x", str(input_files["box_center"][0]),
            "--center_y", str(input_files["box_center"][1]),
            "--center_z", str(input_files["box_center"][2]),
            "--size_x", str(input_files["box_size"][0]),
            "--size_y", str(input_files["box_size"][1]),
            "--size_z", str(input_files["box_size"][2]),
            "--exhaustiveness", str(parameters.get("exhaustiveness", 8)),
            "--num_modes", str(parameters.get("num_modes", 9)),
        ]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode != 0:
            raise RuntimeError(
                f"gnina failed: {(result.stdout + result.stderr)[-600:]}")
        return {
            "output_pdbqt": str(out_pdbqt),
            "log": result.stdout + result.stderr,
            "poses": VinaAdapter._parse_output_pdbqt(self, out_pdbqt),
        }


class EngineRegistry:
    """Registry of engine adapters with honest availability reporting."""

    def __init__(self):
        self._adapters = {}

    def get_adapter(self, engine_type: EngineType):
        if engine_type not in self._adapters:
            if engine_type == EngineType.VINA:
                self._adapters[engine_type] = VinaAdapter()
            elif engine_type == EngineType.GNINA:
                self._adapters[engine_type] = GninaAdapter()
            else:
                raise ValueError(f"Unsupported engine: {engine_type}")
        return self._adapters[engine_type]

    def get_engine_status(self) -> Dict[str, Dict[str, Any]]:
        status = {}
        for engine_type in (EngineType.VINA, EngineType.GNINA):
            try:
                adapter = self.get_adapter(engine_type)
                status[engine_type.value] = {
                    "available": adapter.is_available(),
                    "engine_type": engine_type.value,
                }
            except ValueError:
                continue
        return status
