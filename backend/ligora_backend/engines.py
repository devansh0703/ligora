"""
Engine adapter framework and adapters for external simulation engines.

Provides:
- Common engine adapter interface
- Vina adapter for molecular docking
- gnina adapter (optional)
- Engine discovery and management
- Session/job orchestration for engine runs
"""

import subprocess
import tempfile
import shutil
import os
import json
import time
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple, Callable
from dataclasses import dataclass, field

import numpy as np

from .schemas import (
    Structure,
    Ligand,
    Atom,
    DockingPose,
    DockingResult,
    GeometryCleanupResult,
    EngineType,
    Job,
    JobStatus,
)
from .config import get_config


class EngineAdapter:
    """
    Base class for engine adapters.

    Each engine adapter implements:
    - prepare: Prepare input files for the engine
    - run: Execute the engine
    - parse_results: Parse engine output
    - cleanup: Clean up temporary files
    """

    def __init__(self, engine_type: EngineType):
        """Initialize the adapter."""
        self.engine_type = engine_type
        self.config = get_config()
        self._discovered = False
        self._executable_path: Optional[str] = None

    @property
    def name(self) -> str:
        """Get the engine name."""
        return self.engine_type.value

    @property
    def executable(self) -> Optional[str]:
        """Get the path to the engine executable."""
        if self._executable_path:
            return self._executable_path
        # Try to discover the executable
        discovered = self.discover()
        return self._executable_path if discovered else None

    def discover(self) -> bool:
        """
        Discover the engine executable on the system.

        Returns:
            True if the engine was found.
        """
        # Override in subclass
        self._discovered = False
        self._executable_path = None
        return False

    def is_available(self) -> bool:
        """Check if the engine is available and ready to use."""
        if not self._discovered:
            self.discover()
        return self._discovered and self._executable_path is not None

    def prepare(
        self,
        structure: Structure,
        ligand: Ligand,
        output_dir: Path,
    ) -> Dict[str, Any]:
        """
        Prepare input files for the engine.

        Args:
            structure: The structure.
            ligand: The ligand to dock.
            output_dir: Directory for output files.

        Returns:
            Dictionary with paths and parameters.
        """
        raise NotImplementedError

    def run(
        self,
        input_files: Dict[str, Any],
        parameters: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Run the engine.

        Args:
            input_files: Dictionary with input file paths.
            parameters: Engine parameters.

        Returns:
            Dictionary with output file paths and metadata.
        """
        raise NotImplementedError

    def parse_results(
        self,
        output_files: Dict[str, Any],
        parameters: Dict[str, Any],
    ) -> Any:
        """
        Parse engine output.

        Args:
            output_files: Dictionary with output file paths.
            parameters: Engine parameters.

        Returns:
            Parsed results.
        """
        raise NotImplementedError

    def cleanup(self, files: Dict[str, Any]):
        """Clean up temporary files."""
        for path in files.values():
            if path and isinstance(path, (str, Path)):
                try:
                    p = Path(path)
                    if p.is_file():
                        p.unlink()
                    elif p.is_dir():
                        shutil.rmtree(p)
                except OSError:
                    pass


class VinaAdapter(EngineAdapter):
    """
    Adapter for AutoDock Vina.

    Vina is an open-source molecular docking program.
    This adapter handles:
    - Preparing protein and ligand PDBQT files
    - Running Vina with appropriate parameters
    - Parsing Vina output (poses, affinities)
    """

    def __init__(self):
        """Initialize the Vina adapter."""
        super().__init__(EngineType.VINA)
        self._exhaustiveness = self.config.default_docking_exhaustiveness
        self._num_modes = self.config.default_docking_num_modes
        self._energy_range = self.config.default_docking_energy_range

    def discover(self) -> bool:
        """Discover AutoDock Vina executable."""
        # Check environment variable
        vina_path = os.environ.get('VINA_PATH')
        if vina_path and Path(vina_path).exists():
            self._executable_path = vina_path
            self._discovered = True
            return True

        # Check common installation paths
        possible_paths = [
            'vina',
            '/usr/bin/vina',
            '/usr/local/bin/vina',
            os.path.expanduser('~/miniconda3/bin/vina'),
            os.path.expanduser('~/anaconda3/bin/vina'),
            os.path.expanduser('~/.local/bin/vina'),
        ]

        for path in possible_paths:
            if path == 'vina':
                # Check if it's in PATH
                if shutil.which('vina'):
                    self._executable_path = 'vina'
                    self._discovered = True
                    return True
            elif Path(path).exists():
                self._executable_path = path
                self._discovered = True
                return True

        self._discovered = False
        return False

    def prepare(
        self,
        structure: Structure,
        ligand: Ligand,
        output_dir: Path,
    ) -> Dict[str, Any]:
        """
        Prepare files for Vina docking.

        Args:
            structure: Protein structure.
            ligand: Ligand to dock.
            output_dir: Output directory.

        Returns:
            Dictionary with paths to prepared files.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write protein PDBQT
        protein_pdbqt = output_dir / 'protein.pdbqt'
        self._write_protein_pdbqt(structure, protein_pdbqt)

        # Write ligand PDBQT
        ligand_pdbqt = output_dir / 'ligand.pdbqt'
        self._write_ligand_pdbqt(ligand, ligand_pdbqt)

        # Compute binding box from pocket
        box_info = self._compute_binding_box(structure, ligand)

        return {
            'protein_pdbqt': str(protein_pdbqt),
            'ligand_pdbqt': str(ligand_pdbqt),
            'box_center': box_info['center'],
            'box_size': box_info['size'],
        }

    def _write_protein_pdbqt(self, structure: Structure, output_path: Path):
        """Write protein in PDBQT format for Vina."""
        with open(output_path, 'w') as f:
            # Vina PDBQT format
            atom_serial = 1

            for chain in structure.chains:
                if not chain.is_polymer:
                    continue
                for residue in chain.residues:
                    for atom in residue.atoms:
                        # PDBQT format: ATOM serial name resName chainID resSeq x y z 0.0000 0.00 charge
                        elem = (atom.element or 'C').upper()

                        # Approximate partial charge
                        charge = self._estimate_charge(elem)

                        f.write(
                            f"ATOM  {atom_serial:5d}  {atom.name:<4}"
                            f"{atom.residue_name:<3}{atom.chain_id:<1}"
                            f"{atom.residue_id:4d}    "
                            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
                            f"  0.0000  0.00 {charge:8.3f}\n"
                        )
                        atom_serial += 1

    def _write_ligand_pdbqt(self, ligand: Ligand, output_path: Path):
        """Write ligand in PDBQT format for Vina."""
        with open(output_path, 'w') as f:
            atom_serial = 1

            for atom in ligand.atoms:
                elem = (atom.element or 'C').upper()
                charge = self._estimate_charge(elem)

                f.write(
                    f"HETATM{atom_serial:5d}  {atom.name:<4}"
                    f"{atom.residue_name:<3}{atom.chain_id:<1}"
                    f"{atom.residue_id:4d}    "
                    f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
                    f"  0.0000  0.00 {charge:8.3f}\n"
                )
                atom_serial += 1

    def _estimate_charge(self, element: str) -> float:
        """Estimate partial charge based on element (simplified)."""
        charges = {
            'C': 0.0, 'H': 0.0, 'N': -0.3, 'O': -0.4,
            'S': 0.0, 'P': 0.0, 'F': -0.2, 'CL': -0.2,
            'BR': -0.2, 'I': -0.2, 'NA': 1.0, 'K': 1.0,
            'CA': 2.0, 'MG': 2.0, 'ZN': 2.0,
        }
        return charges.get(element, 0.0)

    def _compute_binding_box(
        self,
        structure: Structure,
        ligand: Ligand,
    ) -> Dict[str, Any]:
        """Compute the docking box around the ligand."""
        if not ligand.atoms:
            return {
                'center': [0.0, 0.0, 0.0],
                'size': [20.0, 20.0, 20.0],
            }

        positions = np.array([
            [a.x, a.y, a.z] for a in ligand.atoms
        ])

        center = np.mean(positions, axis=0).tolist()

        # Box size: ligand extent + padding
        if len(positions) > 1:
            extent = np.ptp(positions, axis=0)
            size = (extent + 8.0).tolist()  # 8A padding
        else:
            size = [20.0, 20.0, 20.0]

        return {
            'center': [round(c, 3) for c in center],
            'size': [round(s, 3) for s in size],
        }

    def _make_config_file(self, output_dir: Path, parameters: Dict[str, Any]) -> Path:
        """Create Vina configuration file."""
        config_path = output_dir / ' Vina.conf'

        with open(config_path, 'w') as f:
            f.write(f"receptor = {parameters['protein_pdbqt']}\n")
            f.write(f"ligand = {parameters['ligand_pdbqt']}\n")
            f.write(f"out = {parameters.get('output', 'out.pdbqt')}\n")
            f.write(f"center_x = {parameters['box_center'][0]}\n")
            f.write(f"center_y = {parameters['box_center'][1]}\n")
            f.write(f"center_z = {parameters['box_center'][2]}\n")
            f.write(f"size_x = {parameters['box_size'][0]}\n")
            f.write(f"size_y = {parameters['box_size'][1]}\n")
            f.write(f"size_z = {parameters['box_size'][2]}\n")
            f.write(f"exhaustiveness = {parameters.get('exhaustiveness', self._exhaustiveness)}\n")
            f.write(f"num_modes = {parameters.get('num_modes', self._num_modes)}\n")
            f.write(f"energy_range = {parameters.get('energy_range', self._energy_range)}\n")

        return config_path

    def run(
        self,
        input_files: Dict[str, Any],
        parameters: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Run Vina docking.

        Args:
            input_files: Dictionary with input file paths.
            parameters: Vina parameters.

        Returns:
            Dictionary with output paths and metadata.
        """
        if not self.is_available():
            raise RuntimeError("Vina not available")

        output_dir = Path(parameters.get('output_dir', tempfile.gettempdir())) / 'vina_run'
        output_dir.mkdir(parents=True, exist_ok=True)

        # Prepare config
        config_params = {
            'protein_pdbqt': input_files['protein_pdbqt'],
            'ligand_pdbqt': input_files['ligand_pdbqt'],
            'output': str(output_dir / 'out.pdbqt'),
            'box_center': input_files['box_center'],
            'box_size': input_files['box_size'],
            'exhaustiveness': parameters.get('exhaustiveness', self._exhaustiveness),
            'num_modes': parameters.get('num_modes', self._num_modes),
            'energy_range': parameters.get('energy_range', self._energy_range),
        }

        config_file = self._make_config_file(output_dir, config_params)

        # Run Vina
        start_time = time.time()

        cmd = [
            self.executable,
            '--config', str(config_file),
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour max
        )

        runtime = time.time() - start_time

        # Parse Vina log output
        log = result.stdout + result.stderr

        # Extract affinity from Vina output
        affinity = None
        n_poses = 0
        for line in log.splitlines():
            if 'Affinity:' in line:
                parts = line.split()
                try:
                    affinity = float(parts[parts.index('Affinity:') + 1])
                except (ValueError, IndexError):
                    pass
            if 'mode' in line.lower() and 'affinity' in line.lower():
                n_poses += 1

        return {
            'output_pdbqt': str(output_dir / 'out.pdbqt'),
            'log': log,
            'runtime_seconds': runtime,
            'return_code': result.returncode,
            'Affinity': affinity,
            'num_poses': n_poses,
        }

    def parse_results(
        self,
        output_files: Dict[str, Any],
        parameters: Dict[str, Any],
    ) -> DockingResult:
        """
        Parse Vina docking results.

        Args:
            output_files: Dictionary with output file paths.
            parameters: Vina parameters.

        Returns:
            DockingResult with poses.
        """
        poses = []

        # Parse the output PDBQT file
        output_pdbqt = output_files.get('output_pdbqt') or parameters.get('output_pdbqt')
        if output_pdbqt and Path(output_pdbqt).exists():
            poses = self._parse_vina_output(output_pdbqt)

        # If no poses from file, parse affinity from log
        if not poses:
            log = output_files.get('log', '')
            affinities = []
            for line in log.splitlines():
                if 'Affinity:' in line:
                    try:
                        affinities.append(float(line.split()[1]))
                    except (ValueError, IndexError):
                        pass
            
            if affinities:
                # Create minimal poses from log affinities
                for i, aff in enumerate(affinities[:9], 1):
                    poses.append(DockingPose(
                        pose_id=i,
                        affinity=aff,
                        ligand_name=parameters.get('ligand_name', 'unknown'),
                    ))

        return DockingResult(
            job_id=parameters.get('job_id', 'unknown'),
            engine=EngineType.VINA,
            poses=poses,
            box_center=output_files.get('box_center'),
            box_size=output_files.get('box_size'),
            runtime_seconds=output_files.get('runtime_seconds', 0),
            output_file=output_files.get('output_pdbqt'),
        )

    def _parse_vina_output(self, pdbqt_path: str) -> List[DockingPose]:
        """Parse Vina PDBQT output file."""
        poses = []
        current_pose = None
        pose_id = 0

        with open(pdbqt_path) as f:
            for line in f:
                if line.startswith('REMARK  VINA RESULT:'):
                    # Extract affinity
                    parts = line.split()
                    try:
                        affinity = float(parts[3])
                    except (ValueError, IndexError):
                        affinity = 0.0

                    if current_pose is not None:
                        poses.append(current_pose)

                    pose_id += 1
                    current_pose = DockingPose(
                        pose_id=pose_id,
                        affinity=affinity,
                        ligand_name='docked_ligand',
                    )
                elif line.startswith('ATOM') or line.startswith('HETATM'):
                    if current_pose is not None:
                        parts = line.split()
                        atom = Atom(
                            id=int(parts[1]) if len(parts) > 1 else 0,
                            name=parts[2] if len(parts) > 2 else '',
                            residue_name=parts[3] if len(parts) > 3 else '',
                            residue_id=int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 0,
                            chain_id=parts[5] if len(parts) > 5 else 'L',
                            x=float(parts[6]) if len(parts) > 6 else 0.0,
                            y=float(parts[7]) if len(parts) > 7 else 0.0,
                            z=float(parts[8]) if len(parts) > 8 else 0.0,
                            element=parts[11] if len(parts) > 11 else None,
                        )
                        current_pose.atoms.append(atom)

            if current_pose is not None:
                poses.append(current_pose)

        return poses


class GninaAdapter(EngineAdapter):
    """
    Adapter for gnina (deep-learning guided docking).

    gnina is a fork of Vina with CNN scoring.
    """

    def __init__(self):
        """Initialize the gnina adapter."""
        super().__init__(EngineType.GNINA)

    def discover(self) -> bool:
        """Discover gnina executable."""
        gnina_path = os.environ.get('GNINA_PATH')
        if gnina_path and Path(gnina_path).exists():
            self._executable_path = gnina_path
            self._discovered = True
            return True

        # Check common paths
        possible_paths = [
            'gnina',
            '/usr/bin/gnina',
            '/usr/local/bin/gnina',
        ]

        for path in possible_paths:
            if path == 'gnina':
                if shutil.which('gnina'):
                    self._executable_path = 'gnina'
                    self._discovered = True
                    return True
            elif Path(path).exists():
                self._executable_path = path
                self._discovered = True
                return True

        self._discovered = False
        return False

    def run(
        self,
        input_files: Dict[str, Any],
        parameters: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Run gnina docking."""
        if not self.is_available():
            raise RuntimeError("gnina not available")

        output_dir = Path(parameters.get('output_dir', tempfile.gettempdir())) / 'gnina_run'
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.executable,
            '--receptor', input_files['protein_pdbqt'],
            '--ligand', input_files['ligand_pdbqt'],
            '--center_x', str(input_files['box_center'][0]),
            '--center_y', str(input_files['box_center'][1]),
            '--center_z', str(input_files['box_center'][2]),
            '--size_x', str(input_files['box_size'][0]),
            '--size_y', str(input_files['box_size'][1]),
            '--size_z', str(input_files['box_size'][2]),
            '--exhaustiveness', str(parameters.get('exhaustiveness', 8)),
            '--num_modes', str(parameters.get('num_modes', 9)),
            '--out', str(output_dir / 'out.pdbqt'),
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200,
        )

        return {
            'output_pdbqt': str(output_dir / 'out.pdbqt'),
            'log': result.stdout + result.stderr,
            'return_code': result.returncode,
        }


class EngineRegistry:
    """
    Registry for engine adapters.

    Manages discovery and availability of all engines.
    """

    def __init__(self):
        """Initialize the registry."""
        self._adapters: Dict[EngineType, EngineAdapter] = {}
        self._discovered = False

    def get_adapter(self, engine_type: EngineType) -> Optional[EngineAdapter]:
        """Get or create an adapter for an engine type."""
        if engine_type not in self._adapters:
            adapter = self._create_adapter(engine_type)
            self._adapters[engine_type] = adapter
        return self._adapters[engine_type]

    def _create_adapter(self, engine_type: EngineType) -> EngineAdapter:
        """Create an adapter for an engine type."""
        if engine_type == EngineType.VINA:
            return VinaAdapter()
        elif engine_type == EngineType.GNINA:
            return GninaAdapter()
        else:
            raise ValueError(f"Unknown engine type: {engine_type}")

    def discover_all(self) -> Dict[str, bool]:
        """
        Discover all registered engines.

        Returns:
            Dictionary of engine name -> available.
        """
        availability = {}
        for engine_type in [EngineType.VINA, EngineType.GNINA]:
            adapter = self.get_adapter(engine_type)
            available = adapter.discover()
            availability[engine_type.value] = available
        self._discovered = True
        return availability

    def is_available(self, engine_type: EngineType) -> bool:
        """Check if an engine is available."""
        adapter = self.get_adapter(engine_type)
        if not self._discovered:
            self.discover_all()
        return adapter.is_available()

    def get_available_engines(self) -> List[EngineType]:
        """Get list of available engine types."""
        if not self._discovered:
            self.discover_all()
        return [etype for etype in [EngineType.VINA, EngineType.GNINA]
                if self.is_available(etype)]

    def get_engine_status(self) -> Dict[str, Dict[str, Any]]:
        """
        Get status of all engines.

        Returns:
            Dictionary of engine name -> status info.
        """
        if not self._discovered:
            self.discover_all()

        status = {}
        for engine_type in [EngineType.VINA, EngineType.GNINA]:
            adapter = self.get_adapter(engine_type)
            status[engine_type.value] = {
                'available': adapter.is_available(),
                'executable': adapter.executable,
                'engine_type': engine_type.value,
            }
        return status
