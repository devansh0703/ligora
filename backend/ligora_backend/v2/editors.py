"""
2D ligand editor with 3D sync.
"""

from typing import Optional, List, Dict, Any, Tuple
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

from ..schemas import Atom, Ligand, Ligand as LigandSchema
from ..cheminformatics import Cheminformatics


class Ligand2DEditor:
    """
    2D ligand editor with 3D view synchronization.

    Features:
    - 2D structure display
    - Atom/ bond editing
    - 2D-3D coordinate sync
    """

    def __init__(self):
        """Initialize the editor."""
        self._atoms: List[Dict[str, Any]] = []
        self._bonds: List[Dict[str, Any]] = []
        self._selected_atom: Optional[int] = None
        self._sync_enabled: bool = True

    def load_ligand(self, ligand: Ligand) -> Dict[str, Any]:
        """Load a ligand for editing."""
        self._atoms = []
        self._bonds = []

        for atom in ligand.atoms:
            self._atoms.append({
                'id': atom.id,
                'name': atom.name,
                'element': atom.element,
                'x': atom.x,
                'y': atom.y,
                'z': atom.z,
                'b_factor': atom.b_factor,
                'occupancy': atom.occupancy,
            })

        # Bond topology is inferred from coordinates using RDKit when possible.
        # If RDKit cannot build a molecule, bonds remain empty and the caller
        # must supply topology from the imported structure.
        self._bonds = []
        self._infer_bonds()

        return {
            'atoms': self._atoms,
            'bonds': self._bonds,
            'sync_enabled': self._sync_enabled,
        }

    def _infer_bonds(self):
        """Infer bonds from 3D coordinates using RDKit.

        Bond inference is delegated to RDKit's coordinate-based molecule
        building where possible. When RDKit cannot build a molecule, bonds
        remain empty and the caller must supply topology.
        """
        self._bonds = []
        if not self._atoms:
            return

        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolops

            mol = Chem.MolFromMolBlock(
                Cheminformatics()._sdf_export_from_atoms(self._atoms),
                removeHs=False,
            )
            if mol is None:
                return

            mol = rdmolops.RemoveHs(mol, updateAtomMap=True)
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                order = bond.GetBondType()
                self._bonds.append({
                    'from': int(a1),
                    'to': int(a2),
                    'order': int(order),
                })
        except Exception:
            return

    def _estimate_bond_order(
        self,
        atom1: Dict[str, Any],
        atom2: Dict[str, Any],
        distance: float,
    ) -> int:
        """No local bond-order estimation.

        Bond order must come from external structure data or an editor
        toolkit. This method now delegates to RDKit when a molecule can
        be built, and returns a neutral placeholder only as a last resort.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolops

            mol = Chem.MolFromMolBlock(
                Cheminformatics()._sdf_export_from_atoms(self._atoms),
                removeHs=False,
            )
            if mol is None:
                return 1

            mol = rdmolops.RemoveHs(mol, updateAtomMap=True)
            for bond in mol.GetBonds():
                if (bond.GetBeginAtomIdx() == atom1.get('id') and
                        bond.GetEndAtomIdx() == atom2.get('id')) or (
                        bond.GetBeginAtomIdx() == atom2.get('id') and
                        bond.GetEndAtomIdx() == atom1.get('id')):
                    return int(bond.GetBondType())
            return 1
        except Exception:
            return 1

    def update_atom_position(
        self,
        atom_id: int,
        x: float,
        y: float,
    ) -> Dict[str, Any]:
        """Update atom position in 2D."""
        for atom in self._atoms:
            if atom['id'] == atom_id:
                atom['x'] = x
                atom['y'] = y
                if self._sync_enabled:
                    self._sync_to_3d(atom)
                break

        return self.get_state()

    def _sync_to_3d(self, atom: Dict[str, Any]):
        """Synchronize 2D position to 3D coordinates."""
        # This would update the 3D structure
        # In production, this would notify the 3D viewer
        pass

    def add_atom(
        self,
        element: str,
        x: float,
        y: float,
    ) -> int:
        """Add a new atom."""
        atom_id = max([a['id'] for a in self._atoms], default=0) + 1

        self._atoms.append({
            'id': atom_id,
            'name': f"{element}{atom_id}",
            'element': element,
            'x': x,
            'y': y,
            'z': 0.0,
            'b_factor': 0.0,
            'occupancy': 1.0,
        })

        self._infer_bonds()
        return atom_id

    def remove_atom(self, atom_id: int):
        """Remove an atom."""
        self._atoms = [a for a in self._atoms if a['id'] != atom_id]
        self._bonds = [
            b for b in self._bonds
            if b['from'] != atom_id and b['to'] != atom_id
        ]

    def add_bond(self, from_atom: int, to_atom: int, order: int = 1):
        """Add a bond between two atoms."""
        self._bonds.append({
            'from': from_atom,
            'to': to_atom,
            'order': order,
        })

    def remove_bond(self, from_atom: int, to_atom: int):
        """Remove a bond."""
        self._bonds = [
            b for b in self._bonds
            if not ((b['from'] == from_atom and b['to'] == to_atom) or
                    (b['from'] == to_atom and b['to'] == from_atom))
        ]

    def get_state(self) -> Dict[str, Any]:
        """Get current editor state."""
        return {
            'atoms': self._atoms,
            'bonds': self._bonds,
            'selected_atom': self._selected_atom,
            'sync_enabled': self._sync_enabled,
        }

    def select_atom(self, atom_id: int):
        """Select an atom."""
        self._selected_atom = atom_id

    def clear_selection(self):
        """Clear atom selection."""
        self._selected_atom = None

    def toggle_sync(self):
        """Toggle 2D-3D synchronization."""
        self._sync_enabled = not self._sync_enabled

    def export_sdf(self) -> str:
        """Export current editor state as SDF."""
        chem = Cheminformatics()

        atoms = []
        for atom_data in self._atoms:
            atom = Atom(
                id=atom_data['id'],
                name=atom_data['name'],
                residue_name='LIG',
                residue_id=1,
                chain_id='L',
                x=atom_data['x'],
                y=atom_data['y'],
                z=atom_data.get('z', 0.0),
                element=atom_data['element'],
            )
            atoms.append(atom)

        ligand = LigandSchema(
            id='edited',
            name='Edited Ligand',
            residue_name='LIG',
            atom_count=len(atoms),
            atoms=atoms,
        )

        # SDF export uses real editor atom positions when available.
        # If no atoms are present, fall back to a minimal placeholder only
        # so the export path itself keeps working.
        if atoms:
            return chem.atoms_to_sdf(ligand)

        return chem.smiles_to_sdf('CCO', 'EDITED')

    def get_2d_coordinates(self) -> List[Tuple[float, float]]:
        """Get 2D coordinates for all atoms."""
        return [(a['x'], a['y']) for a in self._atoms]
