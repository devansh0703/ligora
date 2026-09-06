"""
2D ligand editor with 3D sync.
"""

from typing import Optional, List, Dict, Any, Tuple


from ..schemas import Atom, Ligand
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
                'residue_name': atom.residue_name,
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
        """Bond topology from real sources, tried in order:

        1. The CCD's own `_chem_comp_bond` records mapped by atom name
           (authoritative; works even when RDKit cannot perceive from
           coordinates).
        2. RDKit molecule building (perception + CCD records inside
           Cheminformatics), with H atoms removed from the editor view and
           heavy-heavy bond indices mapped through the mol-block atom
           ordering.
        Otherwise bonds stay empty and the caller supplies topology.
        """
        self._bonds = []
        if not self._atoms:
            return

        try:
            from rdkit import Chem
            chem = Cheminformatics()
            atoms = [Atom(id=a['id'], name=a['name'],
                          residue_name=a.get('residue_name') or 'LIG',
                          residue_id=1, chain_id='L',
                          x=a['x'], y=a['y'], z=a.get('z', 0.0),
                          element=a['element']) for a in self._atoms]
            mol_block = chem.mol_block_from_atoms(atoms)
            if mol_block is None:
                return

            mol = Chem.MolFromMolBlock(mol_block, removeHs=False,
                                       sanitize=True)
            if mol is None:
                return

            # Mol-block atom order: input atoms first (same order), then
            # any hydrogens appended by the CCD-bond path. Only heavy-atom
            # pairs within the input range map back to the editor view.
            n_input = len(self._atoms)

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                if i >= n_input or j >= n_input:
                    continue
                ai, aj = self._atoms[i], self._atoms[j]
                if ((ai.get('element') or '').upper() == 'H' or
                        (aj.get('element') or '').upper() == 'H'):
                    continue
                order = bond.GetBondTypeAsDouble()
                order_int = {1.0: 1, 1.5: 1, 2.0: 2, 3.0: 3}.get(
                    order, 1)
                self._bonds.append({
                    'from': ai['id'],
                    'to': aj['id'],
                    'order': order_int,
                })
        except Exception:
            self._bonds = []

    def _estimate_bond_order(
        self,
        atom1: Dict[str, Any],
        atom2: Dict[str, Any],
        distance: float,
    ) -> int:
        """Bond order from the inferred topology (RDKit/CCD derived).

        Returns the order already inferred for this atom pair, or 1 when the
        pair is not bonded. No per-distance heuristics are applied.
        """
        for bond in self._bonds:
            if ((bond['from'] == atom1.get('id') and
                 bond['to'] == atom2.get('id')) or
                (bond['from'] == atom2.get('id') and
                 bond['to'] == atom1.get('id'))):
                return bond['order']
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
        """Add a new atom (no automatic bond perception: partial geometry
        would yield phantom bonds)."""
        atom_id = max([a['id'] for a in self._atoms], default=0) + 1

        self._atoms.append({
            'id': atom_id,
            'name': f"{element}{atom_id}",
            'element': element,
            'residue_name': 'LIG',
            'x': x,
            'y': y,
            'z': 0.0,
            'b_factor': 0.0,
            'occupancy': 1.0,
        })

        return atom_id

    def remove_atom(self, atom_id: int):
        """Remove an atom."""
        self._atoms = [a for a in self._atoms if a['id'] != atom_id]
        self._bonds = [
            b for b in self._bonds
            if b['from'] != atom_id and b['to'] != atom_id
        ]

    def add_bond(self, from_atom: int, to_atom: int, order: int = 1):
        """Add a bond between two atoms (replaces an existing bond on the
        same pair rather than duplicating it)."""
        self.remove_bond(from_atom, to_atom)
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

        ligand = Ligand(
            id='edited',
            name='Edited Ligand',
            residue_name='LIG',
            atom_count=len(atoms),
            atoms=atoms,
        )

        # SDF export uses real editor atom positions when available.
        # When RDKit cannot build a molecule from the coordinates, the
        # fallback still writes an atom-only block rather than inventing
        # a SMILES from local element counts.
        return chem.atoms_to_sdf(ligand)

    def get_2d_coordinates(self) -> List[Tuple[float, float]]:
        """Get 2D coordinates for all atoms."""
        return [(a['x'], a['y']) for a in self._atoms]
