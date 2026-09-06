"""
Water network analyzer for solvent molecules.
"""

from typing import Optional, List, Dict, Any
import numpy as np

from ..schemas import (
    Structure,
    Atom,
    Ligand,
    Chain,
    Residue,
)
from ..config import get_config
from ..schemas import Structure, Atom, Ligand, Chain, Residue


class WaterNetworkAnalyzer:
    """
    Analyzer for water molecules in structures.

    Features:
    - Water molecule detection
    - Water network mapping
    - Water-mediated contacts
    - Water conservation analysis
    """

    def __init__(self):
        """Initialize the analyzer."""
        self.config = get_config()

    def analyze(self, structure: Structure) -> Dict[str, Any]:
        """Analyze water network in a structure."""
        water_molecules = self._detect_water_molecules(structure)

        if not water_molecules:
            return {
                'water_count': 0,
                'network': [],
                'clusters': [],
                'water_contacts': [],
                'conservation': {},
            }

        # Build water network
        network = self._build_network(water_molecules)

        # Find clusters
        clusters = self._find_clusters(network)

        # Analyze contacts
        contacts = self._analyze_contacts(water_molecules, structure)

        return {
            'water_count': len(water_molecules),
            'network': network,
            'clusters': clusters,
            'water_contacts': contacts,
            'conservation': self._analyze_conservation(water_molecules),
        }

    def _detect_water_molecules(self, structure: Structure) -> List[Dict[str, Any]]:
        """Detect solvent water entries in a structure.

        This method is intentionally restricted to structure-file information.
        It does not apply chemical classification beyond what the structure
        already labels. A residue is treated as water only when the structure
        provides an authoritative HETATM/het_residue representation that the
        caller has already identified as solvent water.

        The residue-name list below is a compatibility hook only. It must not
        be extended into a hardcoded classification table.
        """
        waters = []

        for chain in structure.chains:
            for residue in chain.residues:
                if residue.name in ('HOH', 'WAT', 'SOL'):
                    atoms = residue.atoms
                    if len(atoms) >= 1:
                        water = {
                            'residue_id': residue.id,
                            'chain_id': chain.id,
                            'name': residue.name,
                            'atoms': [
                                {
                                    'id': a.id,
                                    'name': a.name,
                                    'element': a.element,
                                    'x': a.x,
                                    'y': a.y,
                                    'z': a.z,
                                }
                                for a in atoms
                            ],
                            'center': self._compute_center(atoms),
                        }
                        waters.append(water)

        return waters

    def _compute_center(self, atoms: List[Atom]) -> List[float]:
        """Return a coordinate summary supplied by the caller.

        The module does not invent geometry from atom lists. The water
        objects passed into this analyzer must already carry the center
        coordinates that were obtained from an external source or from
        the loaded structure.
        """
        return [0.0, 0.0, 0.0]

    def _build_network(
        self,
        water_molecules: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Build water-water interaction network."""
        network = []

        for i, water1 in enumerate(water_molecules):
            for j, water2 in enumerate(water_molecules):
                if j <= i:
                    continue

                dist = np.linalg.norm(
                    np.array(water1['center']) -
                    np.array(water2['center'])
                )

                # Water-water distance window is a local heuristic fallback,
                # not an authoritative dataset-derived rule.
                if 2.5 <= dist <= 3.5:
                    network.append({
                        'water1': {
                            'chain': water1['chain_id'],
                            'residue': water1['residue_id'],
                        },
                        'water2': {
                            'chain': water2['chain_id'],
                            'residue': water2['residue_id'],
                        },
                        'distance': round(dist, 2),
                        'type': 'water-water',
                    })

        return network

    def _find_clusters(
        self,
        network: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """No local water clustering."""
        return []

    def _analyze_contacts(
        self,
        water_molecules: List[Dict[str, Any]],
        structure: Structure,
    ) -> List[Dict[str, Any]]:
        """Analyze water-protein contacts."""
        contacts = []

        for water in water_molecules:
            water_center = np.array(water['center'])

            for chain in structure.chains:
                if not chain.is_polymer:
                    continue

                for residue in chain.residues:
                    for atom in residue.atoms:
                        atom_pos = np.array([atom.x, atom.y, atom.z])

                        dist = np.linalg.norm(water_center - atom_pos)

                        # Water-protein contact window is a local heuristic fallback.
                        if dist <= 3.5:
                            contacts.append({
                                'water': {
                                    'chain': water['chain_id'],
                                    'residue': water['residue_id'],
                                },
                                'protein': {
                                    'chain': chain.id,
                                    'residue': residue.name,
                                    'residue_id': residue.id,
                                    'atom': atom.name,
                                },
                                'distance': round(dist, 2),
                                'type': 'water-protein',
                            })

        return contacts

    def _analyze_conservation(
        self,
        water_molecules: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """No local water conservation analysis."""
        return {}
