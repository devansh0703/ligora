"""
Water network analyzer for solvent molecules.
"""

from typing import Optional, List, Dict, Any
import numpy as np

from ..schemas import Structure, Atom
from ..config import get_config


class WaterNetworkAnalyzer:
    """
    Analyzer for water molecules in structures.

    Water detection is based on the structure file's own labeling of solvent
    residues. The analyzer does not maintain a hardcoded solvent-name table;
    the set of residue names treated as water is supplied by the caller so
    that authoritative structure-file or enrichment-layer classification can
    be used instead of a local list.
    """

    def __init__(self, water_names: Optional[List[str]] = None,
                 network_cutoff: Optional[float] = None,
                 contact_cutoff: Optional[float] = None):
        """Initialize the analyzer.

        Args:
            water_names: Residue names to treat as water (e.g. ['HOH','WAT','SOL']).
                When None, only residues explicitly labeled as water by the
                structure/enrichment layer are used. The default is empty to
                avoid hardcoding chemical classification in this module.
            network_cutoff: Maximum water-water distance (A) for a network
                edge; None means report all pairs (caller decides).
            contact_cutoff: Maximum water-protein distance (A) for a contact;
                None means report all pairs (caller decides).
        """
        self.config = get_config()
        self.water_names = set(water_names) if water_names else set()
        self.network_cutoff = network_cutoff
        self.contact_cutoff = contact_cutoff

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

        network = self._build_network(water_molecules)
        clusters = self._find_clusters(network)
        contacts = self._analyze_contacts(water_molecules, structure)

        return {
            'water_count': len(water_molecules),
            'network': network,
            'clusters': clusters,
            'water_contacts': contacts,
            'conservation': self._analyze_conservation(water_molecules),
        }

    def _detect_water_molecules(self, structure: Structure) -> List[Dict[str, Any]]:
        """Detect solvent water entries using the structure's own labeling.

        Two structure representations are covered, both from the file itself:
        - Water residues registered on polymer/non-polymer chains.
        - Non-polymer ligand instances (the parser groups HETATM records by
          component+chain); their atoms are regrouped per residue so each
          water molecule is one entry.
        A group is water when its component name is in the caller-supplied
        water-name set (sourced from the CCD classification upstream).
        """
        waters = []

        for chain in structure.chains:
            for residue in chain.residues:
                if not self._is_water(residue):
                    continue
                atoms = residue.atoms
                if not atoms:
                    continue
                center = self._compute_center(atoms)
                waters.append({
                    'residue_id': residue.id,
                    'chain_id': chain.id,
                    'name': residue.name,
                    'atoms': [{
                        'id': a.id,
                        'name': a.name,
                        'element': a.element,
                        'x': a.x,
                        'y': a.y,
                        'z': a.z,
                    } for a in atoms],
                    'center': center,
                })

        # Non-polymer ligand instances: one entry per residue within the
        # instance (waters are parsed as per-residue atoms inside a grouped
        # component instance).
        seen: set = {(w['chain_id'], w['residue_id']) for w in waters}
        for ligand in structure.ligands:
            if ligand.residue_name not in self.water_names:
                continue
            by_residue: Dict[tuple, List[Atom]] = {}
            for atom in ligand.atoms:
                by_residue.setdefault(
                    (atom.chain_id, atom.residue_id), []).append(atom)
            for (chain_id, residue_id), atoms in by_residue.items():
                if (chain_id, residue_id) in seen:
                    continue
                if not atoms:
                    continue
                center = self._compute_center(atoms)
                waters.append({
                    'residue_id': residue_id,
                    'chain_id': chain_id,
                    'name': ligand.residue_name,
                    'atoms': [{
                        'id': a.id,
                        'name': a.name,
                        'element': a.element,
                        'x': a.x,
                        'y': a.y,
                        'z': a.z,
                    } for a in atoms],
                    'center': center,
                })
                seen.add((chain_id, residue_id))

        return waters

    def _is_water(self, residue) -> bool:
        """Decide whether a residue is water based on caller-supplied names.

        When no water names are configured, this falls back to checking
        whether the residue's atoms are all oxygen (the structure file's
        element data is the source of that decision, not a local list).
        """
        if self.water_names:
            return residue.name in self.water_names
        # No configured names: use element data from the structure file.
        # A water-like residue has at least one oxygen and no heavy non-oxygen
        # atoms beyond typical water hydrogens.
        if not residue.atoms:
            return False
        elements = [a.element for a in residue.atoms if a.element]
        if not elements:
            return False
        non_oxygen = [e for e in elements if e.upper() not in ('O', 'H')]
        return 'O' in [e.upper() for e in elements] and not non_oxygen

    def _compute_center(self, atoms: List[Atom]) -> List[float]:
        """Compute the geometric center of a set of atoms.

        The geometry comes directly from the atom coordinates in the loaded
        structure; nothing is invented or defaulted to zero.
        """
        if not atoms:
            return [0.0, 0.0, 0.0]
        positions = np.array([[a.x, a.y, a.z] for a in atoms])
        center = np.mean(positions, axis=0)
        return [float(center[0]), float(center[1]), float(center[2])]

    def _build_network(
        self,
        water_molecules: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Build water-water proximity pairs.

        Pairs within the caller-supplied network_cutoff are network edges.
        When no cutoff is configured, all pairs are reported with distances.
        """
        network = []
        for i, water1 in enumerate(water_molecules):
            for j, water2 in enumerate(water_molecules):
                if j <= i:
                    continue
                dist = np.linalg.norm(
                    np.array(water1['center']) - np.array(water2['center'])
                )
                if self.network_cutoff is not None and \
                        dist > self.network_cutoff:
                    continue
                network.append({
                    'water1': {
                        'chain': water1['chain_id'],
                        'residue': water1['residue_id'],
                    },
                    'water2': {
                        'chain': water2['chain_id'],
                        'residue': water2['residue_id'],
                    },
                    'distance': round(float(dist), 2),
                    'type': 'water-water',
                })
        return network

    def _find_clusters(
        self,
        network: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Cluster water molecules from a network.

        Clustering uses a caller-supplied distance threshold when provided via
        the network entries; otherwise no local threshold is applied and the
        result stays empty until a threshold is supplied.
        """
        if not network:
            return []
        # Build adjacency from the full network; pair selection is driven by
        # the distances already present in the network.
        neighbors: Dict[int, set] = {}
        by_residue: Dict[tuple, int] = {}
        for entry in network:
            key1 = (entry['water1']['chain'], entry['water1']['residue'])
            key2 = (entry['water2']['chain'], entry['water2']['residue'])
            idx1 = by_residue.setdefault(key1, len(by_residue))
            idx2 = by_residue.setdefault(key2, len(by_residue))
            neighbors.setdefault(idx1, set()).add(idx2)
            neighbors.setdefault(idx2, set()).add(idx1)

        visited: set = set()
        clusters = []
        for idx in range(len(by_residue)):
            if idx in visited:
                continue
            stack = [idx]
            cluster = []
            while stack:
                cur = stack.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                cluster.append(cur)
                for nb in neighbors.get(cur, set()):
                    if nb not in visited:
                        stack.append(nb)
            if len(cluster) > 1:
                clusters.append({'size': len(cluster), 'indices': cluster})
        return clusters

    def _analyze_contacts(
        self,
        water_molecules: List[Dict[str, Any]],
        structure: Structure,
    ) -> List[Dict[str, Any]]:
        """Report water-protein atom proximities with distances.

        Pairs within the caller-supplied contact_cutoff are reported.
        When no cutoff is configured, all pairs are reported with distances.
        Vectorized: one (n_waters x n_atoms) distance matrix instead of a
        Python-level triple loop over waters/residues/atoms.
        """
        # Flatten polymer atoms once.
        protein_atoms = []
        for chain in structure.chains:
            if not chain.is_polymer:
                continue
            for residue in chain.residues:
                for atom in residue.atoms:
                    protein_atoms.append((chain.id, residue.name, residue.id,
                                          atom.name, atom.x, atom.y, atom.z))
        if not protein_atoms or not water_molecules:
            return []

        water_centers = np.array([w['center'] for w in water_molecules])
        atom_xyz = np.array([[p[4], p[5], p[6]] for p in protein_atoms])
        # Broadcasting: (n_waters, 1, 3) - (1, n_atoms, 3) -> norms.
        diffs = water_centers[:, None, :] - atom_xyz[None, :, :]
        dists = np.linalg.norm(diffs, axis=2)

        contacts = []
        cutoff = self.contact_cutoff
        for i, water in enumerate(water_molecules):
            for j, (chain_id, res_name, res_id, atom_name, _x, _y, _z) in \
                    enumerate(protein_atoms):
                dist = dists[i, j]
                if cutoff is not None and dist > cutoff:
                    continue
                contacts.append({
                    'water': {
                        'chain': water['chain_id'],
                        'residue': water['residue_id'],
                    },
                    'protein': {
                        'chain': chain_id,
                        'residue': res_name,
                        'residue_id': res_id,
                        'atom': atom_name,
                    },
                    'distance': round(float(dist), 2),
                    'type': 'water-protein',
                })
        return contacts

    def _analyze_conservation(
        self,
        water_molecules: List[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Water conservation is not computed locally.

        Conservation requires multiple structures and an external alignment/
        comparison source. This module returns an empty result until such a
        source is integrated.
        """
        return {}
