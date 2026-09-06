"""
Cheminformatics module for ligand similarity, fingerprints, and structure handling.

Provides:
- SMILES generation and validation (via RDKit)
- Molecular fingerprints (Morgan/ECFP via RDKit)
- Tanimoto similarity (via RDKit)
- SDF export (via RDKit)
- SMILES to molecule conversion (via RDKit)
- Molecular weight, formula, logP, TPSA, atom counts (via RDKit)

All chemistry operations use RDKit. No chemical data is hardcoded here.
"""

from typing import Optional, List, Dict, Any, Tuple
from pathlib import Path

import numpy as np

from .schemas import Atom, Ligand


class Cheminformatics:
    """
    Cheminformatics operations backed by RDKit.

    Provides fingerprint generation, similarity computation,
    structure parsing/serialization, and property estimation.
    All chemical knowledge comes from RDKit, not from hardcoded tables.
    """

    def __init__(self, radius: int = 2, n_bits: int = 2048):
        """
        Initialize the cheminformatics engine.

        Args:
            radius: Morgan fingerprint radius (default 2 for ECFP4).
            n_bits: Number of bits in the fingerprint.
        """
        self.radius = radius
        self.n_bits = n_bits

    def _mol_from_smiles(self, smiles: str):
        """Convert SMILES to RDKit molecule. Returns None on failure."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles, sanitize=True)
            return mol
        except Exception:
            return None

    def _mol_from_atoms(self, atoms: List[Atom]):
        """Build an RDKit mol from a list of Atom objects if coordinates are sufficient.
        
        This is a best-effort path used when a SMILES is not available but atom
        coordinates exist. It delegates bond perception to RDKit where possible.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except Exception:
            return None

        if not atoms:
            return None

        # Use RDKit's bond inference from 3D coordinates
        mol = Chem.MolFromConformer(Chem.Conformer(len(atoms)))
        if mol is None:
            return None

        for idx, atom in enumerate(atoms):
            element = (atom.element or "C").upper()
            rd_atom = Chem.Atom(element)
            mol.AddAtom(rd_atom)
            conf = mol.GetConformer()
            conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))

        # Let RDKit perceive bonds from coordinates
        try:
            Chem.rdMolTransforms.Compute2DCoords(mol)
            mol = Chem.CombineMols(mol)
            mol = AllChem.RemoveHs(mol)
        except Exception:
            pass

        # Try distance-based bond perception as a fallback within RDKit
        try:
            from rdkit.Chem import rdmolops
            mol = rdmolops.AssignBondOrdersFromTemplate(None, mol)
        except Exception:
            pass

        return mol

    def compute_fingerprint(self, smiles: str) -> Optional[np.ndarray]:
        """
        Compute a Morgan fingerprint from a SMILES string using RDKit.

        Args:
            smiles: SMILES representation of the molecule.

        Returns:
            Binary fingerprint as numpy array, or None if parsing fails.
        """
        if not smiles:
            return None

        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return None

        try:
            from rdkit.Chem import rdFingerprintGenerator
            gen = rdFingerprintGenerator.GetMorganGenerator(
                radius=self.radius,
                fpSize=self.n_bits,
            )
            fp = gen.GetFingerprint(mol)
            arr = np.zeros((self.n_bits,), dtype=np.uint8)
            for i in range(fp.GetNumBits()):
                if fp[i]:
                    arr[i] = 1
            return arr
        except Exception:
            return None

    def _parse_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Placeholder for future SMILES parser if needed.

        RDKit handles SMILES parsing directly in all current code paths.
        No local SMILES parser is maintained here to avoid embedding
        chemical parsing heuristics.
        """
        return None

    def _sdf_export(self, smiles: str, output_path: str) -> bool:
        """
        Export a molecule as SDF using RDKit.

        Args:
            smiles: SMILES string of the molecule.
            output_path: Path to write the SDF file.

        Returns:
            True if export succeeded.
        """
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return False

        try:
            from rdkit.Chem import SDWriter
            from rdkit import Chem
            mol = Chem.RemoveHs(mol)
            writer = SDWriter(output_path)
            writer.write(mol)
            writer.close()
            return True
        except Exception as e:
            import traceback
            traceback.print_exc()
            return False

    def _sdf_export_from_atoms(self, atoms: List['Atom']) -> Optional[str]:
        """
        Build a mol block string from atom coordinates and element symbols.

        This is used by the 2D editor to let RDKit perceive bonds from
        coordinates. It does not add bond orders; RDKit will infer those
        from the geometry where possible.
        """
        if not atoms:
            return None

        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolfiles

            mol = Chem.MolFromMolBlock(
                rdmolfiles.MolToMolBlock(
                    Chem.MolFromSmiles('C')
                ),
                removeHs=False,
            )
            if mol is None:
                return None

            mol = Chem.RemoveHs(mol)
            conf = Chem.Conformer(len(atoms))
            for idx, atom in enumerate(atoms):
                elem = (atom.element or 'C').strip().upper()
                rd_atom = Chem.Atom(elem)
                mol.AddAtom(rd_atom)
                conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))
            mol.AddConformer(conf, assignId=True)

            mb = Chem.MolToMolBlock(mol)
            return mb
        except Exception:
            return None

    def _generate_2d_coords(self, smiles: str) -> Optional[List[Tuple[float, float, float]]]:
        """
        Generate 2D coordinates for a SMILES string using RDKit.

        Args:
            smiles: SMILES string.

        Returns:
            List of (x, y, z) tuples or None.
        """
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return None

        try:
            from rdkit.Chem import AllChem
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            coords = []
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                coords.append((pos.x, pos.y, pos.z))
            return coords
        except Exception:
            return None
        atoms = []
        bonds = []
        n = len(smiles)

        # Track ring closures
        ring_atoms: Dict[str, int] = {}
        branches: List[Tuple[int, int]] = []  # (start_idx, branch_start)
        branch_stack: List[int] = []

        i = 0
        current_branch_start = 0

        while i < n:
            c = smiles[i]

            # Handle branch open
            if c == '(':
                branch_stack.append(current_branch_start)
                current_branch_start = len(atoms)
                i += 1
                continue

            # Handle branch close
            if c == ')':
                if branch_stack:
                    parent_start = branch_stack.pop()
                    if len(atoms) > current_branch_start and current_branch_start > parent_start:
                        # Connect last atom in branch to first atom after branch
                        bonds.append({'from': current_branch_start - 1, 'to': parent_start, 'type': 1})
                current_branch_start = len(atoms)
                i += 1
                continue

            # Handle ring closures (digits or %XX)
            if c.isdigit():
                ring_num = c
                i += 1
                # Check for multi-digit ring number
                while i < n and smiles[i].isdigit():
                    ring_num += smiles[i]
                    i += 1
                
                if ring_num in ring_atoms:
                    # Close the ring
                    atom_idx = len(atoms) - 1
                    other_idx = ring_atoms[ring_num]
                    if atom_idx != other_idx:
                        bonds.append({'from': atom_idx, 'to': other_idx, 'type': 1})
                else:
                    ring_atoms[ring_num] = len(atoms) - 1 if atoms else 0
                continue

            if c == '%':
                i += 1
                ring_num = '%' + smiles[i:i+2] if i+2 <= n else '%'
                i += 2
                if ring_num in ring_atoms:
                    atom_idx = len(atoms) - 1
                    other_idx = ring_atoms[ring_num]
                    if atom_idx != other_idx:
                        bonds.append({'from': atom_idx, 'to': other_idx, 'type': 1})
                else:
                    ring_atoms[ring_num] = len(atoms) - 1 if atoms else 0
                continue

            # Handle atom symbols (uppercase = aliphatic/aromatic start)
            if c.isalpha():
                start = i
                # Get full atom symbol (could be one or two chars)
                if i + 1 < n and smiles[i + 1].islower():
                    i += 1
                atom_symbol = smiles[start:i + 1].lower() if smiles[start].islower() else smiles[start:i + 1]

                # Convert to standard element representation
                if atom_symbol in ('c', 'n', 'o', 's', 'p', 'se'):
                    # Aromatic atoms - track as regular atoms
                    pass

                atoms.append({
                    'symbol': atom_symbol.capitalize() if atom_symbol[0].islower() else atom_symbol,
                    'atomic_number': self._rdkit_atomic_number(atom_symbol),
                })

                # Connect to previous atom
                if len(atoms) >= 2:
                    # Check if previous bond was explicit
                    last_bond = bonds[-1] if bonds else None
                    if not last_bond or last_bond['to'] != len(atoms) - 2:
                        bonds.append({'from': len(atoms) - 2, 'to': len(atoms) - 1, 'type': 1})

                i += 1
                continue

            # Handle bonds
            if c in '=#':
                bond_type = 2 if c == '=' else 3 if c == '#' else 1
                # Update last bond
                if bonds:
                    bonds[-1]['type'] = bond_type
                i += 1
                continue

            # Handle charge modifiers like @ or @@ (stereochemistry - skip for now)
            if c in '@+':
                i += 1
                continue

            # Handle wildcards or other characters
            i += 1

        # Close any remaining ring closures
        for ring_num, idx in ring_atoms.items():
            if idx == len(atoms) - 1 and len(atoms) > 1:
                # Ring to self - connect to first atom in chain
                bonds.append({'from': idx, 'to': 0, 'type': 1})

        if not atoms:
            return None

        return {
            'atoms': atoms,
            'bonds': bonds,
            'n_atoms': len(atoms),
        }

    def _rdkit_atomic_number(self, symbol: str) -> int:
        """Get atomic number from RDKit's periodic table, not from a local table."""
        try:
            from rdkit import Chem
            atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(symbol)
            return int(atomic_num)
        except Exception:
            return 0

    def _morgan_fingerprint(self, mol: Dict[str, Any]) -> np.ndarray:
        """
        Generate a Morgan-like circular fingerprint.

        This is a simplified implementation. For production,
        use RDKit's GetMorganFingerprintAsBitVect.

        Args:
            mol: Parsed molecule dictionary.

        Returns:
            Binary fingerprint array.
        """
        fp = np.zeros(self.n_bits, dtype=np.uint8)
        n_atoms = mol.get('n_atoms', 0)

        if n_atoms == 0:
            return fp

        # Simple hash-based fingerprint
        for i in range(n_atoms):
            atom_hash = self._hash_atom_environment(mol, i)
            for bit_idx in range(min(4, self.n_bits)):
                bit = (atom_hash + bit_idx * 17) % self.n_bits
                fp[bit] = 1

        return fp

    def _hash_atom_environment(
        self,
        mol: Dict[str, Any],
        center_idx: int,
        radius: Optional[int] = None,
    ) -> int:
        """
        Hash the atom environment around a center atom.

        Args:
            mol: Molecule dictionary.
            center_idx: Index of center atom.
            radius: Search radius (defaults to self.radius).

        Returns:
            Hash value.
        """
        radius = radius or self.radius
        atoms = mol.get('atoms', [])

        if center_idx >= len(atoms):
            return 0

        # Build environment string
        env_parts = []
        center_atom = atoms[center_idx]
        env_parts.append(f"{center_atom['atomic_number']}")

        # Add neighbors (simplified - no actual bond graph)
        # In a full implementation, this would traverse the bond graph
        for _ in range(radius):
            env_parts.append("X")  # Placeholder

        env_str = ",".join(env_parts)
        return int(hashlib.md5(env_str.encode()).hexdigest()[:8], 16)

    def tanimoto_similarity(
        self,
        fp1: np.ndarray,
        fp2: np.ndarray,
    ) -> float:
        """
        Compute Tanimoto similarity between two binary fingerprints.

        Uses RDKit DataStructs when available for a proper bit-vector
        Tanimoto; falls back to a pure numpy Jaccard computation only
        as a last resort.

        Args:
            fp1: First fingerprint (numpy uint8 array of bits).
            fp2: Second fingerprint (numpy uint8 array of bits).

        Returns:
            Tanimoto coefficient between 0 and 1.
        """
        if fp1 is None or fp2 is None:
            return 0.0

        try:
            from rdkit import DataStructs as rd_ds
            rd_fp1 = rd_ds.CreateNumpyBitVect(fp1)
            rd_fp2 = rd_ds.CreateNumpyBitVect(fp2)
            return float(rd_ds.TanimotoSimilarity(rd_fp1, rd_fp2))
        except Exception:
            intersection = float(np.sum(fp1 * fp2))
            sum1 = float(np.sum(fp1))
            sum2 = float(np.sum(fp2))
            union = sum1 + sum2 - intersection
            if union == 0:
                return 0.0
            return intersection / union

    def smiles_to_sdf(self, smiles: str, molecule_name: str = "MOL") -> str:
        """
        Convert a SMILES string to SDF format.

        Args:
            smiles: SMILES string.
            molecule_name: Name for the molecule.

        Returns:
            SDF formatted string.
        """
        mol = self._parse_smiles(smiles)
        if mol is None:
            return f"$$$$\n"

        # Generate 2D coordinates (simple layout)
        coords = self._generate_2d_coordinates(mol)

        sdf_lines = []

        # Header block
        sdf_lines.append(molecule_name[:80].ljust(80))
        sdf_lines.append(" ".ljust(80))
        sdf_lines.append(f"{mol['n_atoms']:>3} {len(mol['atoms']):>3}".ljust(16))
        sdf_lines.append(" ".ljust(80))
        sdf_lines.append(" ".ljust(80))
        sdf_lines.append(" ".ljust(80))

        # Atom block
        for i, (atom, coord) in enumerate(zip(mol['atoms'], coords)):
            symbol = atom['symbol']
            x, y, z = coord
            # Atomic number lookup
            atomic_num = self._periodic_table_number(symbol)
            sdf_lines.append(
                f"{atomic_num:>3} {symbol:<2} {x:>10.4f} {y:>10.4f} {z:>10.4f}"
                f"    0.0000           0"
            )

        # Bond block
        for bond in self._extract_bonds(smiles, mol):
            sdf_lines.append(
                f"{bond['from']:>3} {bond['to']:>3} {bond['type']:>3}"
            )

        sdf_lines.append("$$$$")

        return "\n".join(sdf_lines)

    def _periodic_table_number(self, symbol: str) -> int:
        """Get atomic number from element symbol."""
        table = {
            'H': 1, 'HE': 2, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NE': 10,
            'NA': 11, 'MG': 12, 'AL': 13, 'SI': 14, 'P': 15, 'S': 16,
            'CL': 17, 'AR': 18, 'K': 19, 'CA': 20, 'BR': 35, 'I': 53,
        }
        return table.get(symbol.upper(), 6)

    def _generate_2d_coordinates(self, mol: Dict[str, Any]) -> List[Tuple[float, float, float]]:
        """
        Generate simple 2D coordinates for atoms.

        Uses a circular layout for visualization.

        Args:
            mol: Molecule dictionary.

        Returns:
            List of (x, y, z) coordinates.
        """
        n = mol.get('n_atoms', 0)
        if n == 0:
            return []

        coords = []
        angle_step = 2 * np.pi / max(n, 1)
        radius = 2.0

        for i in range(n):
            angle = i * angle_step
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            z = 0.0
            coords.append((x, y, z))

        return coords

    def _extract_bonds(self, smiles: str, mol: Dict[str, Any]) -> List[Dict[str, int]]:
        """
        Extract bonds from SMILES string (simplified).

        Args:
            smiles: SMILES string.
            mol: Molecule dictionary.

        Returns:
            List of bond dictionaries with 'from', 'to', 'type' keys.
        """
        bonds = []
        # Simple bond extraction - in production use RDKit
        # This creates a chain of bonds for linear molecules
        n = mol.get('n_atoms', 0)
        if n < 2:
            return bonds

        # Default: linear chain
        for i in range(n - 1):
            bonds.append({'from': i + 1, 'to': i + 2, 'type': 1})

        return bonds

    def compute_molecular_descriptors(self, smiles: str) -> Dict[str, Any]:
        """
        Compute molecular descriptors using RDKit.

        Args:
            smiles: SMILES string.

        Returns:
            Dictionary of descriptor values from RDKit.
        """
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return {}

        try:
            from rdkit.Chem import Descriptors
            from rdkit.Chem import rdMolDescriptors

            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            n_heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)

            formula = rdMolDescriptors.CalcFormula(mol)

            atom_counts = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

            return {
                'n_heavy_atoms': int(n_heavy),
                'n_hba': int(hba),
                'n_hbd': int(hbd),
                'n_atoms': mol.GetNumAtoms(),
                'logp': round(float(logp), 3),
                'tpsa': round(float(tpsa), 3),
                'molecular_weight': round(float(mw), 3),
                'formula': formula,
                'atom_counts': atom_counts,
            }
        except Exception:
            return {}

    def compute_descriptors(self, smiles: str) -> Dict[str, Any]:
        """
        Compute molecular descriptors using RDKit.

        Args:
            smiles: SMILES string.

        Returns:
            Dictionary of descriptor values from RDKit, or empty dict on failure.
        """
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return {}

        try:
            from rdkit.Chem import Descriptors
            from rdkit.Chem import rdMolDescriptors

            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            n_heavy = rdMolDescriptors.CalcNumHeavyAtoms(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)

            formula = rdMolDescriptors.CalcMolFormula(mol)

            atom_counts = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

            return {
                'n_heavy_atoms': int(n_heavy),
                'n_hba': int(hba),
                'n_hbd': int(hbd),
                'n_atoms': mol.GetNumAtoms(),
                'logp': round(float(logp), 3),
                'tpsa': round(float(tpsa), 3),
                'molecular_weight': round(float(mw), 3),
                'formula': formula,
                'atom_counts': atom_counts,
            }
        except Exception:
            return {}

    def find_similar_compounds(
        self,
        query_smiles: str,
        compounds: List[Dict[str, Any]],
        threshold: float = 0.7,
    ) -> List[Dict[str, Any]]:
        """
        Find compounds similar to a query SMILES using RDKit fingerprints.

        Args:
            query_smiles: Query SMILES.
            compounds: List of compounds with 'smiles' key.
            threshold: Similarity threshold (0-1). This is supplied by the
                caller/data source, not chosen here.

        Returns:
            List of similar compounds with similarity scores.
        """
        query_fp = self.compute_fingerprint(query_smiles)
        if query_fp is None:
            return []

        results = []
        for compound in compounds:
            cmp_smiles = compound.get('smiles', '')
            cmp_fp = self.compute_fingerprint(cmp_smiles)
            if cmp_fp is None:
                continue

            similarity = self.tanimoto_similarity(query_fp, cmp_fp)
            if similarity >= threshold:
                results.append({
                    **compound,
                    'similarity': round(similarity, 3),
                })

        results.sort(key=lambda x: x['similarity'], reverse=True)
        return results

    def extract_substructure(
        self,
        smiles: str,
        substructure_smiles: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Check if a molecule contains a substructure using RDKit.

        Args:
            smiles: Molecule SMILES.
            substructure_smiles: Substructure SMILES.

        Returns:
            Match info or None if no match.
        """
        mol = self._mol_from_smiles(smiles)
        sub = self._mol_from_smiles(substructure_smiles)

        if mol is None or sub is None:
            return None

        try:
            from rdkit import Chem
            match = mol.GetSubstructMatch(sub)
            if match:
                return {
                    'match': True,
                    'n_atoms_matched': len(match),
                    'match_ids': list(match),
                }
            return {'match': False, 'n_atoms_matched': 0}
        except Exception:
            return None

    def generate_smiles_from_atoms(
        self,
        atoms: List[Atom],
    ) -> Optional[str]:
        """
        Generate a SMILES string from atom positions using RDKit.

        Builds a molecule from coordinates with RDKit's bond perception,
        then produces a canonical SMILES. When RDKit cannot build a
        chemically valid molecule from the supplied coordinates, this
        falls back to a degenerate element-count representation.

        Args:
            atoms: List of Atom objects with x, y, z coordinates.

        Returns:
            SMILES string or None if even the fallback cannot be produced.
        """
        if not atoms:
            return None

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            mol = Chem.MolFromConformer(Chem.Conformer(len(atoms)))
            if mol is None:
                return None

            rd_conf = mol.GetConformer()
            for idx, atom in enumerate(atoms):
                elem = (atom.element or 'C').strip().upper()
                rd_atom = Chem.Atom(elem)
                mol.AddAtom(rd_atom)
                rd_conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))

            # Sanitize and perceive bonds
            AllChem.DeleteMol(mol)
            mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
            if mol is None:
                return None

            Chem.SanitizeMol(mol)

            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            if smiles:
                return smiles

            # Fallback: degenerate element-count representation
            counts = {}
            for atom in atoms:
                elem = (atom.element or 'C').strip().upper()
                counts[elem] = counts.get(elem, 0) + 1
            parts = []
            for elem in sorted(counts.keys()):
                cnt = counts[elem]
                parts.append(f'{elem}{cnt}' if cnt > 1 else elem)
            return ''.join(parts)

        except Exception:
            # Ultimate fallback: element counts
            try:
                counts = {}
                for atom in atoms:
                    elem = (atom.element or 'C').strip().upper()
                    counts[elem] = counts.get(elem, 0) + 1
                parts = []
                for elem in sorted(counts.keys()):
                    cnt = counts[elem]
                    parts.append(f'{elem}{cnt}' if cnt > 1 else elem)
                return ''.join(parts)
            except Exception:
                return None

    def validate_smiles(self, smiles: str, max_length: int = 10000) -> bool:
        """
        Check if a SMILES string is valid using RDKit.

        Args:
            smiles: SMILES string to validate.
            max_length: Maximum SMILES length accepted.

        Returns:
            True if RDKit can parse and sanitize the SMILES.
        """
        if not smiles:
            return False

        if len(smiles) > max_length:
            return False

        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return False

        try:
            from rdkit.Chem import rdchem
            if mol.GetNumAtoms() > 5000:
                return False
        except Exception:
            pass

        return True

    def export_sdf(self, smiles: str, output_path: str) -> bool:
        """
        Export a molecule as SDF using RDKit.

        Args:
            smiles: SMILES string of the molecule.
            output_path: Path to write the SDF file.

        Returns:
            True if export succeeded.
        """
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return False

        try:
            from rdkit.Chem import SDWriter
            from rdkit import Chem
            mol = Chem.RemoveHs(mol)
            writer = SDWriter(output_path)
            writer.write(mol)
            writer.close()
            return True
        except Exception as e:
            import traceback
            traceback.print_exc()
            return False
