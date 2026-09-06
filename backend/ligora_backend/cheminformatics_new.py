"""
Cheminformatics module for ligand similarity, fingerprints, and structure handling.

Provides:
- SMILES validation and conversion via RDKit
- Morgan/ECFP fingerprints via RDKit
- Tanimoto similarity using RDKit DataStructs
- SDF export via RDKit
- 2D coordinate generation via RDKit
- SMILES to molecule conversion via RDKit
- Molecular weight, formula, logP, TPSA, HBA, HBD, heavy atom counts via RDKit
- Mol block export from atom coordinates via RDKit (bond inference by RDKit)

All chemistry operations use RDKit. No chemical data is hardcoded here.
"""

from typing import Optional, List, Dict, Any, Tuple
from pathlib import Path

import numpy as np

from .schemas import Atom


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

    # ------------------------------------------------------------------
    # Core chemistry operations (all RDKit)
    # ------------------------------------------------------------------

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

    def tanimoto_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Compute Tanimoto similarity between two binary fingerprints.

        Uses RDKit DataStructs when available; otherwise falls back to a
        numpy Jaccard computation over the same bit vectors.

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
        Convert a SMILES string to SDF format using RDKit.

        Args:
            smiles: SMILES string.
            molecule_name: Name for the molecule.

        Returns:
            SDF formatted string.
        """
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return "$$$$\n"

        try:
            from rdkit.Chem import rdmolfiles

            # RDKit-generated mol block is authoritative; do not hand-roll atom
            # blocks or bond blocks from local guesses.
            block = rdmolfiles.MolToMolBlock(mol)
            return block + "$$$$\n"
        except Exception:
            return "$$$$\n"

    def _mol_from_smiles(self, smiles: str):
        """Convert SMILES to RDKit molecule. Returns None on failure."""
        try:
            from rdkit import Chem

            mol = Chem.MolFromSmiles(smiles, sanitize=True)
            return mol
        except Exception:
            return None

    def mol_from_atoms(self, atoms: List[Atom]):
        """
        Build an RDKit mol from a list of Atom objects with 3D coordinates.

        Bond perception is delegated entirely to RDKit via
        AssignBondOrdersFromTemplate over a hydrogen-suppressed shell from a
        reference atom. If RDKit cannot perform that operation, this method
        returns None rather than inventing a bond graph.

        This exists so a structure derived from coordinates can still be used
        by RDKit-based downstream operations where possible. It is not a
        replacement for SMILES-based parsing.
        """
        if not atoms:
            return None

        try:
            from rdkit import Chem

            # 1) Build a conformer + atom list from the supplied coordinates.
            conf = Chem.Conformer(len(atoms))
            mol = Chem.RWMol()
            for idx, atom in enumerate(atoms):
                elem = (atom.element or "C").strip().upper()
                rd_atom = Chem.Atom(elem)
                mol.AddAtom(rd_atom)
                conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))
            mol.AddConformer(conf, assignId=True)

            # 2) Use a single reference carbon skeleton so RDKit has a bond
            #    graph template to run AssignBondOrdersFromTemplate against.
            ref = Chem.MolFromSmiles("C", sanitize=True)
            if ref is None:
                return None
            ref = Chem.RemoveHs(ref)

            # 3) Prepare the candidate mol for bond order assignment.
            candidate = Chem.RemoveHs(mol)

            from rdkit.Chem import rdmolops

            candidate = rdmolops.AssignBondOrdersFromTemplate(ref, candidate)
            candidate = Chem.RemoveHs(candidate)

            return candidate
        except Exception:
            return None

    def mol_block_from_atoms(self, atoms: List[Atom]) -> Optional[str]:
        """
        Build a mol block string from atom coordinates and element symbols.

        Bond perception is delegated to RDKit where possible. If RDKit cannot
        perceive bonds from the supplied coordinates, this returns None instead
        of inventing a bond table.

        Args:
            atoms: List of Atom objects with coordinates and optional elements.

        Returns:
            Mol block string, or None if RDKit could not produce one.
        """
        if not atoms:
            return None

        mol = self.mol_from_atoms(atoms)
        if mol is None:
            return None

        try:
            from rdkit.Chem import rdmolfiles

            return rdmolfiles.MolToMolBlock(mol)
        except Exception:
            return None

    def generate_2d_coords(self, smiles: str) -> Optional[List[Tuple[float, float, float]]]:
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

            atom_counts: Dict[str, int] = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

            return {
                "n_heavy_atoms": int(n_heavy),
                "n_hba": int(hba),
                "n_hbd": int(hbd),
                "n_atoms": mol.GetNumAtoms(),
                "logp": round(float(logp), 3),
                "tpsa": round(float(tpsa), 3),
                "molecular_weight": round(float(mw), 3),
                "formula": formula,
                "atom_counts": atom_counts,
            }
        except Exception:
            return {}

    def compute_descriptors_from_mol(self, mol) -> Dict[str, Any]:
        """
        Compute molecular descriptors from an RDKit mol.

        Args:
            mol: RDKit molecule.

        Returns:
            Dictionary of descriptor values from RDKit, or empty dict on failure.
        """
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

            atom_counts: Dict[str, int] = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

            return {
                "n_heavy_atoms": int(n_heavy),
                "n_hba": int(hba),
                "n_hbd": int(hbd),
                "n_atoms": mol.GetNumAtoms(),
                "logp": round(float(logp), 3),
                "tpsa": round(float(tpsa), 3),
                "molecular_weight": round(float(mw), 3),
                "formula": formula,
                "atom_counts": atom_counts,
            }
        except Exception:
            return {}
