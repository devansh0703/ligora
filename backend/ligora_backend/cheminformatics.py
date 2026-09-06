"""
Cheminformatics module for ligand similarity, fingerprints, and structure handling.

Provides:
- SMILES validation via RDKit
- Morgan/ECFP fingerprints via RDKit
- Tanimoto similarity using RDKit DataStructs
- SDF / mol block export via RDKit
- 2D coordinate generation via RDKit
- SMILES generation only when RDKit can produce a chemically valid one from
  supplied coordinates; otherwise no invented SMILES is produced
- Molecular weight, formula, logP, TPSA, HBA, HBD, heavy atom counts via RDKit

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
        self.radius = radius
        self.n_bits = n_bits

    # ------------------------------------------------------------------
    # Core chemistry operations (all RDKit)
    # ------------------------------------------------------------------

    def compute_fingerprint(self, smiles: str) -> Optional[np.ndarray]:
        if not smiles:
            return None
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return None
        try:
            from rdkit.Chem import rdFingerprintGenerator

            gen = rdFingerprintGenerator.GetMorganGenerator(
                radius=self.radius, fpSize=self.n_bits
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
        if fp1 is None or fp2 is None:
            return 0.0
        try:
            from rdkit import DataStructs as rd_ds

            return float(
                rd_ds.TanimotoSimilarity(
                    rd_ds.CreateNumpyBitVect(fp1),
                    rd_ds.CreateNumpyBitVect(fp2),
                )
            )
        except Exception:
            intersection = float(np.sum(fp1 * fp2))
            sum1 = float(np.sum(fp1))
            sum2 = float(np.sum(fp2))
            union = sum1 + sum2 - intersection
            if union == 0:
                return 0.0
            return intersection / union

    def smiles_to_sdf(self, smiles: str, molecule_name: str = "MOL") -> str:
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return "$$$$\n"
        try:
            from rdkit.Chem import rdmolfiles

            return rdmolfiles.MolToMolBlock(mol) + "$$$$\n"
        except Exception:
            return "$$$$\n"

    def mol_block_from_atoms(self, atoms: List[Atom]) -> Optional[str]:
        mol = self._mol_from_atoms(atoms)
        if mol is None:
            return None
        try:
            from rdkit.Chem import rdmolfiles

            return rdmolfiles.MolToMolBlock(mol)
        except Exception:
            return None

    def generate_2d_coords(self, smiles: str) -> Optional[List[Tuple[float, float, float]]]:
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return None
        try:
            from rdkit.Chem import AllChem

            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            return [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z) for i in range(mol.GetNumAtoms())]
        except Exception:
            return None

    def compute_descriptors(self, smiles: str) -> Dict[str, Any]:
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return {}
        try:
            from rdkit.Chem import Descriptors
            from rdkit.Chem import rdMolDescriptors

            return {
                "n_heavy_atoms": int(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
                "n_hba": int(rdMolDescriptors.CalcNumHBA(mol)),
                "n_hbd": int(rdMolDescriptors.CalcNumHBD(mol)),
                "n_atoms": mol.GetNumAtoms(),
                "logp": round(float(Descriptors.MolLogP(mol)), 3),
                "tpsa": round(float(Descriptors.TPSA(mol)), 3),
                "molecular_weight": round(float(Descriptors.MolWt(mol)), 3),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "atom_counts": {
                    a.GetSymbol(): sum(1 for x in mol.GetAtoms() if x.GetSymbol() == a.GetSymbol())
                    for a in mol.GetAtoms()
                },
            }
        except Exception:
            return {}

    def compute_descriptors_from_mol(self, mol) -> Dict[str, Any]:
        if mol is None:
            return {}
        try:
            from rdkit.Chem import Descriptors
            from rdkit.Chem import rdMolDescriptors

            return {
                "n_heavy_atoms": int(rdMolDescriptors.CalcNumHeavyAtoms(mol)),
                "n_hba": int(rdMolDescriptors.CalcNumHBA(mol)),
                "n_hbd": int(rdMolDescriptors.CalcNumHBD(mol)),
                "n_atoms": mol.GetNumAtoms(),
                "logp": round(float(Descriptors.MolLogP(mol)), 3),
                "tpsa": round(float(Descriptors.TPSA(mol)), 3),
                "molecular_weight": round(float(Descriptors.MolWt(mol)), 3),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "atom_counts": {
                    a.GetSymbol(): sum(1 for x in mol.GetAtoms() if x.GetSymbol() == a.GetSymbol())
                    for a in mol.GetAtoms()
                },
            }
        except Exception:
            return {}

    def generate_smiles_from_atoms(self, atoms: List[Atom]) -> Optional[str]:
        """
        Generate a SMILES only when RDKit can produce a chemically valid one
        from the supplied coordinates. If RDKit cannot, this returns None
        rather than inventing a degenerate element-count or any other
        fabricated SMILES.
        """
        if not atoms:
            return None
        mol = self._mol_from_atoms(atoms)
        if mol is None:
            return None
        try:
            from rdkit import Chem

            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            if smiles:
                return smiles
            return None
        except Exception:
            return None

    def validate_smiles(self, smiles: str, max_length: int = 10000) -> bool:
        if not smiles or len(smiles) > max_length:
            return False
        return self._mol_from_smiles(smiles) is not None

    def find_similar_compounds(
        self,
        query_smiles: str,
        compounds: List[Dict[str, Any]],
        threshold: float = 0.7,
    ) -> List[Dict[str, Any]]:
        query_fp = self.compute_fingerprint(query_smiles)
        if query_fp is None:
            return []
        results: List[Dict[str, Any]] = []
        for compound in compounds:
            cmp_smiles = compound.get("smiles", "")
            cmp_fp = self.compute_fingerprint(cmp_smiles)
            if cmp_fp is None:
                continue
            similarity = self.tanimoto_similarity(query_fp, cmp_fp)
            if similarity >= threshold:
                entry = dict(compound)
                entry["similarity"] = round(similarity, 3)
                results.append(entry)
        results.sort(key=lambda x: x["similarity"], reverse=True)
        return results

    def extract_substructure(
        self,
        smiles: str,
        substructure_smiles: str,
    ) -> Optional[Dict[str, Any]]:
        mol = self._mol_from_smiles(smiles)
        sub = self._mol_from_smiles(substructure_smiles)
        if mol is None or sub is None:
            return None
        try:
            from rdkit import Chem

            match = mol.GetSubstructMatch(sub)
            if match:
                return {"match": True, "n_atoms_matched": len(match), "match_ids": list(match)}
            return {"match": False, "n_atoms_matched": 0}
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Internal helpers (RDKit only)
    # ------------------------------------------------------------------

    def _mol_from_smiles(self, smiles: str):
        try:
            from rdkit import Chem

            return Chem.MolFromSmiles(smiles, sanitize=True)
        except Exception:
            return None

    def _mol_from_atoms(self, atoms: List[Atom]):
        if not atoms:
            return None
        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolops

            conf = Chem.Conformer(len(atoms))
            mol = Chem.RWMol()
            for idx, atom in enumerate(atoms):
                elem = (atom.element or "C").strip().upper()
                mol.AddAtom(Chem.Atom(elem))
                conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))
            mol.AddConformer(conf, assignId=True)

            ref = Chem.MolFromSmiles("C", sanitize=True)
            if ref is None:
                return None
            ref = Chem.RemoveHs(ref)
            candidate = Chem.RemoveHs(mol)
            candidate = rdmolops.AssignBondOrdersFromTemplate(ref, candidate)
            candidate = Chem.RemoveHs(candidate)
            return candidate
        except Exception:
            return None
