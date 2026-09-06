"""
Cheminformatics operations backed by RDKit and Open Babel.

All chemistry is delegated to established open-source toolkits:
- RDKit: SMILES parsing/validation, Morgan fingerprints, Tanimoto
  similarity, descriptors, 2D coordinate generation, SDF/mol blocks,
  bond perception from 3D coordinates (rdDetermineBonds).
- Open Babel (obabel CLI): 3D coordinate generation and force-field
  geometry cleanup when requested.

No chemical data (atomic numbers, weights, bond orders, thresholds) is
embedded in this module. When a molecule cannot be built from the supplied
coordinates, the failure is reported - nothing is invented.
"""

from typing import Optional, List, Dict, Any, Tuple
from pathlib import Path

import numpy as np

from .schemas import Atom


class Cheminformatics:
    """Cheminformatics operations backed by RDKit (and Open Babel for 3D)."""

    def __init__(self, radius: Optional[int] = None,
                 n_bits: Optional[int] = None):
        from .config import get_config
        config = get_config()
        self.radius = radius if radius is not None else config.fingerprint_radius
        self.n_bits = n_bits if n_bits is not None else config.fingerprint_nbits

    # ------------------------------------------------------------------
    # Molecule construction (RDKit)
    # ------------------------------------------------------------------

    def _mol_from_smiles(self, smiles: str):
        if not smiles:
            return None
        try:
            from rdkit import Chem
            return Chem.MolFromSmiles(smiles, sanitize=True)
        except Exception:
            return None

    def _mol_from_atoms(self, atoms: List[Atom]):
        """
        Build an RDKit molecule from atoms with 3D coordinates.

        Bond orders come from real sources only, tried in order:
        1. RDKit bond perception (rdDetermineBonds) on the deposited
           coordinates when hydrogens are present.
        2. The CCD component's own `_chem_comp_bond` records, mapped onto
           the deposited coordinates by atom name (hydrogen positions
           superposed from the CCD's ideal frame via Kabsch alignment).
        Components without CCD coverage fail honestly - no bonds are
        guessed from element counts or distances.
        """
        if not atoms:
            return None
        if any(not a.element for a in atoms):
            # Element data missing: refuse to guess. Callers should complete
            # elements from the CCD before calling.
            return None

        # 1) Direct perception on the deposited coordinates.
        mol = self._build_rdmol(atoms)
        if mol is not None:
            try:
                from rdkit.Chem import rdDetermineBonds
                charged = self._mol_charge(atoms)
                rdDetermineBonds.DetermineBonds(mol, charge=charged)
                return mol
            except Exception:
                pass

        # 2) CCD bond records mapped by atom name.
        return self._mol_from_ccd_bonds(atoms)

    def _build_rdmol(self, atoms: List[Atom]):
        """Build an RDKit molecule (no bonds) from atoms + coordinates."""
        try:
            from rdkit import Chem
            mol = Chem.RWMol()
            conf = Chem.Conformer(len(atoms))
            for idx, atom in enumerate(atoms):
                mol.AddAtom(Chem.Atom(atom.element.strip().upper()))
                conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))
            mol.AddConformer(conf, assignId=True)
            return mol.GetMol()
        except Exception:
            return None

    def _mol_charge(self, atoms: List[Atom]) -> int:
        """Total formal charge from the CCD's per-atom charges (0 w/o CCD)."""
        from .ligand import ccd_atom_charges
        comp_ids = {a.residue_name for a in atoms}
        if len(comp_ids) != 1:
            return 0
        charges = ccd_atom_charges(next(iter(comp_ids)))
        return sum(charges.get(a.name, 0) for a in atoms)

    def _mol_from_ccd_bonds(self, atoms: List[Atom]):
        """Molecule from the CCD's own bond records, mapped by atom name.

        Hydrogens absent from the deposited coordinates get their positions
        from the CCD's ideal frame, Kabsch-superposed onto the deposited
        heavy atoms so the molecule stays geometrically consistent.
        """
        import numpy as np
        try:
            from rdkit import Chem
            from .ligand import (ccd_bond_list, ccd_ideal_coordinates,
                                 ccd_atom_charges)

            comp_ids = {a.residue_name for a in atoms}
            if len(comp_ids) != 1:
                return None
            comp_id = next(iter(comp_ids))
            bonds = ccd_bond_list(comp_id)
            if not bonds:
                return None

            name_to_idx = {a.name: i for i, a in enumerate(atoms)}

            # Kabsch: superpose CCD ideal heavy atoms onto deposited heavy
            # atoms to place missing hydrogens in the deposited frame.
            ideal = ccd_ideal_coordinates(comp_id)
            h_positions: Dict[str, tuple] = {}
            if ideal:
                common = [n for n in name_to_idx
                          if n in ideal and ideal[n][0].upper() != "H"]
                missing_h = [n for n, v in ideal.items()
                             if v[0].upper() == "H" and n not in name_to_idx]
                if common and missing_h:
                    P = np.array([[ideal[n][1], ideal[n][2], ideal[n][3]]
                                  for n in common])
                    Q = np.array([[atoms[name_to_idx[n]].x,
                                   atoms[name_to_idx[n]].y,
                                   atoms[name_to_idx[n]].z] for n in common])
                    p_cent = P.mean(axis=0)
                    q_cent = Q.mean(axis=0)
                    H = (P - p_cent).T @ (Q - q_cent)
                    U, S, Vt = np.linalg.svd(H)
                    d = np.sign(np.linalg.det(Vt.T @ U.T))
                    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
                    for n in missing_h:
                        v = np.array([ideal[n][1], ideal[n][2], ideal[n][3]])
                        placed = R @ (v - p_cent) + q_cent
                        h_positions[n] = tuple(placed)

            charges = ccd_atom_charges(comp_id)

            mol = Chem.RWMol()
            conf = Chem.Conformer(len(atoms) + len(h_positions))
            for idx, atom in enumerate(atoms):
                a = Chem.Atom(atom.element.strip().upper())
                a.SetFormalCharge(charges.get(atom.name, 0))
                mol.AddAtom(a)
                conf.SetAtomPosition(idx, (atom.x, atom.y, atom.z))
            # Append missing hydrogens (CCD frame, superposed).
            h_idx: Dict[str, int] = {}
            for name, (x, y, z) in h_positions.items():
                h_idx[name] = mol.GetNumAtoms()
                a = Chem.Atom("H")
                a.SetFormalCharge(charges.get(name, 0))
                mol.AddAtom(a)
                conf.SetAtomPosition(mol.GetNumAtoms() - 1, (x, y, z))
            mol.AddConformer(conf, assignId=True)

            order_map = {"SING": Chem.BondType.SINGLE,
                         "DOUB": Chem.BondType.DOUBLE,
                         "TRIP": Chem.BondType.TRIPLE,
                         "QUAD": Chem.BondType.QUADRUPLE}
            for bond in bonds:
                i = name_to_idx.get(bond["atom_id_1"],
                                    h_idx.get(bond["atom_id_1"]))
                j = name_to_idx.get(bond["atom_id_2"],
                                    h_idx.get(bond["atom_id_2"]))
                if i is None or j is None:
                    continue
                bt = order_map.get(bond["value_order"], Chem.BondType.SINGLE)
                mol.AddBond(i, j, bt)

            out = mol.GetMol()
            from rdkit.Chem import rdmolops, SanitizeFlags
            rdmolops.SanitizeMol(out, SanitizeFlags.SANITIZE_ALL ^
                                 SanitizeFlags.SANITIZE_PROPERTIES)
            return out
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Fingerprints and similarity
    # ------------------------------------------------------------------

    def compute_fingerprint(self, smiles: str) -> Optional[np.ndarray]:
        """Morgan fingerprint via RDKit; None when the SMILES is invalid."""
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
            on_bits = fp.GetOnBits()
            for i in on_bits:
                arr[i] = 1
            return arr
        except Exception:
            return None

    def tanimoto_similarity(self, fp1: Optional[np.ndarray],
                            fp2: Optional[np.ndarray]) -> float:
        """Tanimoto similarity via RDKit DataStructs; numpy fallback."""
        if fp1 is None or fp2 is None:
            return 0.0
        try:
            from rdkit import DataStructs as rd_ds
            from rdkit.DataStructs import ExplicitBitVect

            bv1 = ExplicitBitVect(self.n_bits)
            bv2 = ExplicitBitVect(self.n_bits)
            for i in np.nonzero(fp1)[0]:
                bv1.SetBit(int(i))
            for i in np.nonzero(fp2)[0]:
                bv2.SetBit(int(i))
            return float(rd_ds.TanimotoSimilarity(bv1, bv2))
        except Exception:
            intersection = float(np.sum(fp1 * fp2))
            sum1 = float(np.sum(fp1))
            sum2 = float(np.sum(fp2))
            union = sum1 + sum2 - intersection
            if union == 0:
                return 0.0
            return intersection / union

    # ------------------------------------------------------------------
    # SMILES / SDF / mol blocks
    # ------------------------------------------------------------------

    def validate_smiles(self, smiles: str, max_length: int = 10000) -> bool:
        """A SMILES is valid when RDKit can parse and sanitize it."""
        if not smiles or len(smiles) > max_length:
            return False
        return self._mol_from_smiles(smiles) is not None

    def smiles_to_sdf(self, smiles: str, molecule_name: str = "MOL") -> str:
        """SDF block from SMILES via RDKit; empty SDF when invalid."""
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return "$$$$\n"
        try:
            from rdkit.Chem import rdmolfiles
            mol.SetProp("_Name", molecule_name)
            return rdmolfiles.MolToMolBlock(mol) + "$$$$\n"
        except Exception:
            return "$$$$\n"

    def mol_block_from_atoms(self, atoms: List[Atom]) -> Optional[str]:
        """Mol block from 3D coordinates with RDKit bond perception."""
        mol = self._mol_from_atoms(atoms)
        if mol is None:
            return None
        try:
            from rdkit.Chem import rdmolfiles
            return rdmolfiles.MolToMolBlock(mol)
        except Exception:
            return None

    def sdf_export_from_atoms(self, atoms: List[Atom]) -> Optional[str]:
        """
        SDF from atom coordinates via RDKit bond perception.

        Returns None when RDKit cannot build a valid molecule (e.g. missing
        element data, unperceivable bonds). No bondless pseudo-SDF is
        fabricated, because such files misrepresent the chemistry.
        """
        block = self.mol_block_from_atoms(atoms)
        if block is None:
            return None
        return block + "$$$$\n"

    def atoms_to_sdf(self, ligand) -> str:
        """Export a ligand object as SDF; raises when chemistry is unknown."""
        sdf = self.sdf_export_from_atoms(ligand.atoms)
        if sdf is None:
            raise ValueError(
                f"Cannot build a valid molecule for ligand "
                f"{ligand.residue_name}: element data missing or bond "
                f"perception failed. Complete element data (e.g. from the "
                f"CCD) before exporting SDF."
            )
        return sdf

    def generate_2d_coords(
        self, smiles: str) -> Optional[List[Tuple[float, float, float]]]:
        """2D coordinates for a SMILES via RDKit."""
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return None
        try:
            from rdkit.Chem import AllChem
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            return [(conf.GetAtomPosition(i).x,
                     conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z)
                    for i in range(mol.GetNumAtoms())]
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Descriptors (RDKit)
    # ------------------------------------------------------------------

    def compute_descriptors(self, smiles: str) -> Dict[str, Any]:
        """Molecular descriptors via RDKit; empty when SMILES invalid."""
        mol = self._mol_from_smiles(smiles)
        if mol is None:
            return {}
        return self.compute_descriptors_from_mol(mol)

    def compute_descriptors_from_mol(self, mol) -> Dict[str, Any]:
        if mol is None:
            return {}
        try:
            from rdkit.Chem import Descriptors
            from rdkit.Chem import rdMolDescriptors

            return {
                "n_heavy_atoms": int(
                    rdMolDescriptors.CalcNumHeavyAtoms(mol)),
                "n_hba": int(rdMolDescriptors.CalcNumHBA(mol)),
                "n_hbd": int(rdMolDescriptors.CalcNumHBD(mol)),
                "n_atoms": mol.GetNumAtoms(),
                "logp": round(float(Descriptors.MolLogP(mol)), 3),
                "tpsa": round(float(Descriptors.TPSA(mol)), 3),
                "molecular_weight": round(
                    float(Descriptors.MolWt(mol)), 3),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "atom_counts": {
                    a.GetSymbol(): sum(
                        1 for x in mol.GetAtoms()
                        if x.GetSymbol() == a.GetSymbol())
                    for a in mol.GetAtoms()
                },
            }
        except Exception:
            return {}

    def generate_smiles_from_atoms(self, atoms: List[Atom]) -> Optional[str]:
        """SMILES from 3D coordinates; None when RDKit cannot build one."""
        if not atoms:
            return None
        mol = self._mol_from_atoms(atoms)
        if mol is None:
            return None
        try:
            from rdkit import Chem
            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            return smiles or None
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Similarity search over supplied compounds
    # ------------------------------------------------------------------

    def find_similar_compounds(
        self,
        query_smiles: str,
        compounds: List[Dict[str, Any]],
        threshold: float = 0.7,
    ) -> List[Dict[str, Any]]:
        """Rank compounds by Tanimoto similarity to the query."""
        query_fp = self.compute_fingerprint(query_smiles)
        if query_fp is None:
            return []
        results: List[Dict[str, Any]] = []
        for compound in compounds:
            cmp_fp = self.compute_fingerprint(compound.get("smiles", ""))
            if cmp_fp is None:
                continue
            similarity = self.tanimoto_similarity(query_fp, cmp_fp)
            if similarity >= threshold:
                entry = dict(compound)
                entry["similarity"] = round(similarity, 3)
                results.append(entry)
        results.sort(key=lambda x: x["similarity"], reverse=True)
        return results

    def extract_substructure(self, smiles: str,
                             substructure_smiles: str
                             ) -> Optional[Dict[str, Any]]:
        """Substructure match via RDKit."""
        mol = self._mol_from_smiles(smiles)
        sub = self._mol_from_smiles(substructure_smiles)
        if mol is None or sub is None:
            return None
        try:
            match = mol.GetSubstructMatch(sub)
            if match:
                return {"match": True, "n_atoms_matched": len(match),
                        "match_ids": list(match)}
            return {"match": False, "n_atoms_matched": 0}
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Open Babel (3D generation, geometry cleanup)
    # ------------------------------------------------------------------

    def _obabel_path(self) -> Optional[str]:
        """Resolve the Open Babel executable."""
        import shutil
        from .config import get_config
        config = get_config()
        if config.obabel_executable:
            return config.obabel_executable
        return shutil.which("obabel")

    def is_obabel_available(self) -> bool:
        return self._obabel_path() is not None

    def generate_3d_coords(self, smiles: str) -> Optional[List[Atom]]:
        """
        Generate 3D coordinates from SMILES using Open Babel's gen3d.

        Returns Atom records (element + coordinates) or None when Open Babel
        is unavailable or the SMILES is invalid. No coordinates are made up
        locally.
        """
        obabel = self._obabel_path()
        if not obabel or not smiles:
            return None
        import subprocess
        import tempfile
        try:
            with tempfile.TemporaryDirectory(prefix="ligora_obabel_") as td:
                inp = Path(td) / "in.smi"
                outp = Path(td) / "out.xyz"
                inp.write_text(smiles + "\n")
                result = subprocess.run(
                    [obabel, str(inp), "-O", str(outp), "--gen3d"],
                    capture_output=True, text=True, timeout=120,
                )
                if result.returncode != 0 or not outp.exists():
                    return None
                return self._parse_xyz(outp.read_text())
        except (OSError, subprocess.TimeoutExpired):
            return None

    @staticmethod
    def _parse_xyz(xyz_text: str) -> Optional[List[Atom]]:
        """Parse an XYZ file (Open Babel output) into Atom records."""
        lines = xyz_text.strip().splitlines()
        if len(lines) < 2:
            return None
        try:
            n = int(lines[0].split()[0])
        except (ValueError, IndexError):
            return None
        atoms: List[Atom] = []
        for i, line in enumerate(lines[2:2 + n]):
            parts = line.split()
            if len(parts) < 4:
                return None
            atoms.append(Atom(
                id=i + 1,
                name=parts[0],
                residue_name="LIG",
                residue_id=1,
                chain_id="L",
                x=float(parts[1]), y=float(parts[2]), z=float(parts[3]),
                element=parts[0].upper(),
            ))
        return atoms if atoms else None

    def cleanup_geometry(self, atoms: List[Atom]) -> Dict[str, Any]:
        """
        Geometry cleanup (force-field minimization) via Open Babel.

        Input is written as XYZ with perceived connectivity (Open Babel
        perceives bonds from distance/element data), minimized with the
        configured force field, and returned with the minimized coordinates
        plus the energy reported by Open Babel. When Open Babel is not
        available, the result reports that honestly.
        """
        obabel = self._obabel_path()
        if not obabel:
            return {"status": "unavailable",
                    "message": "Open Babel (obabel) is not installed"}
        if any(not a.element for a in atoms):
            return {"status": "unavailable",
                    "message": "Element data missing; complete elements "
                               "from the CCD before geometry cleanup"}

        from .config import get_config
        import subprocess
        import tempfile
        config = get_config()

        try:
            with tempfile.TemporaryDirectory(prefix="ligora_geom_") as td:
                inp = Path(td) / "in.xyz"
                outp = Path(td) / "out.xyz"
                inp.write_text(self._write_xyz(atoms))
                result = subprocess.run(
                    [obabel, str(inp), "-O", str(outp),
                     "--minimize",
                     "--ff", config.geometry_cleanup_force_field,
                     "--steps", str(config.geometry_cleanup_steps),
                     "--sd", str(config.geometry_cleanup_tolerance),
                     "--log"],
                    capture_output=True, text=True, timeout=300,
                )
                log = result.stdout + result.stderr
                if not outp.exists():
                    return {"status": "failed", "message": log[-500:]}
                new_atoms = self._parse_xyz(outp.read_text())
                if not new_atoms or len(new_atoms) != len(atoms):
                    return {"status": "failed", "message": "minimization "
                            "changed atom count; refusing result"}

                rmsd = float(np.sqrt(np.mean([
                    (a.x - b.x) ** 2 + (a.y - b.y) ** 2 +
                    (a.z - b.z) ** 2
                    for a, b in zip(atoms, new_atoms)])))

                # Open Babel reports the final energy in its log.
                energy = None
                for line in log.splitlines():
                    if "TOTAL ENERGY" in line.upper():
                        for tok in line.split():
                            try:
                                energy = float(tok)
                            except ValueError:
                                continue
                        break

                return {
                    "status": "ok",
                    "force_field": config.geometry_cleanup_force_field,
                    "steps": config.geometry_cleanup_steps,
                    "atoms": new_atoms,
                    "rmsd": round(rmsd, 4),
                    "energy_kcal_mol": energy,
                }
        except (OSError, subprocess.TimeoutExpired) as e:
            return {"status": "failed", "message": str(e)}

    @staticmethod
    def _write_xyz(atoms: List[Atom]) -> str:
        lines = [str(len(atoms)), "Ligora geometry cleanup input"]
        for a in atoms:
            lines.append(
                f"{a.element.strip().upper()} {a.x:.6f} {a.y:.6f} {a.z:.6f}")
        return "\n".join(lines) + "\n"
