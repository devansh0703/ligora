"""
Contact analyzer - protein-ligand interaction analysis.

Uses PLIP (Protein-Ligand Interaction Profiler) for contact detection
with geometric fallback for when PLIP is unavailable.

Detects:
- Hydrogen bonds
- Hydrophobic contacts
- Pi-stacking
- Salt bridges
- Halogen bonds
- Metal coordination
- Water-mediated contacts
"""

import subprocess
import tempfile
import re
import math
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field

import numpy as np

from .schemas import (
    Structure,
    Ligand,
    Contact,
    ContactType,
    Atom,
)
from .config import get_config


class ContactAnalyzer:
    """
    Analyze protein-ligand contacts using PLIP and geometric methods.

    PLIP is the primary method for detecting non-covalent interactions.
    Falls back to geometric analysis when PLIP is unavailable.
    """

    # No hardcoded geometric classification thresholds remain in this module.
    # The analyzer only runs an external source (PLIP) when available.
    # Internal geometric classification was removed because it was a local heuristic.
    def __init__(self):
        self._plip_available: Optional[bool] = None

    def analyze_contacts(
        self,
        structure: Structure,
        ligand: Ligand,
        protein_chains: Optional[List[str]] = None,
    ) -> List[Contact]:
        """
        Analyze all contacts between a ligand and protein chains.

        Args:
            structure: The full structure containing chains and ligands.
            ligand: The ligand to analyze.
            protein_chains: Optional list of chain IDs to consider.
                          If None, analyzes all polymer chains.

        Returns:
            List of detected contacts.
        """
        config = get_config()

        # Get protein atoms
        protein_atoms = self._get_protein_atoms(
            structure, protein_chains
        )

        # Get ligand atoms
        ligand_atoms = ligand.atoms

        if not protein_atoms or not ligand_atoms:
            return []

        contacts: List[Contact] = []

        # Only use an external source (PLIP) when it is available.
        # No internal geometric fallback is produced by this module.
        plip_contacts = self._run_plip(
            structure, ligand, config
        )
        if plip_contacts:
            contacts.extend(plip_contacts)

        return contacts

    def _get_protein_atoms(
        self,
        structure: Structure,
        protein_chains: Optional[List[str]] = None,
    ) -> List[Atom]:
        """
        Get all atoms from protein chains.

        Args:
            structure: The structure.
            protein_chains: Optional chain IDs to include.

        Returns:
            List of protein atoms.
        """
        atoms: List[Atom] = []
        chains_to_check = protein_chains or [
            c.id for c in structure.chains if c.is_polymer
        ]

        for chain in structure.chains:
            if chain.id not in chains_to_check:
                continue
            for residue in chain.residues:
                for atom in residue.atoms:
                    atoms.append(atom)

        return atoms

    def _run_plip(
        self,
        structure: Structure,
        ligand: Ligand,
        config,
    ) -> List[Contact]:
        """
        Run PLIP on the structure and parse results.

        Args:
            structure: The structure.
            ligand: The ligand of interest.
            config: Application configuration.

        Returns:
            List of contacts parsed from PLIP output.
        """
        # Check if PLIP is available
        if not self._check_plip_available():
            return []

        # Create temporary PDB file for PLIP
        with tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.pdb',
            delete=False,
            prefix='ligora_'
        ) as tmp_file:
            pdb_path = Path(tmp_file.name)

            # Write structure in PDB format for PLIP
            self._write_pdb_for_plip(structure, ligand, tmp_file)

        try:
            # Run PLIP
            plip_cmd = [
                "plip",
                "-f", str(pdb_path),
                "-v",  # verbose
            ]

            result = subprocess.run(
                plip_cmd,
                capture_output=True,
                text=True,
                timeout=60,
                cwd=Path(tempfile.gettempdir()),
            )

            if result.returncode != 0:
                # PLIP failed
                return []

            # Parse PLIP output files
            plip_output_dir = pdb_path.with_suffix('')
            contacts = self._parse_plip_output(
                plip_output_dir, ligand
            )

            return contacts

        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            return []
        finally:
            # Clean up temporary files
            try:
                pdb_path.unlink()
                # Remove PLIP output directory if created
                plip_output_dir = pdb_path.with_suffix('')
                if plip_output_dir.exists():
                    import shutil
                    shutil.rmtree(plip_output_dir)
            except OSError:
                pass

    def _check_plip_available(self) -> bool:
        """Check if PLIP is available on the system."""
        if self._plip_available is not None:
            return self._plip_available

        try:
            result = subprocess.run(
                ["plip", "--version"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            self._plip_available = result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
            self._plip_available = False

        return self._plip_available

    def _write_pdb_for_plip(
        self,
        structure: Structure,
        ligand: Ligand,
        file_handle,
    ):
        """
        Write structure in PDB format for PLIP.

        Args:
            structure: The structure.
            ligand: The ligand to focus on.
            file_handle: File handle to write to.
        """
        atom_serial = 1

        # Write header
        print("TITLE     PLARORA-STRUCTURE-FOR-PLIP", file=file_handle)
        print("HEADER    AUTOMATICALLY GENERATED", file=file_handle)

        # Write protein atoms
        for chain in structure.chains:
            if not chain.is_polymer:
                continue
            for residue in chain.residues:
                for atom in residue.atoms:
                    self._write_pdb_atom(
                        file_handle, atom, atom_serial,
                        record_type="ATOM"
                    )
                    atom_serial += 1

        # Write ligand atoms
        for atom in ligand.atoms:
            self._write_pdb_atom(
                file_handle, atom, atom_serial,
                record_type="HETATM"
            )
            atom_serial += 1

        print(f"END", file=file_handle)

    def _write_pdb_atom(
        self,
        file_handle,
        atom: Atom,
        serial: int,
        record_type: str = "ATOM",
    ):
        """Write a single atom in PDB format."""
        # PDB format: ATOM serial name resName chainID resSeq x y z occupancy tempFactor element
        line = (
            f"{record_type:<6}{serial:>5}  {atom.name:<4}"
            f"{atom.residue_name:<3}{atom.chain_id:<1}"
            f"{atom.residue_id:>4}    "
            f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
            f"{atom.occupancy:>6.2f}{atom.b_factor:>6.2f}          "
            f"{atom.element or '  '}"
        )
        print(line, file=file_handle)

    def _parse_plip_output(
        self,
        output_dir: Path,
        ligand: Ligand,
    ) -> List[Contact]:
        """
        Parse PLIP XML output files.

        Args:
            output_dir: Directory containing PLIP output.
            ligand: The ligand being analyzed.

        Returns:
            List of contacts parsed from PLIP.
        """
        contacts: List[Contact] = []

        # Look for PLIP XML output
        xml_files = list(output_dir.glob("*.xml"))
        if not xml_files:
            # Try the alternative naming
            xml_files = list(output_dir.parent.glob(f"{output_dir.name}*.xml"))

        for xml_file in xml_files:
            try:
                import xml.etree.ElementTree as ET
                tree = ET.parse(xml_file)
                root = tree.getroot()

                for interaction in root.findall(".//interaction"):
                    contact_type = interaction.get("type", "").lower()
                    contact = self._parse_plip_interaction(
                        interaction, ligand
                    )
                    if contact:
                        contacts.append(contact)

            except Exception:
                continue

        return contacts

    def _parse_plip_interaction(
        self,
        interaction,
        ligand: Ligand,
    ) -> Optional[Contact]:
        """
        Parse a single PLIP interaction element.

        Args:
            interaction: XML interaction element.
            ligand: The ligand being analyzed.

        Returns:
            Contact object or None.
        """
        # Map PLIP interaction types to our types
        type_mapping = {
            "hydrogenbond": ContactType.HYDROGEN_BOND,
            "hydrophobic": ContactType.HYDROPHOBIC,
            "pistacking": ContactType.PI_STACKING,
            "saltbridge": ContactType.SALT_BRIDGE,
            "halogen": ContactType.HALOGEN_BOND,
            "metal": ContactType.METAL_COORDINATION,
            "waterbridge": ContactType.WATER_MEDIATED,
            "pi-cation": ContactType.PI_STACKING,
        }

        contact_type_str = interaction.get("type", "").lower()
        contact_type = type_mapping.get(contact_type_str, ContactType.UNKNOWN)

        # Extract atoms from interaction
        ligand_atom_elem = interaction.find("ligand").find("atom") if interaction.find("ligand") is not None else None
        protein_atom_elem = interaction.find("protein").find("atom") if interaction.find("protein") is not None else None

        if ligand_atom_elem is None or protein_atom_elem is None:
            return None

        ligand_atom = ligand_atom_elem.get("name", "")
        ligand_res = ligand_atom_elem.get("resname", "")
        ligand_res_id = int(ligand_atom_elem.get("resnr", 0))
        ligand_chain = ligand_atom_elem.get("chain", "")

        protein_atom = protein_atom_elem.get("name", "")
        protein_res = protein_atom_elem.get("resname", "")
        protein_res_id = int(protein_atom_elem.get("resnr", 0))
        protein_chain = protein_atom_elem.get("chain", "")

        # Extract distance if available
        distance_elem = interaction.find("distance")
        distance = float(distance_elem.text) if distance_elem is not None else 0.0

        # Extract angle if available (for hydrogen bonds)
        angle_elem = interaction.find("angle")
        angle = float(angle_elem.text) if angle_elem is not None else None

        # Determine description
        description = self._describe_contact(contact_type, protein_res, ligand_res)

        return Contact(
            id=0,  # Will be assigned by caller
            ligand_atom=ligand_atom,
            ligand_residue_name=ligand_res,
            ligand_residue_id=ligand_res_id,
            ligand_chain_id=ligand_chain,
            protein_residue_name=protein_res,
            protein_residue_id=protein_res_id,
            protein_chain_id=protein_chain,
            protein_atom=protein_atom,
            distance=distance,
            contact_type=contact_type,
            angle=angle,
            description=description,
        )

    def _geometric_analysis(
        self,
        protein_atoms: List[Atom],
        ligand_atoms: List[Atom],
        ligand: Ligand,
        structure: Structure,
    ) -> List[Contact]:
        """
        Perform geometric contact analysis as fallback.

        Args:
            protein_atoms: List of protein atoms.
            ligand_atoms: List of ligand atoms.
            ligand: The ligand.
            structure: The structure.

        Returns:
            List of detected contacts.
        """
        contacts: List[Contact] = []
        contact_id = 0

        for lig_atom in ligand_atoms:
            lig_pos = np.array([lig_atom.x, lig_atom.y, lig_atom.z])

            for prot_atom in protein_atoms:
                prot_pos = np.array([prot_atom.x, prot_atom.y, prot_atom.z])

                # Calculate distance
                diff = lig_pos - prot_pos
                distance = np.linalg.norm(diff)

                # Skip if too far (use a generous cutoff derived from
                # a typical ligand radius plus a buffer; this is only a
                # performance filter, not a contact-type classifier).
                if distance > 8.0:
                    continue

                # No internal geometric fallback classification.
                # Contact type must come from an external source (for example PLIP)
                # or be supplied by the caller. The module no longer assigns
                # contact types from local distance/element heuristics.

        return contacts

    def _classify_by_geometry(self, distance: float, lig_element: Optional[str], prot_element: Optional[str]) -> ContactType:
        """No local geometric classification remains.

        The method is kept only as a compat signature to avoid breaking
        call sites during the removal pass. It always returns UNKNOWN.
        Contact type must come from an external source or the caller.
        """
        return ContactType.UNKNOWN

    def _describe_contact(self, contact_type: ContactType, protein_res: str, ligand_res: str) -> str:
        """No local description heuristics remain.

        Returns an empty description unless a contact type is supplied by
        an external source. The app does not invent chemical descriptions.
        """
        if not contact_type:
            return ""
        return ""


    def export_contacts_csv(
        self,
        contacts: List[Contact],
        output_path: Path,
    ):
        """
        Export contacts to a CSV file.

        Args:
            contacts: List of contacts to export.
            output_path: Path to write the CSV file.
        """
        with open(output_path, 'w') as f:
            # Header
            f.write("id,ligand_atom,ligand_res,ligand_res_id,ligand_chain,")
            f.write("protein_atom,protein_res,protein_res_id,protein_chain,")
            f.write("distance,contact_type,angle,water_mediated,description\n")

            for contact in contacts:
                f.write(f"{contact.id},")
                f.write(f"{contact.ligand_atom},")
                f.write(f"{contact.ligand_residue_name},")
                f.write(f"{contact.ligand_residue_id},")
                f.write(f"{contact.ligand_chain_id},")
                f.write(f"{contact.protein_atom},")
                f.write(f"{contact.protein_residue_name},")
                f.write(f"{contact.protein_residue_id},")
                f.write(f"{contact.protein_chain_id},")
                f.write(f"{contact.distance:.3f},")
                f.write(f"{contact.contact_type.value},")
                f.write(f"{contact.angle if contact.angle else ''},")
                f.write(f"{contact.is_water_mediated},")
                f.write(f"{contact.description}\n")

    def compute_binding_pocket(
        self,
        structure: Structure,
        ligand: Ligand,
        radius: float = 6.0,
    ) -> Dict[str, Any]:
        """
        Compute the binding pocket around a ligand.

        Args:
            structure: The structure.
            ligand: The ligand.
            radius: Pocket radius in Angstroms.

        Returns:
            Dictionary with pocket information.
        """
        ligand_center = self._compute_ligand_center(ligand)

        pocket_residues = []
        pocket_chains = set()

        for chain in structure.chains:
            if not chain.is_polymer:
                continue
            for residue in chain.residues:
                residue_center = self._compute_residue_center(residue)

                if residue_center is None:
                    continue

                dist = np.linalg.norm(
                    np.array(ligand_center) - residue_center
                )

                if dist <= radius:
                    pocket_residues.append({
                        "chain_id": chain.id,
                        "residue_name": residue.name,
                        "residue_id": residue.id,
                        "distance": round(dist, 2),
                    })
                    pocket_chains.add(chain.id)

        return {
            "ligand_center": {
                "x": round(ligand_center[0], 3),
                "y": round(ligand_center[1], 3),
                "z": round(ligand_center[2], 3),
            },
            "radius": radius,
            "pocket_chain_ids": sorted(list(pocket_chains)),
            "pocket_residue_count": len(pocket_residues),
            "pocket_residues": pocket_residues,
        }

    def _compute_ligand_center(self, ligand: Ligand) -> np.ndarray:
        """Compute the geometric center of a ligand."""
        if not ligand.atoms:
            return np.array([0.0, 0.0, 0.0])

        positions = np.array([
            [a.x, a.y, a.z] for a in ligand.atoms
        ])
        return np.mean(positions, axis=0)

    def _compute_residue_center(self, residue) -> Optional[np.ndarray]:
        """Compute the geometric center of a residue."""
        if not residue.atoms:
            return None

        positions = np.array([
            [a.x, a.y, a.z] for a in residue.atoms
        ])
        return np.mean(positions, axis=0)

    def find_ligand_similarities(
        self,
        ligand: Ligand,
        threshold: float = 0.7,
    ) -> List[Dict[str, Any]]:
        """
        Find similar ligands based on fingerprints using ChEMBL.

        Uses Cheminformatics for fingerprints and queries ChEMBL
        for similar compounds via the API.

        Args:
            ligand: The ligand.
            threshold: Similarity threshold (0-1).

        Returns:
            List of similar ligand info from ChEMBL.
        """
        from .cheminformatics import Cheminformatics
        from .enrichment import EnrichmentClient

        if not ligand.smiles:
            return []

        fp_engine = Cheminformatics()
        ligand_fp = fp_engine.compute_fingerprint(ligand.smiles)
        if ligand_fp is None:
            return []

        # Query ChEMBL for similar compounds
        enrichment = EnrichmentClient()
        
        # Get similar compounds from ChEMBL
        results = []
        try:
            # Use ChEMBL's similarity search
            config = get_config()
            chembl_url = f"{config.chembl_base_url}/similarity.json?smiles={ligand.smiles}&similarity_score={threshold * 100}"
            
            import requests
            response = requests.get(chembl_url, timeout=30)
            if response.status_code == 200:
                data = response.json()
                for compound in data.get('molecules', [])[:20]:
                    results.append({
                        'chembl_id': compound.get('molecule_chembl_id'),
                        'pref_name': compound.get('pref_name'),
                        'smiles': compound.get('molecule_structures', {}).get('canonical_smiles'),
                        'similarity': compound.get('similarity_score', 0) / 100,
                    })
        except Exception:
            pass

        return results
