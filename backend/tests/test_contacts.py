"""
Tests for the contact analysis module.

Contact classification is PLIP's job - the module never assigns interaction
types from local geometry. Unit tests use constructed content; the
end-to-end test runs real PLIP on a real RCSB structure.
"""

import pytest
import numpy as np
from pathlib import Path

from ligora_backend.contacts import ContactAnalyzer
from ligora_backend.schemas import (
    Structure, Chain, Residue, Atom, Ligand, Contact, ContactType
)


def make_structure_with_ligand() -> tuple:
    structure = Structure(id="test", source="local")
    chain = Chain(id="A", name="Protein", is_polymer=True)

    residues_data = [
        ("ALA", 1, [(1, "N", "N", 1.0, 1.0, 1.0),
                    (2, "CA", "C", 2.0, 1.0, 1.0),
                    (3, "C", "C", 3.0, 1.0, 1.0),
                    (4, "O", "O", 4.0, 1.0, 1.0)]),
        ("GLY", 2, [(5, "N", "N", 5.0, 1.0, 1.0),
                    (6, "CA", "C", 6.0, 1.0, 1.0),
                    (7, "C", "C", 7.0, 1.0, 1.0),
                    (8, "O", "O", 8.0, 1.0, 1.0)]),
    ]

    for res_name, res_id, atoms_data in residues_data:
        residue = Residue(id=res_id, name=res_name, chain_id="A",
                          residue_number=res_id)
        for atom_id, atom_name, element, x, y, z in atoms_data:
            residue.atoms.append(Atom(
                id=atom_id, name=atom_name, residue_name=res_name,
                residue_id=res_id, chain_id="A",
                x=x, y=y, z=z, element=element,
                b_factor=30.0, occupancy=1.0))
        chain.residues.append(residue)

    structure.chains.append(chain)

    ligand = Ligand(id="L1", name="Test Ligand", residue_name="LIG",
                    atom_count=3)
    ligand.atoms = [
        Atom(id=1, name="C1", residue_name="LIG", residue_id=101,
             chain_id="A", x=4.0, y=1.0, z=1.0, element="C"),
        Atom(id=2, name="O1", residue_name="LIG", residue_id=101,
             chain_id="A", x=5.0, y=1.0, z=1.0, element="O"),
        Atom(id=3, name="N1", residue_name="LIG", residue_id=101,
             chain_id="A", x=6.0, y=1.0, z=1.0, element="N"),
    ]
    structure.ligands.append(ligand)
    return structure, ligand


class TestContactAnalyzer:
    """Tests for the ContactAnalyzer class."""

    def setup_method(self):
        self.analyzer = ContactAnalyzer()

    def test_plip_availability_reported(self):
        # Availability must be a real check, whatever the environment has.
        assert isinstance(self.analyzer.is_plip_available(), bool)

    def test_compute_binding_pocket(self):
        structure, ligand = make_structure_with_ligand()
        pocket = self.analyzer.compute_binding_pocket(structure, ligand,
                                                      radius=6.0)
        assert "ligand_center" in pocket
        assert "radius" in pocket
        assert "pocket_chain_ids" in pocket
        assert "pocket_residue_count" in pocket
        assert pocket["pocket_residue_count"] >= 1
        center = pocket["ligand_center"]
        assert {"x", "y", "z"} <= set(center)

    def test_compute_binding_pocket_radius_selects(self):
        structure, ligand = make_structure_with_ligand()
        near = self.analyzer.compute_binding_pocket(structure, ligand,
                                                    radius=3.0)
        far = self.analyzer.compute_binding_pocket(structure, ligand,
                                                   radius=0.5)
        assert near["pocket_residue_count"] >= far["pocket_residue_count"]

    def test_export_contacts_csv(self):
        contacts = [
            Contact(
                id=0,
                ligand_atom="C1",
                ligand_residue_name="LIG",
                ligand_residue_id=101,
                ligand_chain_id="L",
                protein_residue_name="ALA",
                protein_residue_id=1,
                protein_chain_id="A",
                protein_atom="O",
                distance=2.5,
                contact_type=ContactType.HYDROGEN_BOND,
                angle=150.0,
                description="donor/acceptor",
            ),
        ]
        output_path = Path("/tmp/test_contacts.csv")
        self.analyzer.export_contacts_csv(contacts, output_path)
        content = output_path.read_text()
        assert "ligand_atom" in content
        assert "protein_atom" in content
        assert "distance" in content
        assert "hydrogen_bond" in content
        output_path.unlink()

    def test_compute_ligand_center(self):
        ligand = Ligand(id="L1", name="Test", residue_name="LIG",
                        atom_count=3)
        ligand.atoms = [
            Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=1.0, y=2.0, z=3.0, element="C"),
            Atom(id=2, name="C2", residue_name="LIG", residue_id=1,
                 chain_id="L", x=3.0, y=4.0, z=5.0, element="C"),
            Atom(id=3, name="C3", residue_name="LIG", residue_id=1,
                 chain_id="L", x=5.0, y=6.0, z=7.0, element="C"),
        ]
        center = self.analyzer._compute_ligand_center(ligand)
        assert np.allclose(center, [3.0, 4.0, 5.0])

    def test_analyze_without_plip_honest(self):
        """No contacts may be fabricated when PLIP is missing."""
        structure, ligand = make_structure_with_ligand()
        if self.analyzer.is_plip_available():
            pytest.skip("PLIP is installed; the unavailable path is moot")
        contacts, available = self.analyzer.analyze_contacts(
            structure, ligand)
        assert contacts == []
        assert available is False

    def test_analyze_real_plip_on_real_structure(self):
        """End-to-end: real RCSB structure, real PLIP run, real contacts."""
        if not self.analyzer.is_plip_available():
            pytest.skip("PLIP is not installed")

        import requests
        from ligora_backend.parser import StructureParser

        response = requests.get(
            "https://files.rcsb.org/download/3W85.cif", timeout=60)
        assert response.status_code == 200
        structure = StructureParser().parse_mmcif(
            response.text, source_id="3W85", source="rcsb")

        ligand = next(
            lig for lig in structure.ligands
            if lig.residue_name == "W85" and len(lig.atoms) == 23)
        for a in ligand.atoms:
            assert a.element  # elements required by the PDB writer

        contacts, available = self.analyzer.analyze_contacts(
            structure, ligand)
        assert available is True
        # W85 makes real hydrogen bonds / hydrophobic contacts in 3W85.
        assert len(contacts) > 0
        types = {c.contact_type for c in contacts}
        known = {t for t in types if t != ContactType.UNKNOWN}
        assert known  # PLIP classified at least one interaction
        for c in contacts:
            assert 0 < c.distance < 10
            assert c.protein_residue_name
            assert c.description


if __name__ == "__main__":
    pytest.main([__file__])
