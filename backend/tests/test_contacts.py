"""
Tests for the contact analysis module.
"""

import pytest
import numpy as np
from pathlib import Path

from ligora_backend.contacts import ContactAnalyzer
from ligora_backend.schemas import (
    Structure, Chain, Residue, Atom, Ligand, Contact, ContactType
)


class TestContactAnalyzer:
    """Tests for the ContactAnalyzer class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.analyzer = ContactAnalyzer()

    def create_test_structure(self):
        """Create a simple test structure."""
        structure = Structure(id="test", source="local")

        # Create a protein chain
        chain = Chain(id="A", name="Protein", is_polymer=True)

        # Add residues
        residues_data = [
            ("ALA", 1, [(1, "N", "N", 1.0, 1.0, 1.0),
                       (2, "CA", "C", 2.0, 1.0, 1.0),
                       (3, "C", "C", 3.0, 1.0, 1.0),
                       (4, "O", "O", 4.0, 1.0, 1.0)]),
            ("ALA", 2, [(5, "N", "N", 5.0, 1.0, 1.0),
                       (6, "CA", "C", 6.0, 1.0, 1.0),
                       (7, "C", "C", 7.0, 1.0, 1.0),
                       (8, "O", "O", 8.0, 1.0, 1.0)]),
            ("GLY", 3, [(9, "N", "N", 9.0, 1.0, 1.0),
                       (10, "CA", "C", 10.0, 1.0, 1.0),
                       (11, "C", "C", 11.0, 1.0, 1.0),
                       (12, "O", "O", 12.0, 1.0, 1.0)]),
        ]

        for res_name, res_id, atoms_data in residues_data:
            residue = Residue(id=res_id, name=res_name, chain_id="A", residue_number=res_id)
            for atom_id, atom_name, element, x, y, z in atoms_data:
                atom = Atom(
                    id=atom_id, name=atom_name, residue_name=res_name,
                    residue_id=res_id, chain_id="A",
                    x=x, y=y, z=z, element=element,
                    b_factor=30.0, occupancy=1.0
                )
                residue.atoms.append(atom)
            chain.residues.append(residue)

        structure.chains.append(chain)

        # Create a ligand near the protein
        ligand = Ligand(
            id="L1",
            name="Test Ligand",
            residue_name="LIG",
            formula="C6H6O",
            atom_count=3,
        )

        ligand.atoms = [
            Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=4.0, y=1.0, z=1.0, element="C"),
            Atom(id=2, name="O1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=5.0, y=1.0, z=1.0, element="O"),
            Atom(id=3, name="N1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=6.0, y=1.0, z=1.0, element="N"),
        ]

        structure.ligands.append(ligand)
        return structure, ligand

    def test_geometric_contacts_require_external_type(self):
        """Test that geometric analysis no longer assigns contact types.

        Contact type must come from an external source (for example PLIP)
        or be supplied by the caller. The module no longer assigns types
        from local distance/element heuristics.
        """
        structure, ligand = self.create_test_structure()

        protein_oxygen = Atom(
            id=100, name="O", residue_name="ALA", residue_id=1,
            chain_id="A", x=4.5, y=1.0, z=1.0, element="O"
        )
        structure.chains[0].residues[0].atoms.append(protein_oxygen)

        contacts = self.analyzer._geometric_analysis(
            structure.chains[0].residues[0].atoms,
            ligand.atoms,
            ligand,
            structure
        )

        # Any returned contact must not be locally classified.
        for contact in contacts:
            assert contact.contact_type == ContactType.UNKNOWN

    def test_geometric_hydrophobic_contact_is_not_locally_classified(self):
        """Test that hydrophobic contacts are not locally classified.

        Hydrophobic contact type must come from an external source.
        """
        structure, ligand = self.create_test_structure()

        protein_carbon = Atom(
            id=100, name="CB", residue_name="ALA", residue_id=1,
            chain_id="A", x=4.5, y=1.5, z=1.5, element="C"
        )
        structure.chains[0].residues[0].atoms.append(protein_carbon)

        contacts = self.analyzer._geometric_analysis(
            structure.chains[0].residues[0].atoms,
            ligand.atoms,
            ligand,
            structure
        )

        # No local hydrophobic classification.
        for contact in contacts:
            assert contact.contact_type == ContactType.UNKNOWN

    def test_compute_binding_pocket(self):
        """Test binding pocket computation."""
        structure, ligand = self.create_test_structure()

        pocket = self.analyzer.compute_binding_pocket(structure, ligand, radius=6.0)

        assert "ligand_center" in pocket
        assert "radius" in pocket
        assert "pocket_chain_ids" in pocket
        assert "pocket_residue_count" in pocket

        # Ligand center should be computed
        center = pocket["ligand_center"]
        assert "x" in center
        assert "y" in center
        assert "z" in center

    def test_export_contacts_csv(self):
        """Test CSV export of contacts."""
        contacts = [
            Contact(
                id=1,
                ligand_atom="C1",
                ligand_residue_name="LIG",
                ligand_residue_id=1,
                ligand_chain_id="L",
                protein_residue_name="ALA",
                protein_residue_id=1,
                protein_chain_id="A",
                protein_atom="O",
                distance=2.5,
                contact_type=ContactType.HYDROGEN_BOND,
                angle=150.0,
                description="H-bond between ALA and LIG",
            ),
            Contact(
                id=2,
                ligand_atom="C2",
                ligand_residue_name="LIG",
                ligand_residue_id=1,
                ligand_chain_id="L",
                protein_residue_name="GLY",
                protein_residue_id=2,
                protein_chain_id="A",
                protein_atom="CA",
                distance=3.5,
                contact_type=ContactType.HYDROPHOBIC,
                description="Hydrophobic contact",
            ),
        ]

        output_path = Path("/tmp/test_contacts.csv")
        self.analyzer.export_contacts_csv(contacts, output_path)

        # Verify file was created
        assert output_path.exists()

        # Read and verify content
        content = output_path.read_text()
        assert "ligand_atom" in content
        assert "protein_atom" in content
        assert "distance" in content

        # Cleanup
        output_path.unlink()

    def test_compute_ligand_center(self):
        """Test ligand center computation."""
        ligand = Ligand(id="L1", name="Test", residue_name="LIG", atom_count=3)
        ligand.atoms = [
            Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=1.0, y=2.0, z=3.0, element="C"),
            Atom(id=2, name="C2", residue_name="LIG", residue_id=1,
                 chain_id="L", x=3.0, y=4.0, z=5.0, element="C"),
            Atom(id=3, name="C3", residue_name="LIG", residue_id=1,
                 chain_id="L", x=5.0, y=6.0, z=7.0, element="C"),
        ]

        center = self.analyzer._compute_ligand_center(ligand)

        # Center should be average of positions
        assert np.allclose(center, [3.0, 4.0, 5.0])

    def test_classify_by_geometry_returns_unknown(self):
        """Test geometry-based classification returns UNKNOWN.

        The module no longer classifies contacts from local geometry.
        Contact type must come from an external source or the caller.
        """
        analyzer = self.analyzer

        # Hydrogen bond (short distance, N/O atoms)
        result = analyzer._classify_by_geometry(2.8, "N", "O")
        assert result == ContactType.UNKNOWN

        # Hydrophobic (carbon atoms, medium distance)
        result = analyzer._classify_by_geometry(3.8, "C", "C")
        assert result == ContactType.UNKNOWN

        # No contact (too far)
        result = analyzer._classify_by_geometry(6.0, "C", "C")
        assert result == ContactType.UNKNOWN

        # Van der Waals (medium distance C-C)
        result = analyzer._classify_by_geometry(3.0, "C", "C")
        assert result == ContactType.UNKNOWN


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
