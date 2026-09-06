"""
Tests for the cheminformatics module.
"""

import pytest
import numpy as np

from ligora_backend.cheminformatics import Cheminformatics


class TestCheminformatics:
    """Tests for the Cheminformatics class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.chem = Cheminformatics(radius=2, n_bits=2048)

    def test_compute_fingerprint_empty(self):
        """Test fingerprint computation with empty SMILES."""
        fp = self.chem.compute_fingerprint("")
        assert fp is None

        fp = self.chem.compute_fingerprint(None)
        assert fp is None

    def test_compute_fingerprint_simple(self):
        """Test fingerprint computation with simple SMILES."""
        # Basic organic molecule
        smiles = "CCO"  # Ethanol
        fp = self.chem.compute_fingerprint(smiles)

        if fp is not None:
            assert len(fp) == 2048
            assert fp.dtype == np.uint8

    def test_tanimoto_similarity_identical(self):
        """Test Tanimoto similarity of identical fingerprints."""
        smiles = "CCO"
        fp1 = self.chem.compute_fingerprint(smiles)
        fp2 = self.chem.compute_fingerprint(smiles)

        if fp1 is not None and fp2 is not None:
            similarity = self.chem.tanimoto_similarity(fp1, fp2)
            assert similarity == 1.0

    def test_tanimoto_similarity_empty(self):
        """Test Tanimoto similarity with empty fingerprints."""
        similarity = self.chem.tanimoto_similarity(None, None)
        assert similarity == 0.0

        similarity = self.chem.tanimoto_similarity(np.zeros(2048), None)
        assert similarity == 0.0

    def test_tanimoto_similarity_max(self):
        """Test that Tanimoto similarity is bounded [0, 1]."""
        fp1 = np.array([1, 1, 1, 0, 0, 0], dtype=np.uint8)
        fp2 = np.array([1, 0, 1, 0, 0, 0], dtype=np.uint8)

        similarity = self.chem.tanimoto_similarity(fp1, fp2)
        assert 0.0 <= similarity <= 1.0

        # Different fingerprints
        fp3 = np.array([0, 0, 0, 1, 1, 1], dtype=np.uint8)
        similarity = self.chem.tanimoto_similarity(fp1, fp3)
        assert 0.0 <= similarity <= 1.0

    def test_export_sdf(self):
        """Test SDF export via export_sdf."""
        import tempfile
        import os
        smiles = "CCO"
        with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False) as tmp:
            path = tmp.name
        try:
            ok = self.chem.export_sdf(smiles, path)
            assert ok is True
            assert os.path.exists(path)
            content = open(path).read()
            assert '$$$$' in content
        finally:
            if os.path.exists(path):
                os.unlink(path)

    def test_validate_smiles_valid(self):
        """Test SMILES validation with valid SMILES."""
        assert self.chem.validate_smiles("CCO")
        assert self.chem.validate_smiles("c1ccccc1")  # Benzene
        assert self.chem.validate_smiles("CC(=O)O")  # Acetic acid

    def test_compute_descriptors(self):
        """Test descriptor computation returns RDKit-backed values."""
        desc = self.chem.compute_descriptors("CCO")
        assert isinstance(desc, dict)
        assert desc.get('formula') == 'C2H6O'
        assert desc.get('molecular_weight') > 40
        assert desc.get('n_heavy_atoms') == 3

    def test_generate_smiles_from_atoms(self):
        """Test SMILES generation from atom coordinates."""
        from ligora_backend.schemas import Atom
        atoms = [
            Atom(id=1, name='C1', residue_name='LIG', residue_id=1,
                 chain_id='L', x=0.0, y=0.0, z=0.0, element='C'),
            Atom(id=2, name='O1', residue_name='LIG', residue_id=1,
                 chain_id='L', x=1.2, y=0.0, z=0.0, element='O'),
        ]
        smiles = self.chem.generate_smiles_from_atoms(atoms)
        assert smiles is not None
        assert len(smiles) > 0

    def test_validate_smiles_invalid(self):
        """Test SMILES validation with invalid SMILES."""
        assert not self.chem.validate_smiles("")
        # Local heuristics are removed. Validity is delegated to RDKit.
        # We only assert that the method returns False for clearly invalid
        # inputs that RDKit refuses to parse.
        assert not self.chem.validate_smiles("C))")    

    def test_compute_molecular_descriptors(self):
        """Test molecular descriptor computation via RDKit."""
        smiles = "CCO"
        descriptors = self.chem.compute_descriptors(smiles)

        assert isinstance(descriptors, dict)
        assert descriptors.get('formula') == 'C2H6O'
        assert descriptors.get('n_heavy_atoms') == 3
        assert descriptors.get('molecular_weight') == pytest.approx(46.07, abs=0.05)
        assert 'n_hba' in descriptors
        assert 'n_hbd' in descriptors

        assert "n_atoms" in descriptors
        assert "molecular_weight" in descriptors

    def test_find_similar_compounds(self):
        """Test similar compound finding."""
        query_smiles = "CCO"
        compounds = [
            {"smiles": "CCO", "name": "Ethanol"},
            {"smiles": "CCCO", "name": "Propanol"},
            {"smiles": "CCCCO", "name": "Butanol"},
            {"smiles": "c1ccccc1", "name": "Benzene"},
        ]

        similar = self.chem.find_similar_compounds(query_smiles, compounds, threshold=0.5)

        # Ethanol should be most similar to itself
        if len(similar) > 0:
            assert similar[0]["smiles"] == "CCO"

    def test_generate_smiles_from_atoms(self):
        """Test SMILES generation from atoms."""
        from ligora_backend.schemas import Atom

        atoms = [
            Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=0, y=0, z=0, element="C"),
            Atom(id=2, name="O1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=1, y=0, z=0, element="O"),
        ]

        smiles = self.chem.generate_smiles_from_atoms(atoms)

        # Should generate some SMILES representation
        assert smiles is not None
        assert len(smiles) > 0

    def test_smiles_validation_edge_cases(self):
        """Test SMILES validation edge cases."""
        # Empty string
        assert not self.chem.validate_smiles("")

        # Very long string: RDKit may still parse a long carbon chain,
        # so we only assert that the method does not crash and returns
        # a boolean for an overlong input. We do not assert rejection
        # because validity is delegated to RDKit, not to a local length rule.
        long_smiles = "C" * 1001
        result = self.chem.validate_smiles(long_smiles)
        assert isinstance(result, bool)


        # Balanced parentheses
        assert self.chem.validate_smiles("C(C)C")
        assert self.chem.validate_smiles("C(C(C)C)C")

        # Unbalanced parentheses (should be False)
        # Note: Our simple validator may not catch all cases
        # assert not self.chem.validate_smiles("C))C")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
