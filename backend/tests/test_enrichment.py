"""
Tests for the enrichment client module.
"""

import pytest
from unittest.mock import Mock, patch
from pathlib import Path

from ligora_backend.enrichment import EnrichmentClient


class TestEnrichmentClient:
    """Tests for the EnrichmentClient class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.client = EnrichmentClient()

    def test_get_structure_metadata_mock(self):
        """Test structure metadata with mocked response."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'title': 'Test Protein Structure',
            'structure_determination_method': {
                'method': 'X-RAY DIFFRACTION',
                'resolution': 2.0,
            },
            'deposited': '2023-01-01',
            'released': '2023-06-01',
            'macromolecule_type': 'Protein',
            'assembly': {'name': 'Monomer'},
            'source_organism': {'scientific_name': 'Homo sapiens'},
        }

        with patch('ligora_backend.enrichment.requests.Session.get', return_value=mock_response):
            metadata = self.client.get_structure_metadata('1ABC')

        assert metadata is not None
        assert metadata['title'] == 'Test Protein Structure'
        assert metadata['experiment_type'] == 'X-RAY DIFFRACTION'
        assert metadata['resolution'] == 2.0

    def test_get_pubchem_info_mock(self):
        """Test PubChem info with mocked response."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'IdentifierList': {
                'CID': [12345]
            }
        }

        mock_props_response = Mock()
        mock_props_response.status_code = 200
        mock_props_response.json.return_value = {
            'PropertyTable': {
                'Properties': [
                    {
                        'Title': 'Test Compound',
                        'CanonicalSMILES': 'CCO',
                        'InChIKey': 'TESTKEY123',
                        'MolecularWeight': 46.07,
                        'XLogP': -0.3,
                        'TSPA': 20.2,
                    }
                ]
            }
        }

        with patch('ligora_backend.enrichment.requests.Session.get') as mock_get:
            # First call for search
            mock_get.side_effect = [mock_response, mock_props_response]

            info = self.client.get_pubchem_info(compound_name="Test Compound")

        assert info is not None
        assert info['cid'] == 12345
        assert info['title'] == 'Test Compound'
        assert info['smiles'] == 'CCO'

    def test_get_chembl_info_mock(self):
        """Test ChEMBL info with mocked response."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'molecule': {
                'pref_name': 'Test Drug',
                'molecule_properties': {
                    'canonical_smiles': 'CCO',
                    'full_mwt': 46.07,
                    'alogp': -0.3,
                    'hba': 1,
                    'hbd': 1,
                }
            }
        }

        with patch('ligora_backend.enrichment.requests.Session.get', return_value=mock_response):
            info = self.client._get_chembl_by_id('CHEMBL123')

        assert info is not None
        assert info['chembl_id'] == 'CHEMBL123'
        assert info['pref_name'] == 'Test Drug'
        assert info['smiles'] == 'CCO'

    def test_get_pdbbind_info_mock(self):
        """Test PDBBind info with mocked response."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'affinity_value': 50.0,
            'affinity_type': 'Kd',
        }

        with patch('ligora_backend.enrichment.requests.Session.get', return_value=mock_response):
            info = self.client.get_pdbbind_info('1ABC', 'LIG')

        assert info is not None
        assert info['pdb_id'] == '1ABC'
        assert info['ligand'] == 'LIG'
        assert info['affinity'] == 50.0

    def test_enrich_ligand_complete(self):
        """Test complete ligand enrichment."""
        ligand = {
            'name': 'Test Compound',
            'residue_name': 'TCOM',
            'formula': 'C2H6O',
        }

        # Mock all external calls
        mock_pubchem = {
            'cid': 12345,
            'title': 'Test Compound',
            'smiles': 'CCO',
            'inchi_key': 'TESTKEY',
            'molecular_weight': 46.07,
            'xlogp': -0.3,
            'tpsa': 20.2,
        }

        mock_chembl = {
            'chembl_id': 'CHEMBL123',
            'pref_name': 'Test Drug',
            'smiles': 'CCO',
        }

        mock_pdbbind = {
            'pdb_id': '1ABC',
            'ligand': 'TCOM',
            'affinity': 50.0,
        }

        with patch.object(self.client, 'get_pubchem_info', return_value=mock_pubchem):
            with patch.object(self.client, 'get_chembl_info', return_value=mock_chembl):
                with patch.object(self.client, 'get_pdbbind_info', return_value=mock_pdbbind):
                    enriched = self.client.enrich_ligand(ligand, pdb_id='1ABC')

        assert enriched['pubchem_cid'] == 12345
        assert enriched['smiles'] == 'CCO'
        assert enriched['chembl_id'] == 'CHEMBL123'
        assert enriched['pdbbind_affinity'] == 50.0
        assert 'evidence' in enriched
        assert len(enriched['evidence']) == 3

    def test_cache_functionality(self):
        """Test caching behavior."""
        # Set a cache entry
        self.client._set_cache('test_key', {'data': 'value'})

        # Should be retrievable
        assert 'test_key' in self.client._cache

        # Should be valid
        assert self.client._is_cache_valid('test_key')

    def test_clear_cache(self):
        """Test cache clearing."""
        self.client._set_cache('key1', 'value1')
        self.client._set_cache('key2', 'value2')

        assert len(self.client._cache) == 2

        self.client.clear_cache()

        assert len(self.client._cache) == 0

    def test_health_check_mock(self):
        """Test health check with mocked responses."""
        mock_response_ok = Mock()
        mock_response_ok.status_code = 200

        with patch('ligora_backend.enrichment.requests.Session.get') as mock_get:
            mock_get.return_value = mock_response_ok

            health = self.client.health_check()

        assert health['rscb'] == 'ok'
        assert health['pubchem'] == 'ok'
        assert health['chembl'] == 'ok'
        assert health['pdbbind'] == 'ok'

    def test_get_all_evidence(self):
        """Test evidence collection."""
        ligand = {
            'pubchem_cid': 12345,
            'chembl_id': 'CHEMBL123',
            'pdbbind_affinity': 50.0,
        }

        evidence = self.client.get_all_evidence(ligand, pdb_id='1ABC')

        assert len(evidence) == 3

        # Check PubChem evidence
        pubchem_ev = next(e for e in evidence if e.source == 'pubchem')
        assert pubchem_ev.field == 'compound_identity'
        assert pubchem_ev.value['cid'] == 12345

        # Check ChEMBL evidence
        chembl_ev = next(e for e in evidence if e.source == 'chembl')
        assert chembl_ev.field == 'bioactivity'
        assert chembl_ev.value['chembl_id'] == 'CHEMBL123'

        # Check PDBBind evidence
        pdbbind_ev = next(e for e in evidence if e.source == 'pdbbind')
        assert pdbbind_ev.field == 'binding_affinity'
        assert pdbbind_ev.value['affinity'] == 50.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
