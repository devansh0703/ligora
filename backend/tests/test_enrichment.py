"""
Tests for the enrichment client module.

All tests run against the live PubChem / ChEMBL / UniChem / RCSB services.
No mocks: if a service is down the test fails, which is the honest result.
"""

import pytest

from ligora_backend.enrichment import EnrichmentClient


class TestEnrichmentClient:
    """Tests for the EnrichmentClient class (live services)."""

    def setup_method(self):
        self.client = EnrichmentClient()

    def test_get_structure_metadata_real(self):
        metadata = self.client.get_structure_metadata('1CRN')
        assert metadata is not None
        assert metadata['pdb_id'] == '1CRN'
        assert metadata['title']
        assert metadata['experiment_type']

    def test_get_structure_metadata_invalid(self):
        assert self.client.get_structure_metadata('XXXX') is None

    def test_get_pubchem_info_real(self):
        """Aspirin is CID 2244 in PubChem."""
        info = self.client.get_pubchem_info(pubchem_cid=2244)
        assert info is not None
        assert info['cid'] == 2244
        assert info['InChIKey'] == 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        assert info['smiles']

    def test_get_pubchem_info_by_name_real(self):
        info = self.client.get_pubchem_info(compound_name='ethanol')
        assert info is not None
        assert info['cid']
        assert info['smiles'] == 'CCO'

    def test_get_pubchem_info_unknown(self):
        assert self.client.get_pubchem_info(
            compound_name='zzz-no-such-compound-zzz') is None

    def test_get_chembl_info_real(self):
        info = self.client._get_chembl_by_id('CHEMBL25')
        assert info is not None
        assert info['chembl_id'] == 'CHEMBL25'

    def test_get_chembl_via_unichem_real(self):
        """Ethanol's InChIKey maps to CHEMBL545 via UniChem."""
        info = self.client.get_chembl_info(
            inchi_key='LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
        assert info is not None
        assert info['chembl_id'] == 'CHEMBL545'

    def test_get_chembl_activity_real(self):
        activities = self.client.get_chembl_activity('CHEMBL25', limit=3)
        assert isinstance(activities, list)
        if activities:
            assert 'standard_type' in activities[0]

    def test_cache_roundtrip(self):
        self.client._cache['k'] = {'v': 1}
        self.client._cache_expiry['k'] = float('inf')
        assert self.client._is_cache_valid('k')
        self.client.clear_cache()
        assert not self.client._is_cache_valid('k')

    def test_health_check_real(self):
        health = self.client.health_check()
        for service in ('rcsb', 'pubchem', 'chembl', 'unichem'):
            assert health.get(service) == 'ok'

    def test_get_all_evidence(self):
        ligand = {
            'pubchem_cid': 2244,
            'chembl_id': 'CHEMBL25',
            'residue_name': 'LIG',
        }
        evidence = self.client.get_all_evidence(ligand, pdb_id='1ABC')
        sources = {e.source for e in evidence}
        assert 'pubchem' in sources
        assert 'chembl' in sources
        for e in evidence:
            assert e.url

    def test_ccd_lookup_real(self):
        comp = self.client.get_compound_from_ccd('HOH')
        assert comp is not None
        assert comp['id'] == 'HOH'
        assert comp['formula'] == 'H2 O'


if __name__ == "__main__":
    pytest.main([__file__])
