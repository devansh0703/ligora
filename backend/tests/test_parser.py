"""
Tests for the structure parser module.

Uses real PDB/mmCIF content and a real RCSB fetch. No service mocks: when
the network is unavailable the live test fails, which is the honest result.
"""

import pytest
import requests

from ligora_backend.parser import StructureParser


PDB_CONTENT = """HEADER    TEST STRUCTURE
TITLE     TEST PROTEIN
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  1.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  1.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  1.00           C
ATOM      4  O   ALA A   1       4.000   2.000   3.000  1.00  1.00           O
ATOM      5  CB  ALA A   1       2.000   3.000   3.000  1.00  1.00           C
HETATM    6  C   LIG A 101       5.000   5.000   5.000  1.00  1.00           C
HETATM    7  N   LIG A 101       6.000   5.000   5.000  1.00  1.00           N
HETATM    8  O   LIG A 101       7.000   5.000   5.000  1.00  1.00           O
HETATM    9  O   HOH A 201       9.000   9.000   9.000  1.00  1.00           O
HETATM   10  O   HOH A 202       9.000   9.000   9.000  1.00  1.00           O
END
"""


@pytest.fixture(scope="module")
def cif_3w85_content() -> str:
    """Real mmCIF content for PDB entry 3W85 (fetched once)."""
    response = requests.get(
        "https://files.rcsb.org/download/3W85.cif", timeout=60)
    assert response.status_code == 200
    return response.text


class TestStructureParser:
    """Tests for the StructureParser class."""

    def setup_method(self):
        self.parser = StructureParser()

    def test_parse_empty_mmcif(self):
        structure = self.parser.parse_mmcif("# Empty file")
        assert structure.id == "unknown"
        assert structure.source == "local"
        assert len(structure.chains) == 0
        assert len(structure.ligands) == 0

    def test_parse_simple_mmcif(self):
        # All rows share auth_asym_id A and are ATOM records, so they form
        # one polymer chain even without an _entity category. The file ends
        # without a trailing '#' to cover loop termination at EOF.
        content = "data_test\n" \
            "loop_\n" \
            "_atom_site.group_PDB\n" \
            "_atom_site.auth_asym_id\n" \
            "_atom_site.auth_comp_id\n" \
            "_atom_site.auth_seq_id\n" \
            "_atom_site.auth_atom_id\n" \
            "_atom_site.type_symbol\n" \
            "_atom_site.Cartn_x\n" \
            "_atom_site.Cartn_y\n" \
            "_atom_site.Cartn_z\n" \
            "ATOM A ALA 1 N N 10.0 10.0 10.0\n" \
            "ATOM A ALA 1 CA C 11.0 10.0 10.0\n" \
            "ATOM A ALA 1 C C 12.0 10.0 10.0\n"
        structure = self.parser.parse_mmcif(content, source_id="test")
        assert structure.id == "test"
        assert structure.source == "local"
        assert len(structure.chains) == 1
        assert structure.chains[0].id == "A"
        assert len(structure.chains[0].residues) == 1
        assert len(structure.chains[0].residues[0].atoms) == 3
        assert [a.name for a in structure.chains[0].residues[0].atoms] == \
            ["N", "CA", "C"]

    def test_parse_pdb_format(self):
        structure = self.parser.parse_pdb(PDB_CONTENT)
        assert structure.id == "unknown"
        assert structure.file_format == "pdb"
        assert structure.title == "TEST PROTEIN"

        chain = next(c for c in structure.chains if c.id == "A")
        assert chain.is_polymer
        assert len(chain.residues) == 1
        assert len(chain.residues[0].atoms) == 5

        ligand = next(lig for lig in structure.ligands if lig.residue_name == "LIG")
        assert ligand.atom_count == 3
        assert [a.element for a in ligand.atoms] == ["C", "N", "O"]

    def test_parse_pdb_water_grouping(self):
        """Waters become non-polymer instances; elements come from the file."""
        structure = self.parser.parse_pdb(PDB_CONTENT)
        waters = [lig for lig in structure.ligands if lig.residue_name == "HOH"]
        assert len(waters) == 1  # grouped per (component, chain)
        assert waters[0].atom_count == 2
        assert all(a.element == "O" for a in waters[0].atoms)

    def test_parse_pdb_keeps_last_residue(self):
        """Regression: the final residue of a PDB chain must not be lost."""
        content = PDB_CONTENT.replace("ALA A   1", "ALA A   1") + \
            "ATOM     11  N   GLY A   5      11.000  12.000  13.000" \
            "  1.00  1.00           N\n"
        structure = self.parser.parse_pdb(content)
        chain = next(c for c in structure.chains if c.id == "A")
        residue_ids = {r.id for r in chain.residues}
        assert residue_ids == {1, 5}

    def test_parse_pdb_altloc_first_kept(self):
        # Alternate locations: only altLoc 'A' (or blank) is kept.
        # Column layout: name in cols 13-16, altLoc char in col 17.
        # PDB columns (1-based): atom name 13-16, altLoc 17. With a
        # 1-char atom name, 'B' lands in the altLoc column (17).
        alt_record = ("HETATM    8  O B LIG A 101       7.000   5.000"
                      "   5.000  1.00  1.00           O")
        assert alt_record[13] == "O" and alt_record[15] == "B"
        content = PDB_CONTENT.replace(
            "HETATM    8  O   LIG A 101", alt_record)
        structure = self.parser.parse_pdb(content)
        ligand = next(lig for lig in structure.ligands if lig.residue_name == "LIG")
        names = [a.name for a in ligand.atoms]
        assert names == ["C", "N", "O"], names  # altLoc B atom dropped

    def test_fetch_from_rcsb_real(self):
        """Live fetch of a real small structure (crambin, 1CRN)."""
        structure = self.parser.fetch_from_rcsb("1CRN", format="mmCIF")
        assert structure.id == "1CRN"
        assert structure.source == "rcsb"
        assert any(c.is_polymer for c in structure.chains)
        assert any(r.name == "ALA" for c in structure.chains
                   for r in c.residues)

    def test_parse_mmcif_3w85_real(self, cif_3w85_content):
        """Real structure: chains, ligand instances, elements from the file."""
        structure = self.parser.parse_mmcif(
            cif_3w85_content, source_id="3W85", source="rcsb")
        assert structure.id == "3W85"
        assert structure.title and structure.title != "Unknown structure"
        assert structure.resolution is not None
        assert structure.experiment_type

        polymer_chains = [c for c in structure.chains if c.is_polymer]
        assert len(polymer_chains) >= 1

        w85 = [lig for lig in structure.ligands if lig.residue_name == "W85"]
        assert len(w85) == 2  # one instance per binding chain
        assert all(len(lig.atoms) == 23 for lig in w85)
        assert all(a.element for lig in w85 for a in lig.atoms)

    def test_complete_missing_elements_from_ccd_real(
            self, cif_3w85_content):
        """CCD-backed element completion uses the live component files."""
        # Waters in 3W85 carry elements from the file; use a PDB fixture
        # with element columns stripped to exercise the CCD path.
        stripped = "\n".join(
            line[:76].ljust(78) if line.startswith(("ATOM", "HETATM"))
            else line
            for line in PDB_CONTENT.splitlines()) + "\n"
        structure = self.parser.parse_pdb(stripped)
        ligand = next(lig for lig in structure.ligands if lig.residue_name == "LIG")
        assert all(a.element is None for a in ligand.atoms)
        updated = self.parser.complete_missing_elements(structure)
        # LIG is not a real CCD component: nothing may be invented for it.
        assert all(a.element is None for a in ligand.atoms)
        assert isinstance(updated, int)

    def test_complete_missing_elements_from_ccd_real_component(
            self, cif_3w85_content):
        """A real CCD component (GOL) gets authoritative elements."""
        content = """HETATM    1  C1  GOL A 101       0.000   0.000   0.000  1.00  1.00
HETATM    2  C2  GOL A 101       1.500   0.000   0.000  1.00  1.00
HETATM    3  C3  GOL A 101       3.000   0.000   0.000  1.00  1.00
HETATM    4  O1  GOL A 101       4.000   1.200   0.000  1.00  1.00
HETATM    5  O2  GOL A 101       1.900  -1.300   0.500  1.00  1.00
HETATM    6  O3  GOL A 101       3.600  -1.400  -0.300  1.00  1.00
END
"""
        structure = self.parser.parse_pdb(content)
        gol = next(lig for lig in structure.ligands if lig.residue_name == "GOL")
        assert all(a.element is None for a in gol.atoms)
        self.parser.complete_missing_elements(structure)
        assert [a.element for a in gol.atoms] == [
            "C", "C", "C", "O", "O", "O"]


if __name__ == "__main__":
    pytest.main([__file__])
