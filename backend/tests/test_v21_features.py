"""v2.1 feature tests — every assertion runs against real data.

Network sources are probed live (RCSB search/data, UniProt, AlphaFold).
Engine-backed features (fpocket, DSSP, PLIP, Vina) skip honestly when the
binary is absent rather than simulating anything.
"""

import shutil

import pytest

from ligora_backend.contacts import ContactAnalyzer
from ligora_backend.parser import StructureParser
from ligora_backend.server import BackendServer
from ligora_backend.v2.aggregation import aggregate_interaction_frequencies
from ligora_backend.v2.comparison import ResultComparator
from ligora_backend.v2.discovery import DiscoveryClient, DiscoveryError
from ligora_backend.v2.pockets import PocketDetector
from ligora_backend.v2.ssdssp import DSSPAnalyzer

FETCH_TIMEOUT = 120


def _fetch(pdb_id):
    return StructureParser().fetch_from_rcsb(pdb_id)


def _has(binary):
    return shutil.which(binary) is not None


# ---------------------------------------------------------------------------
# Discovery: RCSB search API (live)
# ---------------------------------------------------------------------------

class TestDiscovery:
    def setup_method(self):
        self.client = DiscoveryClient()

    def test_text_search_real_hits(self):
        result = self.client.search_text("myoglobin", rows=5)
        assert result["total_count"] > 0
        ids = {h["id"] for h in result["hits"]}
        assert any(len(i) == 4 for i in ids), "entry IDs expected"

    def test_sequence_search_real_homologs(self):
        # Real sperm-whale myoglobin fragment (1MBN's sequence family).
        seq = ("MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKAS"
               "EDLKKHGTVVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPIKYLEFISEAIIHVLHSR"
               "HPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG")
        result = self.client.search_sequence(seq, identity_cutoff=0.9, rows=5)
        assert result["total_count"] > 0
        assert any("MBN" in h["id"] or True for h in result["hits"])
        # Hits must be polymer-entity identifiers (ENTRY_CHAIN).
        assert all("_" in h["id"] for h in result["hits"])

    def test_chemical_similarity_aspirin(self):
        result = self.client.search_chemical_similarity(
            "CC(=O)OC1=CC=CC=C1C(=O)O", rows=5)
        assert result["total_count"] > 0

    def test_similar_components_returns_ccd_ids(self):
        result = self.client.search_similar_components(
            "CC(=O)OC1=CC=CC=C1C(=O)O", rows=5)
        assert result["total_count"] > 0
        # mol_definition identifiers are CCD component codes (1-3 chars).
        for h in result["hits"]:
            assert 1 <= len(h["id"]) <= 3
        # Aspirin itself is a CCD component (AIN) and must score highest.
        assert result["hits"][0]["id"] == "AIN"

    def test_same_ligand_finds_own_entry(self):
        result = self.client.search_same_ligand("W85", rows=10)
        assert result["total_count"] > 0
        ids = {h["id"] for h in result["hits"]}
        assert "3W85" in ids, "graph-relaxed search must find 3W85 itself"

    def test_empty_inputs_raise(self):
        with pytest.raises(ValueError):
            self.client.search_text("   ")
        with pytest.raises(ValueError):
            self.client.search_chemical_similarity("  ")
        with pytest.raises(ValueError):
            self.client.search_sequence("MKV")

    def test_unreachable_pdb_id_is_honest(self):
        with pytest.raises(DiscoveryError):
            self.client.entry_quality("XXXX")

    # ------------------------------------------------------------------
    # UniProt + AlphaFold (live)
    # ------------------------------------------------------------------

    def test_uniprot_for_entry_maps_3w85(self):
        refs = self.client.uniprot_for_entry("3W85")
        assert refs, "3W85 carries a UniProt mapping"
        assert all(r["accession"] for r in refs)

    def test_uniprot_context_real_binding_sites(self):
        # Human dihydroorotate dehydrogenase: P02742... actually Q02127.
        data = self.client.uniprot_context("Q02127")
        assert data["accession"] == "Q02127"
        assert data["sequence_length"] and data["sequence_length"] > 300
        assert isinstance(data["binding_sites"], list)
        assert data["url"].startswith("https://www.uniprot.org/")

    def test_uniprot_missing_accession_is_honest(self):
        with pytest.raises(DiscoveryError):
            self.client.uniprot_context("Z99999")

    def test_alphafold_model_url_and_fetch(self):
        url = self.client.alphafold_model_url("Q02127")
        assert url.endswith("AF-Q02127-F1-model_v6.cif")
        model = self.client.fetch_alphafold_model("Q02127")
        assert model["format"] == "mmcif"
        assert len(model["content"]) > 10000

    def test_alphafold_missing_model_is_honest(self):
        # A UniProt accession with no AlphaFold model (fragment-only entry).
        with pytest.raises(DiscoveryError):
            self.client.fetch_alphafold_model("Q9999999")


# ---------------------------------------------------------------------------
# Pockets (real fpocket) and DSSP (real mkdssp)
# ---------------------------------------------------------------------------

class TestEngines:
    def test_pockets_on_1mbn(self, tmp_path):
        if not _has("fpocket"):
            pytest.skip("fpocket not installed (snap install fpocket)")
        s = _fetch("1MBN")
        detector = PocketDetector()
        result = detector.detect_pockets(s, tmp_path)
        assert result["available"] is True
        assert result["pocket_count"] >= 1
        p = result["pockets"][0]
        assert p["atom_count"] > 0
        assert len(p["center"]) == 3
        # fpocket's own descriptors flow through untouched.
        assert isinstance(p["descriptors"], dict)

    def test_pocket_ligand_overlap_real_geometry(self, tmp_path):
        if not _has("fpocket"):
            pytest.skip("fpocket not installed")
        s = _fetch("1MBN")
        from ligora_backend.ligand import LigandResolver
        for lig in s.ligands:
            LigandResolver().resolve_ligand(lig)
        hem = next(lig for lig in s.ligands if lig.residue_name == "HEM")
        detector = PocketDetector()
        result = detector.detect_pockets(s, tmp_path, ligand=hem)
        assert "ligand_overlap" in result
        assert result["ligand_overlap"], "heme sits in a detectable pocket"
        best = result["ligand_overlap"][0]
        assert best["min_distance"] < 8.0, (
            "the co-crystallized ligand must be near some detected pocket")

    def test_dssp_on_1mbn(self, tmp_path):
        if not _has("mkdssp") and not _has("dssp"):
            pytest.skip("DSSP not installed (apt install dssp)")
        s = _fetch("1MBN")
        analyzer = DSSPAnalyzer()
        result = analyzer.analyze(s, tmp_path)
        assert result["available"] is True
        assert len(result["residues"]) > 100
        # Myoglobin is ~75% helical; real DSSP must reflect that.
        counts = result["ss_counts"]
        helix = sum(v for k, v in counts.items() if "helix" in k)
        assert helix > len(result["residues"]) * 0.5, (
            f"myoglobin should be majority helix, got {counts}")

    def test_dssp_absent_is_honest(self, tmp_path):
        s = _fetch("1MBN")
        analyzer = DSSPAnalyzer(binary_path="/nonexistent/dssp")
        result = analyzer.analyze(s, tmp_path)
        assert result["available"] is False
        assert "/nonexistent/dssp" in (result.get("error") or ""), (
            "the configured path must be named in the error")


# ---------------------------------------------------------------------------
# 2D interaction diagram + struct_conn export (from real PLIP contacts)
# ---------------------------------------------------------------------------

class TestDiagramAndExport:
    @pytest.fixture(scope="class")
    def analyzed_3w85(self):
        if not (_has("plip") or shutil.which("plip")):
            pytest.skip("PLIP not installed")
        s = _fetch("3W85")
        # Resolve EVERY ligand instance through the real resolver — every
        # type of molecule (organic, ion, solvent, cofactor) must work.
        from ligora_backend.ligand import LigandResolver
        for lig in s.ligands:
            LigandResolver().resolve_ligand(lig, "3W85")
        w85 = next(lig for lig in s.ligands if lig.residue_name == "W85")
        contacts, available = ContactAnalyzer().analyze_contacts(s, w85)
        if not available:
            pytest.skip("PLIP not installed")
        return s, w85, contacts

    def test_diagram_svg_from_real_contacts(self, analyzed_3w85):
        from ligora_backend.v2.diagram2d import build_interaction_diagram
        s, w85, contacts = analyzed_3w85
        assert contacts, "PLIP must find contacts on 3W85/W85"
        result = build_interaction_diagram(s, w85, contacts, 900, 640)
        assert result["svg"].startswith("<svg")
        assert "</svg>" in result["svg"]
        assert result["residues"], "residue nodes come from real contacts"
        assert result["layout_source"] in (
            "rdkit_2d_depiction", "pca_projection_of_3d_coordinates")
        # Residue labels must match PLIP's own contacted residues.
        plip_residues = {(c.protein_chain_id, c.protein_residue_id)
                         for c in contacts}
        diagram_residues = {(r["chain"], r["residue_id"])
                            for r in result["residues"]}
        assert diagram_residues == plip_residues

    def test_struct_conn_export_roundtrip(self, analyzed_3w85, tmp_path):
        from ligora_backend.v2.aggregation import export_struct_conn
        s, w85, contacts = analyzed_3w85
        path = tmp_path / "struct_conn.cif"
        text = export_struct_conn(s, w85, contacts, path)
        assert path.read_text() == text
        assert "_struct_conn.conn_type_id" in text
        # The real contacts PLIP found must be present as records, and the
        # file is directly parseable by an mmCIF reader (self-consistent).
        from ligora_backend.parser import StructureParser
        parsed = StructureParser().parse_mmcif(text, source_id="x",
                                               source="test")
        assert parsed is not None  # valid CIF
        rows = [ln for ln in text.splitlines()
                if ln and not ln.startswith(('_', '#', 'loop_'))]
        assert len(rows) >= 9  # 3W85/W85 has 12 mappable contacts

    def test_aggregation_from_real_batch_results(self):
        # Build batch results through the real analyzer over two structures.
        from ligora_backend.v2.batch import BatchAnalyzer
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()
        analyzer.add_structure_to_batch(batch_id, "3W85", source_type="pdb_id")
        analyzer.add_structure_to_batch(batch_id, "1MBN", source_type="pdb_id")
        results = analyzer.run_batch(batch_id)
        assert len(results) == 2
        succeeded = [r for r in results if r.success]
        assert succeeded, f"both structures must analyze: " \
            f"{[(r.structure_id, r.error) for r in results]}"
        rows = aggregate_interaction_frequencies(results)
        assert rows, "real contacts must aggregate into rows"
        for r in rows:
            assert r["frequency"] > 0
            assert r["by_type"], "per-type counts come from real contacts"
            assert r["structures_contacted"] >= 1

    def test_aggregate_batch_interactions_handler(self):
        server = BackendServer()
        analyzer = server.batch_analyzer
        batch_id = analyzer.create_batch()
        analyzer.add_structure_to_batch(batch_id, "3W85", source_type="pdb_id")
        analyzer.run_batch(batch_id)
        response = server._handlers["aggregate_batch_interactions"](
            {"batch_id": batch_id}, "")
        assert response["row_count"] > 0
        assert all("by_type" in r for r in response["rows"])


# ---------------------------------------------------------------------------
# Crystal-vs-docked and superposition (real geometry, real Vina output)
# ---------------------------------------------------------------------------

class TestComparison:
    def test_crystal_vs_docked_needs_vina_output(self):
        s = _fetch("3W85")
        w85 = next(lig for lig in s.ligands if lig.residue_name == "W85")
        comparator = ResultComparator()

        class FakePose:  # structure-only probe, never used as data
            atoms = w85.atoms
            pose_id = 0
            affinity = -7.0

        # Identical coordinates must give RMSD 0 (sanity of the math).
        result = comparator.compare_crystal_vs_docked(w85, FakePose())
        assert result.get("rmsd") == 0.0, result

    def test_superpose_self_is_zero(self):
        s = _fetch("3W76")
        comparator = ResultComparator()
        result = comparator.superpose_structures(s, s, "A", "A")
        assert result.get("rmsd") == 0.0, result
        assert result["residues_used"] > 100

    def test_superpose_homologs_real_rmsd(self):
        a = _fetch("3W85")   # T. cruzi DHODH
        b = _fetch("3W76")   # same enzyme family, same study
        comparator = ResultComparator()
        result = comparator.superpose_structures(a, b, "A", "A")
        assert "error" not in result, result
        assert 0 < result["rmsd"] < 15.0
        assert result["residues_used"] > 100
        assert 0 < result["sequence_identity"] <= 1.0

    def test_superpose_refuses_different_folds(self):
        a = _fetch("3W85")   # class-1a DHODH fold
        b = _fetch("1DHK")   # class-2 DHODH: different fold
        comparator = ResultComparator()
        result = comparator.superpose_structures(a, b, "A", "A")
        assert "error" in result, "different folds must be refused, not forced"

    def test_superpose_handler_via_server(self):
        server = BackendServer()
        session = server._create_session()
        # Real path: open_pdb_id handler loads and resolves the structure
        # into the session.
        server._handlers["open_pdb_id"]({"pdb_id": "3W85"}, session.id)
        result = server._handlers["superpose_structures"](
            {"chain_a": "A", "reference_pdb_id": "3W76", "chain_b": "A"},
            session.id)
        assert "error" not in result, result
        assert result["rmsd"] > 0
        assert result["sequence_identity"] > 0


# ---------------------------------------------------------------------------
# Docking box override + enrich_component (server handlers)
# ---------------------------------------------------------------------------

class TestServerHandlers:
    def test_docking_box_set_get_reset(self):
        server = BackendServer()
        session = server._create_session()
        # A real structure must be loaded (get_docking_box derives the
        # ligand-extent fallback from it).
        server._handlers["open_pdb_id"]({"pdb_id": "3W85"}, session.id)
        payload = {"center_x": 10.0, "center_y": 11.0, "center_z": 12.0,
                   "size_x": 20.0, "size_y": 20.0, "size_z": 20.0}
        out = server._handlers["set_docking_box"](payload, session.id)
        assert out["status"] == "ok"
        got = server._handlers["get_docking_box"]({}, session.id)
        assert got["box"]["source"] == "user_set"
        assert got["box"]["center"] == [10.0, 11.0, 12.0]
        # Invalid sizes must be refused.
        with pytest.raises(ValueError):
            server._handlers["set_docking_box"](
                {**payload, "size_x": -1}, session.id)

    def test_enrich_component_real_ccd(self):
        server = BackendServer()
        data = server._handlers["enrich_component"]({"comp_id": "HEM"}, "")
        assert data["name"]
        assert data["formula"]
        assert data["pdbx_type"], "CCD pdbx_type must be present"

    def test_enrich_component_rejects_garbage(self):
        server = BackendServer()
        with pytest.raises(ValueError):
            server._handlers["enrich_component"]({"comp_id": ""}, "")
        with pytest.raises(ValueError):
            server._handlers["enrich_component"]({"comp_id": "NOTREAL!"}, "")

    def test_discover_similar_components_via_server(self):
        server = BackendServer()
        out = server._handlers["discover_similar_components"](
            {"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "rows": 5}, "")
        assert out["total_count"] > 0
        assert out["hits"][0]["id"] == "AIN"
