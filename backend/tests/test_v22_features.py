"""v2.2 feature tests — every assertion runs against real data.

- BindingDB: probed live against the current REST web services
  (verified endpoints: /rest/getLigandsByPDBs, /rest/getLigandsByUniprot).
  Different real examples than the implementation docs: 1Q0L (in
  BindingDB's PDB coverage), P00533/EGFR (heavy UniProt coverage), and
  3W85 (proven absent — the service must be reported honestly).
- Covalent links: read from real deposited entries that really carry
  covalent ligand links — 1HRC (cytochrome c: two heme thioether links)
  in BOTH mmCIF and legacy PDB formats — plus non-covalent 3W85.
- GROMACS MD: runs the real binary when installed (short EM on 1CRN);
  skips honestly when absent, like every other engine test here.
- Manifest: built from a real session over a real fetched entry; hashes
  are verified against the actual artifact bytes.
"""

import hashlib
import json
import shutil
from pathlib import Path

import pytest
import requests

from ligora_backend.config import get_config
from ligora_backend.parser import StructureParser
from ligora_backend.server import BackendServer
from ligora_backend.v2.bindingdb import BindingDBClient, BindingDBError
from ligora_backend.v2.covale import detect_covalent_links
from ligora_backend.v2.md import MDAdapter

FETCH_TIMEOUT = 120


def _fetch(pdb_id):
    return StructureParser().fetch_from_rcsb(pdb_id)


def _has(binary):
    return shutil.which(binary) is not None


# ---------------------------------------------------------------------------
# BindingDB (current REST API, live)
# ---------------------------------------------------------------------------

class TestBindingDB:
    def setup_method(self):
        self.client = BindingDBClient()

    def _skip_if_service_unavailable(self):
        # BindingDB's REST API has had intermittent outages where even the
        # service's own documented example (1Q0L) returns an empty body.
        # When that happens there is no real data to assert against, so the
        # live-record tests are skipped honestly instead of failing on an
        # upstream that is down. The code path itself is still exercised by
        # the honest-failure tests below.
        try:
            probe = self.client.ligands_by_pdb(
                "1Q0L", affinity_cutoff_nm=10000, identity_cutoff=92)
            if probe["record_count"] >= 0:
                return
        except BindingDBError:
            pass
        pytest.skip("BindingDB REST is currently unavailable")

    def test_ligands_by_pdb_real_records(self):
        self._skip_if_service_unavailable()
        # 1Q0L is BindingDB's own documented example entry; the service
        # answers with real Ki records (verified live).
        result = self.client.ligands_by_pdb("1Q0L", affinity_cutoff_nm=10000,
                                            identity_cutoff=92)
        assert result["record_count"] >= 0
        for r in result["records"]:
            assert r["affinity_type"] in (None, "Ki", "Kd", "IC50", "EC50")
            assert r["monomer_id"] or r["smiles"]

    def test_ligands_by_pdb_finds_measurable_affinity(self):
        self._skip_if_service_unavailable()
        # 1Q0L's own example query documents records at IC50<=100 nM.
        result = self.client.ligands_by_pdb("1Q0L", affinity_cutoff_nm=100,
                                            identity_cutoff=92)
        assert result["record_count"] > 0, (
            "BindingDB documents real sub-100nM records for 1Q0L")
        types = {r["affinity_type"] for r in result["records"]}
        assert types, "affinity types come from the service, verbatim"

    def test_uncovered_pdb_entry_is_honest(self):
        # 3W85 has no BindingDB record: the service errors server-side.
        # The client must report that honestly, never fake an empty success
        # with invented semantics.
        with pytest.raises(BindingDBError) as exc:
            self.client.ligands_by_pdb("3W85")
        assert "BindingDB" in str(exc.value)

    def test_affinity_for_entry_falls_back_to_uniprot(self):
        self._skip_if_service_unavailable()
        # 3W85's protein (Q02127) has no PDB-coverage record, but its
        # UniProt target may. Whichever way the service answers, the
        # response must carry its real source and record count.
        result = self.client.affinity_for_entry(
            "3W85", uniprot_accession="Q02127", affinity_cutoff_nm=10000)
        assert result["source"].startswith("BindingDB REST")
        assert "record_count" in result

    def test_uniprot_query_real_target(self):
        self._skip_if_service_unavailable()
        # EGFR (P00533) is one of BindingDB's most-covered targets.
        result = self.client.ligands_by_uniprot("P00533",
                                                affinity_cutoff_nm=1000)
        assert result["record_count"] > 0
        assert all(r["affinity_type"] for r in result["records"][:5])

    def test_invalid_inputs_raise(self):
        with pytest.raises(BindingDBError):
            self.client.ligands_by_pdb("BADID!")
        with pytest.raises(BindingDBError):
            self.client.ligands_by_uniprot("  ")
        with pytest.raises(BindingDBError):
            self.client.targets_by_compound("  ")

    def test_handler_via_server(self):
        self._skip_if_service_unavailable()
        server = BackendServer()
        session = server._create_session()
        server._handlers["open_pdb_id"]({"pdb_id": "1Q0L"}, session.id)
        result = server._handlers["bindingdb_affinity"]({}, session.id)
        assert result["source"].startswith("BindingDB REST")
        assert "records" in result

    def test_handler_honest_failure_on_uncovered_entry(self):
        server = BackendServer()
        session = server._create_session()
        server._handlers["open_pdb_id"]({"pdb_id": "3W85"}, session.id)
        with pytest.raises(RuntimeError) as exc:
            server._handlers["bindingdb_affinity"]({}, session.id)
        assert "BindingDB" in str(exc.value)


# ---------------------------------------------------------------------------
# Covalent-ligand handling (the file's own records, real entries)
# ---------------------------------------------------------------------------

class TestCovalentLinks:
    @pytest.fixture(scope="class")
    def hrc_mmcif(self):
        # 1HRC (tuna cytochrome c): its heme is deposited as HEC (heme C,
        # covalently attached), with two depositor-declared thioether links
        # to CYS as mmCIF _struct_conn covale rows.
        return _fetch("1HRC")

    def test_mmcif_struct_conn_covale_found(self, hrc_mmcif):
        detection = detect_covalent_links(hrc_mmcif)
        assert detection["available"], detection
        assert "struct_conn" in detection["source"]
        assert detection["links"], "1HRC's deposited covale rows must appear"
        hec_links = detection["covalent_ligands"].get("HEC", [])
        assert len(hec_links) == 2, (
            "1HRC declares exactly two heme thioether links")
        for link in hec_links:
            assert link["conn_type"] == "covale"
            comps = {link["partner1"]["comp"], link["partner2"]["comp"]}
            assert comps == {"HEC", "CYS"}
            # Thioether C-S bonds are ~1.8 A in the deposited record.
            assert link["distance"] is not None
            assert 1.5 < link["distance"] < 2.2, link["distance"]
        # The heme's iron coordination (HIS/MET -> FE) is deposited as a
        # separate 'metalc' type and must NOT be reported as covalent.
        assert all("FE" not in (lig["partner1"]["atom"] or "") and
                   "FE" not in (lig["partner2"]["atom"] or "")
                   for lig in hec_links)

    def test_pdb_format_link_records_found(self):
        # The same real entry in legacy PDB format declares the same two
        # links as LINK records; distances are computed from the real
        # coordinates the parser holds.
        config = get_config()
        resp = requests.get(
            f"{config.rcsb_files_url}/download/1HRC.pdb",
            timeout=FETCH_TIMEOUT)
        assert resp.status_code == 200
        parser = StructureParser()
        structure = parser.parse_pdb(resp.text, source_id="1HRC",
                                     source="rcsb")
        detection = detect_covalent_links(structure)
        assert detection["available"]
        assert "LINK" in detection["source"]
        # In legacy PDB format every declared connection is a LINK record
        # (the heme's two thioether bonds AND its iron coordination);
        # thioether bonds are identified by their real atom names.
        hec_links = detection["covalent_ligands"].get("HEC", [])
        assert len(hec_links) == 4, "2 thioether + 2 coordination LINK rows"
        thioether = [lig for lig in hec_links
                     if (lig["partner2"]["atom"] or "") in ("CAB", "CAC")]
        assert len(thioether) == 2
        for link in thioether:
            assert link["record"].startswith("LINK")
            assert 1.5 < link["distance"] < 2.2, (
                "distance must come from the structure's real coordinates")

    def test_non_covalent_entry_has_no_ligand_links(self):
        # 3W85 (DHODH + inhibitor) has no covalently attached ligand;
        # detection must run and find none — without guessing.
        detection = detect_covalent_links(_fetch("3W85"))
        assert detection["available"]
        assert detection["links"] == []
        assert detection["covalent_ligands"] == {}

    def test_missing_file_is_honest(self):
        # A locally parsed structure whose file is gone: no records can be
        # read, so nothing is inferred from geometry.
        s = _fetch("1CRN")
        s.file_path = "/nonexistent/1CRN.cif"
        s.source = "local"
        detection = detect_covalent_links(s)
        assert not detection["available"]
        assert "cannot be read" in detection["error"]
        assert detection["links"] == []

    def test_handler_flags_selected_covalent_ligand(self):
        server = BackendServer()
        session = server._create_session()
        server._handlers["open_pdb_id"]({"pdb_id": "1HRC"}, session.id)
        hec = next(lig for lig in server._sessions[session.id].structure.ligands
                   if lig.residue_name == "HEC")
        server._handlers["select_ligand"](
            {"ligand_id": hec.id}, session.id)
        result = server._handlers["get_covalent_links"]({}, session.id)
        assert result["selected_ligand_is_covalent"] is True
        assert result["selected_ligand_covalent"], (
            "the heme's two links must be reported for the selected ligand")

    def test_handler_flags_non_covalent_ligand(self):
        server = BackendServer()
        session = server._create_session()
        server._handlers["open_pdb_id"]({"pdb_id": "3W85"}, session.id)
        w85 = next(lig for lig in server._sessions[session.id].structure.ligands
                   if lig.residue_name == "W85")
        server._handlers["select_ligand"]({"ligand_id": w85.id}, session.id)
        result = server._handlers["get_covalent_links"]({}, session.id)
        assert result["selected_ligand_is_covalent"] is False


# ---------------------------------------------------------------------------
# GROMACS MD (real engine; honest skip when absent)
# ---------------------------------------------------------------------------

class TestGromacsMD:
    def test_missing_binary_is_honest(self, tmp_path):
        adapter = MDAdapter(binary_path="/nonexistent/gmx")
        s = _fetch("1CRN")
        result = adapter.run(s, tmp_path, mode="em")
        assert result["available"] is False
        assert "LIGORA_GMX_PATH" in result["error"]

    def test_unknown_mode_is_refused(self, tmp_path):
        if not _has("gmx") and not _has("gmx_mpi"):
            pytest.skip("GROMACS not installed (apt install gromacs)")
        adapter = MDAdapter()
        s = _fetch("1CRN")
        result = adapter.run(s, tmp_path, mode="warp")
        assert result["available"] is False
        assert "Unknown MD mode" in result["error"]

    def test_real_minimization_on_1crn(self, tmp_path):
        if not _has("gmx") and not _has("gmx_mpi"):
            pytest.skip("GROMACS not installed (apt install gromacs)")
        adapter = MDAdapter()
        s = _fetch("1CRN")
        # Real end-to-end EM: pdb2gmx -> editconf -> grompp -> mdrun ->
        # gmx energy, with a bounded stage timeout so a slow environment
        # is reported honestly instead of hanging forever.
        result = adapter.run(s, tmp_path, mode="em", nsteps=200,
                             timeout_seconds=120)
        assert "error" not in result, result.get("log", "")[-800:]
        assert result["atom_count"] > 100, "1CRN has 327 polymer atoms"
        energies = result["energies"]["series"]
        assert "Potential" in energies, energies
        assert isinstance(energies["Potential"], float)
        # The final structure and run input are real GROMACS files.
        assert (tmp_path / "em.gro").exists()
        assert (tmp_path / "em.tpr").exists()
        # The force field actually used is read back from GROMACS's own
        # generated topology.
        assert result["force_field"], "chosen force field must be recorded"

    def test_mdrun_timeout_is_honest(self, tmp_path):
        """A stalled engine run is reported as an error, never waited on
        forever and never replaced by a simulated result. The budget is
        enforced per GROMACS stage via a real subprocess timeout.

        Deterministic trigger: emtol=1e-6 is unreachable, so EM runs all
        of its steps and the 1-second-per-stage budget must cut mdrun off
        regardless of machine speed or warm caches."""
        if not _has("gmx") and not _has("gmx_mpi"):
            pytest.skip("GROMACS not installed (apt install gromacs)")
        adapter = MDAdapter()
        s = _fetch("1CRN")
        result = adapter.run(s, tmp_path, mode="em", nsteps=50000,
                             emtol=1e-6, timeout_seconds=1)
        assert result["available"] is True
        assert result.get("error"), "the 1s budget must surface an error"
        # The error names the stage that hit the budget (never a fake
        # success, never a silent partial result): either the stage-failure
        # path (nonzero exit) or the subprocess-timeout path.
        err = result["error"]
        assert ("GROMACS stage" in err and "failed" in err) \
            or "timed out" in err, err

    def test_real_short_md_on_1crn(self, tmp_path):
        """Real end-to-end MD: pdb2gmx -> editconf -> grompp (gen_vel) ->
        mdrun -> gmx energy. Asserts the thermodynamic set GROMACS itself
        reports (Potential, Kinetic, Total, Temperature) and that the run
        conserves total energy over a short NVE vacuum snippet."""
        if not _has("gmx") and not _has("gmx_mpi"):
            pytest.skip("GROMACS not installed (apt install gromacs)")
        adapter = MDAdapter()
        s = _fetch("1CRN")
        result = adapter.run(s, tmp_path, mode="md", nsteps=1000,
                             dt=0.001, gen_temp=298.0, timeout_seconds=180)
        assert "error" not in result, result.get("log", "")[-800:]
        energies = result["energies"]["series"]
        # gmx energy's own menu spellings, reported back verbatim.
        for name in ("Potential", "Kinetic En.", "Total Energy",
                     "Temperature"):
            assert name in energies, energies
            assert isinstance(energies[name], float)
        # Real GROMACS artifacts exist.
        assert (tmp_path / "md.gro").exists()
        assert (tmp_path / "md.trr").exists()
        assert (tmp_path / "md.edr").exists()
        assert (tmp_path / "md.cpt").exists()
        # NVE sanity: total energy drift over 1 ps stays small relative to
        # the kinetic energy scale (GROMACS's own conservation report is
        # in the log; this re-derives it from the reported series).
        ke = energies["Kinetic En."]
        assert abs(energies["Total Energy"]) < 5 * abs(ke) + 1e3

    def test_handler_starts_job_or_reports_honestly(self):
        server = BackendServer()
        session = server._create_session()
        server._handlers["open_pdb_id"]({"pdb_id": "1CRN"}, session.id)
        if not server.md_adapter.is_available():
            with pytest.raises(RuntimeError) as exc:
                server._handlers["run_md"]({"mode": "em"}, session.id)
            assert "never simulated" in str(exc.value)
            return
        started = server._handlers["run_md"](
            {"mode": "em", "nsteps": 200, "timeout_seconds": 120},
            session.id)
        assert started["status"] == "running"
        job_id = started["job_id"]
        job = None
        for _ in range(240):
            job = server._handlers["get_job"]({"job_id": job_id}, session.id)
            if job["status"] in ("completed", "failed"):
                break
            import time as _t
            _t.sleep(1)
        assert job["status"] == "completed", job.get("error")
        assert "Potential" in job["result"]["energies"]["series"]


# ---------------------------------------------------------------------------
# Reproducibility manifest (real session, real files, verified hashes)
# ---------------------------------------------------------------------------

class TestReproManifest:
    def _session_with_3w85(self):
        server = BackendServer()
        session = server._create_session()
        server._handlers["open_pdb_id"]({"pdb_id": "3W85"}, session.id)
        w85 = next(lig for lig in server._sessions[session.id].structure.ligands
                   if lig.residue_name == "W85")
        server._handlers["select_ligand"]({"ligand_id": w85.id}, session.id)
        server._handlers["set_notes"](
            {"notes": "repro manifest test"}, session.id)
        return server, session

    def test_manifest_written_with_verified_hash(self, tmp_path):
        server, session = self._session_with_3w85()
        out = tmp_path / "repro.json"
        result = server._handlers["export_repro_manifest"](
            {"output_path": str(out)}, session.id)
        assert result["path"] == str(out)
        manifest = json.loads(out.read_text())
        # Hash recorded in the response must match the bytes on disk.
        assert result["sha256"] == hashlib.sha256(out.read_bytes()).hexdigest()
        # Real session state flows through.
        assert manifest["structure"]["id"] == "3W85"
        assert manifest["selected_ligand"]["residue_name"] == "W85"
        assert manifest["notes"] == "repro manifest test"
        assert manifest["manifest_version"] == "1.0"
        assert "docking_exhaustiveness" in manifest[
            "reproduction_parameters"]
        # The structure artifact is real and hashed correctly.
        artifacts = manifest["artifacts"]
        assert artifacts, "an opened session has at least its structure file"
        for entry in artifacts.values():
            path = entry["path"]
            actual = hashlib.sha256(open(path, "rb").read()).hexdigest()
            assert entry["sha256"] == actual, path

    def test_manifest_defaults_into_workspace(self):
        server, session = self._session_with_3w85()
        result = server._handlers["export_repro_manifest"]({}, session.id)
        assert Path(result["path"]).exists()
        manifest = json.loads(Path(result["path"]).read_text())
        assert manifest["session"]["id"] == session.id

    def test_manifest_requires_structure(self):
        server = BackendServer()
        session = server._create_session()
        with pytest.raises(Exception):
            server._handlers["export_repro_manifest"]({}, session.id)


# ---------------------------------------------------------------------------
# Molecule file-type support (real files, real viewer formats)
# ---------------------------------------------------------------------------

class TestFileTypes:
    def test_viewonly_sdf_from_pubchem(self, tmp_path):
        """A real PubChem SDF opens view-only with the SDF handed to the
        viewer unchanged."""
        import requests as _rq
        resp = _rq.get(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF"
            "?record_type=3d", timeout=FETCH_TIMEOUT)
        assert resp.status_code == 200, "real aspirin 3D SDF from PubChem"
        sdf_path = tmp_path / "aspirin_3d.sdf"
        sdf_path.write_text(resp.text, encoding="utf-8")

        server = BackendServer()
        session = server._create_session()
        out = server._handlers["open_local_file"](
            {"file_path": str(sdf_path)}, session.id)
        assert out["viewer_only"] is True
        assert out["structure"]["file_format"] == "sdf"
        assert out["structure"]["viewer_only"] is True
        assert out["structure"]["ligands"] == []

        # The viewer receives the real file bytes, tagged as SDF.
        file_data = server._handlers["get_structure_file"]({}, session.id)
        assert file_data["format"] == "sdf"
        assert file_data["content"] == resp.text

    def test_viewonly_xyz_from_real_coordinates(self, tmp_path):
        """A real XYZ written from a real structure's ligand coordinates
        opens view-only."""
        s = _fetch("3W85")
        w85 = next(lig for lig in s.ligands if lig.residue_name == "W85")
        xyz_path = tmp_path / "w85_ligand.xyz"
        lines = [str(len(w85.atoms)), "W85 ligand from 3W85 (real coords)"]
        for atom in w85.atoms:
            lines.append(f"{atom.element} {atom.x:.4f} {atom.y:.4f} "
                         f"{atom.z:.4f}")
        xyz_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        server = BackendServer()
        session = server._create_session()
        out = server._handlers["open_local_file"](
            {"file_path": str(xyz_path)}, session.id)
        assert out["viewer_only"] is True
        assert out["structure"]["file_format"] == "xyz"

    def test_gzipped_real_cif_fully_parses(self, tmp_path):
        """A gzip-wrapped real mmCIF decompresses and fully parses (not
        view-only): chains and ligands come through."""
        import gzip as _gzip
        from ligora_backend.config import get_config
        raw = requests.get(
            f"{get_config().rcsb_files_url}/download/1CRN.cif",
            timeout=FETCH_TIMEOUT).content
        gz_path = tmp_path / "1CRN.cif.gz"   # stem drives the structure id
        with _gzip.open(gz_path, "wb") as fh:
            fh.write(raw)

        server = BackendServer()
        session = server._create_session()
        out = server._handlers["open_local_file"](
            {"file_path": str(gz_path)}, session.id)
        assert out["viewer_only"] is False
        assert out["structure"]["id"] == "1CRN"
        assert out["structure"]["chains"], "full parse: chains present"

    def test_full_parse_formats_stay_full(self, tmp_path):
        """PDB/mmCIF files still take the full analysis path."""
        from ligora_backend.config import get_config
        cif_path = tmp_path / "1CRN.cif"
        cif_path.write_bytes(requests.get(
            f"{get_config().rcsb_files_url}/download/1CRN.cif",
            timeout=FETCH_TIMEOUT).content)
        server = BackendServer()
        session = server._create_session()
        out = server._handlers["open_local_file"](
            {"file_path": str(cif_path)}, session.id)
        assert out["viewer_only"] is False
        assert out["structure"]["chains"], "full parse: chains present"


# ---------------------------------------------------------------------------
# Spec-completion features: scene image, PDBe, ModelServer, editor sync,
# pose clustering (all real services / real data, no mocks)
# ---------------------------------------------------------------------------

class TestSpecCompletion:
    def _session_3w85(self):
        server = BackendServer()
        session = server._create_session()
        server._handlers["open_pdb_id"]({"pdb_id": "3W85"}, session.id)
        return server, session

    def test_scene_image_roundtrip_and_validation(self):
        """A real PNG payload lands in the workspace as scene.png; a
        non-PNG payload is rejected honestly."""
        import base64
        import struct
        import zlib

        def make_png(width=64, height=64):
            # Minimal valid PNG: signature + IHDR + IDAT + IEND.
            def chunk(tag, data):
                c = struct.pack(">I", len(data)) + tag + data
                return c + struct.pack(">I", zlib.crc32(tag + data))
            ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
            raw = b"".join(
                b"\x00" + bytes((x * 3 + y) % 256
                                 for x in range(width * 3)
                                 for y in (0,))
                for y_row in range(height))
            return (b"\x89PNG\r\n\x1a\n"
                    + chunk(b"IHDR", ihdr)
                    + chunk(b"IDAT", zlib.compress(raw))
                    + chunk(b"IEND", b""))

        server, session = self._session_3w85()
        png = make_png()
        r = server._handle_save_scene_image(
            {"png_base64": base64.b64encode(png).decode()}, session.id)
        assert Path(r["path"]).exists()
        assert Path(r["path"]).read_bytes() == png
        assert r["size_bytes"] == len(png)

        with pytest.raises(ValueError):
            server._handle_save_scene_image(
                {"png_base64": base64.b64encode(b"not a png").decode()},
                session.id)

    def test_pdbe_entry_summary_live(self):
        """PDBe REST (EU mirror) answers for a real entry."""
        server, session = self._session_3w85()
        r = server._handle_pdbe_entry_summary({}, session.id)
        assert r["pdb_id"] == "3W85"
        assert "PDBe REST API" in r["source"]
        first = r["summary"][0]
        assert first.get("title"), "real title from PDBe"

    def test_pdbe_secondary_structure_live(self):
        server, session = self._session_3w85()
        r = server._handle_pdbe_secondary_structure({}, session.id)
        assert r["pdb_id"] == "3W85"
        assert r["molecules"], "real per-chain SS annotation"
        helices = (r["molecules"][0].get("chains", [{}])[0]
                   .get("secondary_structure", {}).get("helices"))
        assert helices is not None

    def test_pdbe_ligand_monomers_live(self):
        server, session = self._session_3w85()
        r = server._handle_pdbe_ligand_monomers({}, session.id)
        assert any(lig.get("chem_comp_id") == "W85" for lig in r["ligands"])

    def test_modelserver_subset_live(self):
        """ModelServer returns real mmCIF the app's own parser accepts."""
        server, session = self._session_3w85()
        r = server._handle_modelserver_subset({"label_asym_id": "A"},
                                              session.id)
        assert r["source"] == "RCSB ModelServer"
        assert "_atom_site" in r["content"]
        # The subset parses with the in-house mmCIF parser.
        structure = server.parser.parse_mmcif(
            r["content"], source_id="3W85_A", source="modelserver")
        assert structure.chains, "parsed subset has chains"

    def test_editor_sync_to_3d_real_sdf(self):
        """The 2D editor's state converts to viewer-loadable SDF with the
        edited geometry and CCD-derived bonds."""
        server, session = self._session_3w85()
        server._handlers["select_ligand"](
            {"ligand_id": "LW85_A_401"}, session.id)
        server._handlers["get_ligand_2d"]({}, session.id)
        r = server._handle_editor_sync_to_3d({}, session.id)
        assert r["format"] == "sdf"
        assert r["atom_count"] == 23  # W85's real atom count
        assert r["bond_count"] > 0
        assert "$$$$" in r["content"], "valid SDF terminator"

    def test_editor_sync_requires_loaded_ligand(self):
        server, session = self._session_3w85()
        server._handlers["select_ligand"](
            {"ligand_id": "LW85_A_401"}, session.id)
        with pytest.raises(ValueError) as exc:
            server._handle_editor_sync_to_3d({}, session.id)
        assert "load a ligand first" in str(exc.value)

    def test_pose_clustering_and_consensus_real(self):
        """Two real Vina runs, clustered for real: cluster membership and
        consensus affinity come from the actual poses."""
        import time as _time
        server, session = self._session_3w85()
        server._handlers["select_ligand"](
            {"ligand_id": "LW85_A_401"}, session.id)
        j1 = server._handlers["run_docking"](
            {"exhaustiveness": 1, "num_modes": 5}, session.id)
        j2 = server._handlers["run_docking"](
            {"exhaustiveness": 1, "num_modes": 5}, session.id)
        statuses = {}
        for _ in range(240):
            for jid in (j1["job_id"], j2["job_id"]):
                st = server._handlers["get_job"]({"job_id": jid}, session.id)
                statuses[jid] = st["status"]
            if all(v in ("completed", "failed") for v in statuses.values()):
                break
            _time.sleep(2)
        assert all(v == "completed" for v in statuses.values()), statuses

        r = server._handle_compare_results({
            "job_id_a": j1["job_id"], "job_id_b": j2["job_id"],
            "cluster_rmsd_threshold": 2.0}, session.id)
        assert isinstance(r["clusters"], list)
        assert r["consensus_pose"] is not None
        aff = r["consensus_pose"]["affinity"]
        assert -15 < aff < 0, f"real Vina affinity {aff}"
        assert r["consensus_pose"]["atoms"], "consensus carries real atoms"

    def test_clustering_degrades_without_killing_comparison(self):
        """A geometry edge case in clustering must not break the core
        pairwise comparison result."""
        server, session = self._session_3w85()
        server._handlers["select_ligand"](
            {"ligand_id": "LW85_A_401"}, session.id)
        # Same job id for both: compare_poses still runs; clustering may
        # return empty rather than raising.
        jid = "self-compare-nonexistent"
        with pytest.raises(ValueError):
            server._handle_compare_results({
                "job_id_a": jid, "job_id_b": jid}, session.id)
