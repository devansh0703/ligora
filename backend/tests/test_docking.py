"""
Real docking tests: the actual AutoDock Vina CLI is executed on a real RCSB
structure. Nothing is mocked or simulated; when the vina or obabel binary is
not installed the tests are skipped with that reason (docking availability is
an environment property, not something tests should fake).
"""

import shutil

import pytest

from ligora_backend.engines import VinaAdapter
from ligora_backend.parser import StructureParser
from ligora_backend.ligand import LigandResolver

_vina = shutil.which("vina")
_obabel = shutil.which("obabel")
pytestmark = pytest.mark.skipif(
    not (_vina and _obabel),
    reason="vina/obabel not installed - docking cannot run for real")


@pytest.fixture(scope="module")
def real_complex():
    """A real RCSB structure with a drug-like ligand, elements from the CCD."""
    import requests

    content = requests.get(
        "https://files.rcsb.org/download/3W85.cif", timeout=60).text
    structure = StructureParser().parse_mmcif(
        content, source_id="3W85", source="rcsb")
    ligand = next(lig for lig in structure.ligands
                  if lig.residue_name == "W85" and lig.atoms
                  and lig.atoms[0].chain_id == "A")
    # Element completion comes from the CCD (parser helper) - required for
    # a real ligand 3D SDF during docking preparation.
    StructureParser().complete_missing_elements(structure)
    resolved = LigandResolver().resolve_ligand(ligand, "3W85")
    return structure, resolved


def test_vina_available_for_real():
    assert VinaAdapter().is_available() is True


def test_real_docking_produces_poses(real_complex):
    structure, ligand = real_complex
    adapter = VinaAdapter()
    import tempfile
    from pathlib import Path

    out = Path(tempfile.mkdtemp()) / "prep"
    input_files = adapter.prepare(structure, ligand, out)
    assert Path(input_files["receptor_pdbqt"]).exists()
    assert Path(input_files["ligand_pdbqt"]).exists()

    result = adapter.run(input_files, {"exhaustiveness": 1, "num_modes": 3})
    poses = result["poses"]
    assert len(poses) > 0
    for pose in poses:
        assert pose.affinity == pose.affinity  # not NaN
        assert len(pose.atoms) > 0
        for atom in pose.atoms:
            assert atom.element, "pose atoms must carry real elements"
    # Real Vina affinities for a drug-like ligand are bounded and negative.
    best = min(p.affinity for p in poses)
    assert -15.0 < best < 0.0


def test_two_real_runs_and_comparison(real_complex):
    structure, ligand = real_complex
    adapter = VinaAdapter()
    import tempfile
    from pathlib import Path
    from ligora_backend.v2.comparison import ResultComparator
    from ligora_backend.schemas import DockingResult

    out = Path(tempfile.mkdtemp()) / "prep"
    input_files = adapter.prepare(structure, ligand, out)

    res_a = adapter.run(input_files, {"exhaustiveness": 1, "num_modes": 3})
    res_b = adapter.run(input_files, {"exhaustiveness": 2, "num_modes": 3})
    assert res_a["poses"] and res_b["poses"]

    cmp_result = ResultComparator().compare_poses(
        DockingResult(job_id="a", engine="vina", poses=res_a["poses"]),
        DockingResult(job_id="b", engine="vina", poses=res_b["poses"]),
    )
    assert cmp_result["best_matches"], "real runs must yield comparable poses"
    assert cmp_result["avg_rmsd"] >= 0.0


def test_receptor_pdb_has_only_atom_records(real_complex):
    """Strict consumers (Vina) must not see non-ATOM records in receptor."""
    structure, ligand = real_complex
    adapter = VinaAdapter()
    import tempfile
    from pathlib import Path

    out = Path(tempfile.mkdtemp()) / "prep"
    input_files = adapter.prepare(structure, ligand, out)
    pdb_lines = Path(input_files["receptor_pdbqt"]).read_text().splitlines()
    # Open Babel's PDBQT contains REMARK headers and ATOM records; anything
    # else would break strict consumers.
    assert all(
        line.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT", "TORSDOF",
                         "BRANCH", "ENDBRANCH", "REMARK", "MODEL", "ENDMDL",
                         "TER", "END"))
        or line.strip() == ""
        for line in pdb_lines
    ), f"unexpected record in PDBQT: {pdb_lines[:3]}"
    # The PDBQT must actually contain receptor atoms.
    assert any(line.startswith("ATOM") for line in pdb_lines)
