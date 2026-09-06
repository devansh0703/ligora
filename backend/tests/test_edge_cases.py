"""
Edge-case tests across genuinely different real structures, all live from
RCSB. No mocks: each case fetches the real file and asserts honest behavior.

False-positive guard: water/ions must never be picked as "the ligand".
False-negative guard: real drug-like/cofactor ligands must be found.
"""

import requests

from ligora_backend.parser import StructureParser
from ligora_backend.ligand import LigandResolver
from ligora_backend.v2.water import WaterNetworkAnalyzer

FETCH_TIMEOUT = 90


def fetch(pdb_id: str):
    content = requests.get(
        f"https://files.rcsb.org/download/{pdb_id}.cif",
        timeout=FETCH_TIMEOUT).text
    return StructureParser().parse_mmcif(
        content, source_id=pdb_id, source="rcsb")


def resolve_representatives(structure):
    """Run the app's real resolution pipeline once per unique component.

    The parser deliberately leaves classification as None (it comes from the
    live CCD via LigandResolver, never invented in the parser), so tests that
    assert on classification must drive the same pipeline the app uses.
    """
    resolver = LigandResolver()
    resolved = {}
    for lig in structure.ligands:
        if lig.residue_name in resolved:
            continue
        resolved[lig.residue_name] = resolver.resolve_ligand(lig)
    return resolved


def drug_like_ligands(structure):
    """Ligands whose CCD classification is neither solvent nor ion."""
    return [l for l in structure.ligands
            if l.classification_hint not in ("HETAS", "HETAI")
            and l.classification_hint]


def test_1crn_protein_only_no_false_positive_ligand():
    """1CRN has no ligands at all - nothing may be picked as one."""
    s = fetch("1CRN")
    assert len(s.chains) >= 1
    assert all(c.is_polymer for c in s.chains) or s.ligands == []
    resolve_representatives(s)
    assert drug_like_ligands(s) == []
    # Every ligand entry (if any parsed as non-polymer) must be solvent-class.
    for lig in s.ligands:
        assert lig.classification_hint in ("HETAS", "HETAI", "non-polymer")


def test_1ubq_water_never_the_ligand():
    """1UBQ: protein + waters only; water must not be selected as ligand."""
    s = fetch("1UBQ")
    resolve_representatives(s)
    candidates = [l for l in s.ligands
                  if l.classification_hint not in ("HETAS", "HETAI")]
    assert candidates == [], (
        f"false positive: {[(l.residue_name, l.classification_hint) for l in candidates]}")
    # Waters themselves are classified as solvent by the CCD.
    waters = [l for l in s.ligands if l.classification_hint == "HETAS"]
    assert waters, "CCD must classify HOH as HETAS solvent"
    assert all(l.residue_name == "HOH" for l in waters)


def test_1fnb_fad_cofactor_is_found():
    """1FNB: ferredoxin-NADP+ reductase with a real FAD cofactor entity -
    no false negative. (1A4Y is the RNH-angiogenin complex and genuinely
    carries no FAD; 1FNB verified live to contain FAD entity 4.)"""
    s = fetch("1FNB")
    resolve_representatives(s)
    names = {l.residue_name for l in drug_like_ligands(s)}
    assert "FAD" in names, f"FAD cofactor missed; got {names}"
    fad = next(l for l in s.ligands if l.residue_name == "FAD")
    assert fad.classification_hint not in ("HETAS", "HETAI")


def test_3w85_multiple_ligand_types_all_classified():
    """3W85: inhibitor + FMN cofactor + sulfate + glycerol + waters."""
    s = fetch("3W85")
    resolve_representatives(s)
    by_name = {}
    for lig in s.ligands:
        by_name.setdefault(lig.residue_name, []).append(lig)
    assert "W85" in by_name and "FMN" in by_name
    assert by_name["HOH"][0].classification_hint == "HETAS"
    # The inhibitor and cofactor are never solvent.
    assert by_name["W85"][0].classification_hint not in ("HETAS", "HETAI")
    assert by_name["FMN"][0].classification_hint not in ("HETAS", "HETAI")


def test_1ehz_trna_rna_system_parses():
    """1EHZ: the 2.2 A tRNA structure - RNA + ions + spermine + water."""
    s = fetch("1EHZ")
    assert s.chains, "chains must parse for an all-RNA entry"
    resolve_representatives(s)
    # Waters must be present and classified as solvent (huge solvent arena).
    waters = [l for l in s.ligands if l.residue_name == "HOH"]
    assert waters, "1EHZ must contain waters"
    assert waters[0].classification_hint == "HETAS"
    # Spermine (SPM) is a real polyamine ligand, not solvent.
    names = {l.residue_name for l in s.ligands}
    if "SPM" in names:
        spm = next(l for l in s.ligands if l.residue_name == "SPM")
        assert spm.classification_hint not in ("HETAS", "HETAI")


def test_2m6q_nmr_entry_no_resolution_field():
    """2M6Q: NMR ensemble - parses without X-ray resolution and no crash."""
    s = fetch("2M6Q")
    assert s.chains, "NMR entry must still parse chains"
    assert s.resolution is None or s.resolution == 0 or s.resolution != s.resolution


def test_1mbn_heme_found_not_water():
    """1MBN: myoglobin with heme (HEM) + sulfate; heme must not be water."""
    s = fetch("1MBN")
    resolve_representatives(s)
    names = {l.residue_name for l in drug_like_ligands(s)}
    assert "HEM" in names, f"heme missed; ligands: {sorted(names)[:10]}"
    hem = next(l for l in s.ligands if l.residue_name == "HEM")
    assert hem.classification_hint not in ("HETAS", "HETAI")


def test_water_network_on_protein_only_structure():
    """Water analysis on 1UBQ (has waters): counts must be > 0 and honest."""
    s = fetch("1UBQ")
    analyzer = WaterNetworkAnalyzer(
        water_names={"HOH"}, network_cutoff=3.5, contact_cutoff=3.5)
    result = analyzer.analyze(s)
    # 1UBQ has ~100+ waters; if classification failed this is 0.
    assert result["water_count"] > 20, (
        f"false negative: only {result['water_count']} waters detected")
    assert len(result["water_contacts"]) > 0


def test_water_network_on_structure_without_water():
    """A structure with no waters reports 0 honestly (no invention)."""
    # 3TGX: uploaded without waters? Safer: strip - use 1CRN (no waters at all).
    s = fetch("1CRN")
    analyzer = WaterNetworkAnalyzer(
        water_names={"HOH"}, network_cutoff=3.5, contact_cutoff=3.5)
    result = analyzer.analyze(s)
    assert result["water_count"] == 0
    assert result["network"] == []
    assert result["water_contacts"] == []


def test_ligand_resolution_real_cids_across_types():
    """Resolution returns real PubChem CIDs for known ligands (no fallback)."""
    cases = {
        "3W85": ("W85", 73167555),
    }
    for pdb_id, (comp, expected_cid) in cases.items():
        s = fetch(pdb_id)
        lig = next(l for l in s.ligands if l.residue_name == comp)
        resolved = LigandResolver().resolve_ligand(lig, pdb_id)
        assert resolved.pubchem_cid == expected_cid, (
            f"{comp}: got {resolved.pubchem_cid}, expected {expected_cid}")
