"""
Protein-ligand interaction analysis using PLIP.

PLIP (Protein-Ligand Interaction Profiler) is the analysis engine: it is
executed on a PDB rendering of the loaded complex and its XML report is
parsed into contact records. All contact classification comes from PLIP.
When PLIP is not installed, the app reports that contact analysis is
unavailable - no substitute classifier is invented locally.

Geometric helpers in this module (binding pocket, ligand center) are pure
coordinate math on the user's structure, not chemical classification.
"""

import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

import numpy as np

from .schemas import (
    Structure,
    Ligand,
    Contact,
    ContactType,
)
from .config import get_config

# PLIP XML interaction *section* tag -> contact type. Section names are
# PLIP's own output schema (plural container tags); this is format mapping,
# not chemical classification.
_PLIP_TAG_MAP = {
    "hydrogen_bonds": ContactType.HYDROGEN_BOND,
    "hydrophobic_interactions": ContactType.HYDROPHOBIC,
    "pi_stacks": ContactType.PI_STACKING,
    "pi_stacking": ContactType.PI_STACKING,
    "salt_bridges": ContactType.SALT_BRIDGE,
    "halogen_bonds": ContactType.HALOGEN_BOND,
    "metal_complexes": ContactType.METAL_COORDINATION,
    "water_bridges": ContactType.WATER_MEDIATED,
    "pi_cation_interactions": ContactType.PI_STACKING,
    "cation_pi_interactions": ContactType.PI_STACKING,
    # singular variants tolerated for robustness across PLIP versions
    "hydrogen_bond": ContactType.HYDROGEN_BOND,
    "hydrophobic_interaction": ContactType.HYDROPHOBIC,
    "pi_stack": ContactType.PI_STACKING,
    "salt_bridge": ContactType.SALT_BRIDGE,
    "halogen_bond": ContactType.HALOGEN_BOND,
    "metal_complex": ContactType.METAL_COORDINATION,
    "water_bridge": ContactType.WATER_MEDIATED,
    "pi_cation_interaction": ContactType.PI_STACKING,
}


class ContactAnalyzer:
    """
    Analyze protein-ligand contacts by running PLIP and parsing its report.
    """

    def __init__(self):
        self._plip_available: Optional[bool] = None

    # ------------------------------------------------------------------
    # Availability
    # ------------------------------------------------------------------

    def _plip_command(self) -> Optional[str]:
        """Resolve the PLIP executable (explicit config or PATH)."""
        config = get_config()
        if config.plip_executable:
            return config.plip_executable
        return shutil.which("plip")

    def is_plip_available(self) -> bool:
        """Check PLIP availability by actually locating the executable."""
        if self._plip_available is None:
            self._plip_available = self._plip_command() is not None
        return self._plip_available

    # ------------------------------------------------------------------
    # Main entry
    # ------------------------------------------------------------------

    def analyze_contacts(
        self,
        structure: Structure,
        ligand: Ligand,
        protein_chains: Optional[List[str]] = None,
    ) -> Tuple[List[Contact], bool]:
        """
        Analyze contacts between a ligand and the protein using PLIP.

        Returns:
            (contacts, plip_available). When PLIP is unavailable the list is
            empty and plip_available is False; the caller must surface that
            to the user rather than presenting 'no contacts'.
        """
        if not self.is_plip_available():
            return [], False
        if not ligand.atoms:
            return [], True

        contacts = self._run_plip(structure, ligand, protein_chains)
        for i, c in enumerate(contacts):
            c.id = i
        return contacts, True

    # ------------------------------------------------------------------
    # PLIP execution
    # ------------------------------------------------------------------

    def _run_plip(self, structure: Structure, ligand: Ligand,
                  protein_chains: Optional[List[str]] = None) -> List[Contact]:
        """Write the complex to PDB, run PLIP, parse its XML report."""
        config = get_config()
        with tempfile.TemporaryDirectory(prefix="ligora_plip_") as tmpdir:
            pdb_path = Path(tmpdir) / "complex.pdb"
            outdir = Path(tmpdir) / "out"
            outdir.mkdir()
            self._write_pdb_for_plip(structure, ligand, pdb_path,
                                     protein_chains)

            cmd = [
                self._plip_command(),
                "-f", str(pdb_path),
                "-x",                 # XML report
                "-o", str(outdir),
                "-q",                 # quiet
                "--name", "report",
            ]
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=config.plip_timeout,
                    cwd=str(tmpdir),
                )
            except (subprocess.TimeoutExpired, OSError):
                return []

            if result.returncode != 0:
                return []

            xml_path = outdir / "report_report.xml"
            if not xml_path.exists():
                # PLIP names output as <name>_<inputstem>.xml; find any xml.
                candidates = sorted(outdir.glob("*.xml"))
                if not candidates:
                    return []
                xml_path = candidates[0]

            # PLIP's *idx fields are PDB serial numbers from the file we
            # wrote; build the serial -> atom-name map from that same file.
            # Sections without serials (salt bridges, pi interactions) carry
            # group coordinates instead, resolved against the same file.
            serial_names = self._serial_atom_names(pdb_path)
            return self._parse_plip_xml(xml_path, serial_names, pdb_path)

    @staticmethod
    def _serial_atom_names(pdb_path: Path) -> Dict[int, str]:
        """Map PDB serial numbers to atom names from the complex file."""
        mapping: Dict[int, str] = {}
        try:
            for line in pdb_path.read_text().splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        serial = int(line[6:11])
                    except ValueError:
                        continue
                    mapping[serial] = line[12:16].strip()
        except OSError:
            pass
        return mapping

    @staticmethod
    def _coord_to_atom(coord_text: Optional[str], resname: str,
                       resid: int, chain: str, pdb_path: Path
                       ) -> Optional[str]:
        """Resolve a PLIP coordinate field to an atom name.

        PLIP reports group coordinates (e.g. carboxylate center) rather than
        atom serials for some sections (salt bridges, pi interactions). The
        real participating atom is recovered by finding the nearest atom of
        the interacting residue to the reported coordinate - geometry from
        the file, not a heuristic classification.
        """
        if not coord_text:
            return None
        parts = coord_text.replace("(", " ").replace(")", " ").split()
        try:
            target = np.array([float(p) for p in parts[:3]])
        except (ValueError, TypeError):
            return None
        best_name = None
        best_dist = None
        try:
            for line in pdb_path.read_text().splitlines():
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                if line[17:20].strip() != resname or \
                        line[21].strip() != chain:
                    continue
                try:
                    if int(line[22:26]) != resid:
                        continue
                except ValueError:
                    continue
                pos = np.array([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])])
                d = float(np.linalg.norm(pos - target))
                if best_dist is None or d < best_dist:
                    best_dist = d
                    best_name = line[12:16].strip()
        except OSError:
            return None
        return best_name

    def _parse_plip_xml(self, xml_path: Path,
                        serial_names: Optional[Dict[int, str]] = None,
                        pdb_path: Optional[Path] = None,
                        ) -> List[Contact]:
        """Parse a PLIP XML report into Contact records."""
        contacts: List[Contact] = []
        try:
            tree = ET.parse(str(xml_path))
        except ET.ParseError:
            return []

        for bindingsite in tree.getroot().findall("bindingsite"):
            interactions = bindingsite.find("interactions")
            if interactions is None:
                continue
            for section in interactions:
                tag = section.tag.lower()
                contact_type = _PLIP_TAG_MAP.get(tag)
                if contact_type is None:
                    continue
                for interaction in section:
                    contact = self._parse_interaction(
                        interaction, contact_type, section.tag,
                        serial_names or {}, pdb_path)
                    if contact:
                        contacts.append(contact)
        return contacts

    def _parse_interaction(self, elem, contact_type: ContactType,
                           section_name: str,
                           serial_names: Dict[int, str],
                           pdb_path: Optional[Path] = None,
                           ) -> Optional[Contact]:
        """Parse one PLIP interaction element into a Contact."""

        def _text(tag: str) -> Optional[str]:
            node = elem.find(tag)
            return node.text.strip() if node is not None and node.text else None

        resnr = _text("resnr")
        restype = _text("restype")
        reschain = _text("reschain")
        resnr_lig = _text("resnr_lig")
        restype_lig = _text("restype_lig")
        reschain_lig = _text("reschain_lig")
        if not restype or not resnr:
            return None

        # Distances per PLIP section type (fields from PLIP's report schema).
        distance = None
        angle = None
        for dist_tag in ("dist_h-a", "dist_d-a", "distance", "dist",
                         "mindist", "dist_dist"):
            val = _text(dist_tag)
            if val:
                try:
                    distance = float(val)
                    break
                except ValueError:
                    continue
        for angle_tag in ("don_angle", "angle", "ang"):
            val = _text(angle_tag)
            if val:
                try:
                    angle = float(val)
                    break
                except ValueError:
                    continue

        # Atom identity: PLIP reports the participating atoms as PDB serial
        # numbers in *idx fields (field names differ per section). Resolve
        # them to the real atom names from the file we wrote.
        def _first_serial(tag: str) -> Optional[int]:
            val = _text(tag)
            if not val:
                return None
            try:
                return int(val.split(",")[0].strip())
            except ValueError:
                return None

        protisdon = (_text("protisdon") or "").lower() == "true"
        lig_serial: Optional[int] = None
        prot_serial: Optional[int] = None
        if section_name == "hydrophobic_interactions":
            lig_serial = _first_serial("ligcarbonidx")
            prot_serial = _first_serial("protcarbonidx")
        elif section_name == "hydrogen_bonds":
            donor = _first_serial("donoridx")
            acceptor = _first_serial("acceptoridx")
            if protisdon:
                prot_serial, lig_serial = donor, acceptor
            else:
                lig_serial, prot_serial = donor, acceptor
        elif section_name == "salt_bridges":
            lig_serial = _first_serial("lig_idx_list")
            prot_serial = _first_serial("prot_idx_list")
        elif section_name == "water_bridges":
            donor = _first_serial("donoridx")
            acceptor = _first_serial("acceptoridx")
            if protisdon:
                prot_serial, lig_serial = donor, acceptor
            else:
                lig_serial, prot_serial = donor, acceptor
        elif section_name == "metal_complexes":
            # metalidx is the metal; the counterpart atom comes from the
            # target side reported by PLIP for that section.
            metal = _first_serial("metalidx")
            target = _first_serial("targetidx")
            lig_serial, prot_serial = metal, target

        def _coord_text(tag: str) -> Optional[str]:
            """Read a PLIP coordinate element (<x>/<y>/<z> children)."""
            node = elem.find(tag)
            if node is None:
                return None
            vals = []
            for axis in ('x', 'y', 'z'):
                child = node.find(axis)
                if child is None or not child.text:
                    return None
                vals.append(child.text.strip())
            return ' '.join(vals) if len(vals) == 3 else None

        def _name(serial: Optional[int], fallback_idx, coord_tag: str,
                  resname: str, resid: int, chain: str) -> str:
            if serial is not None and serial in serial_names:
                return serial_names[serial]
            if pdb_path is not None:
                resolved = self._coord_to_atom(
                    _coord_text(coord_tag), resname, resid, chain, pdb_path)
                if resolved:
                    return resolved
            return f"idx:{fallback_idx or '?'}"

        lig_atom_name = _name(
            lig_serial, _text("ligidx") or _text("ligcarbonidx"),
            "ligcoo", restype_lig or "", int(resnr_lig) if resnr_lig else 0,
            reschain_lig or "")
        prot_atom_name = _name(
            prot_serial, _text("protidx") or _text("protcarbonidx"),
            "protcoo", restype, int(resnr), reschain or "")

        description = (
            f"{section_name} between {restype}{resnr}({reschain}) and "
            f"{restype_lig}{resnr_lig}({reschain_lig})"
        )

        is_water = bool(_text("wateridx")) or restype in ("HOH", "WAT")

        return Contact(
            id=0,
            ligand_atom=lig_atom_name,
            ligand_residue_name=restype_lig or "",
            ligand_residue_id=int(resnr_lig) if resnr_lig else 0,
            ligand_chain_id=reschain_lig or "",
            protein_residue_name=restype,
            protein_residue_id=int(resnr),
            protein_chain_id=reschain or "",
            protein_atom=prot_atom_name,
            distance=distance if distance is not None else 0.0,
            contact_type=contact_type,
            angle=angle,
            is_water_mediated=is_water,
            description=description,
        )

    # ------------------------------------------------------------------
    # PDB rendering for PLIP
    # ------------------------------------------------------------------

    def _write_pdb_for_plip(self, structure: Structure, ligand: Ligand,
                            pdb_path: Path,
                            protein_chains: Optional[List[str]] = None):
        """Write the complex (protein + selected ligand) in PDB format.

        Only ATOM/HETATM records: strict parsers must not encounter
        unexpected record types in interaction input files.
        """
        lines: List[str] = []
        serial = 1

        def fmt_atom(record, atom, resname, chain, resid):
            # PDB fixed columns: record 1-6, serial 7-11, name 13-16,
            # altLoc 17, resName 18-20, chainID 22, resSeq 23-26,
            # x/y/z 31-54, occupancy 55-60, tempFactor 61-66,
            # element 77-78.
            element = (atom.element or "").upper()[:2].rjust(2)
            return (
                f"{record:<6}{serial:>5}  "
                f"{atom.name:<4}"
                f"{resname:<3} "
                f"{chain:<1}"
                f"{resid:>4}    "
                f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}"
                f"{atom.occupancy:>6.2f}{atom.b_factor:>6.2f}"
                f"{'':>10}"
                f"{element}"
            )

        for chain in structure.chains:
            if not chain.is_polymer:
                continue
            if protein_chains and chain.id not in protein_chains:
                continue
            for residue in chain.residues:
                for atom in residue.atoms:
                    lines.append(fmt_atom("ATOM", atom, residue.name,
                                          chain.id, residue.id))
                    serial += 1

        for atom in ligand.atoms:
            lines.append(fmt_atom("HETATM", atom, ligand.residue_name,
                                  atom.chain_id or "L",
                                  atom.residue_id or 1))
            serial += 1

        lines.append("END")
        pdb_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ------------------------------------------------------------------
    # Geometric helpers (coordinate math on user data, no classification)
    # ------------------------------------------------------------------

    def compute_binding_pocket(
        self,
        structure: Structure,
        ligand: Ligand,
        radius: float = None,
    ) -> Dict[str, Any]:
        """
        Compute the binding pocket around a ligand by geometric proximity.

        radius is a user setting (default from config, 6 A); it selects which
        residues are near the ligand and implies nothing about interaction
        types.
        """
        if radius is None:
            radius = 6.0

        ligand_center = self._compute_ligand_center(ligand)

        pocket_residues = []
        pocket_chains = set()

        for chain in structure.chains:
            if not chain.is_polymer:
                continue
            for residue in chain.residues:
                residue_center = self._compute_residue_center(residue)
                if residue_center is None:
                    continue
                dist = float(np.linalg.norm(ligand_center - residue_center))
                if dist <= radius:
                    pocket_residues.append({
                        "chain_id": chain.id,
                        "residue_name": residue.name,
                        "residue_id": residue.id,
                        "distance": round(dist, 2),
                    })
                    pocket_chains.add(chain.id)

        pocket_residues.sort(key=lambda r: r["distance"])
        return {
            "ligand_center": {
                "x": round(float(ligand_center[0]), 3),
                "y": round(float(ligand_center[1]), 3),
                "z": round(float(ligand_center[2]), 3),
            },
            "radius": radius,
            "pocket_chain_ids": sorted(pocket_chains),
            "pocket_residue_count": len(pocket_residues),
            "pocket_residues": pocket_residues,
        }

    def _compute_ligand_center(self, ligand: Ligand) -> np.ndarray:
        """Compute the geometric center of a ligand's atoms."""
        if not ligand.atoms:
            raise ValueError("Ligand has no atoms; center is undefined")
        positions = np.array([[a.x, a.y, a.z] for a in ligand.atoms])
        return np.mean(positions, axis=0)

    def _compute_residue_center(self, residue) -> Optional[np.ndarray]:
        """Compute the geometric center of a residue's atoms."""
        if not residue.atoms:
            return None
        positions = np.array([[a.x, a.y, a.z] for a in residue.atoms])
        return np.mean(positions, axis=0)

    def export_contacts_csv(self, contacts: List[Contact],
                            output_path: Path):
        """Export contacts to a CSV file."""
        import csv
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([
                "id", "ligand_atom", "ligand_res", "ligand_res_id",
                "ligand_chain", "protein_atom", "protein_res",
                "protein_res_id", "protein_chain", "distance",
                "contact_type", "angle", "water_mediated", "description",
            ])
            for c in contacts:
                writer.writerow([
                    c.id, c.ligand_atom, c.ligand_residue_name,
                    c.ligand_residue_id, c.ligand_chain_id,
                    c.protein_atom, c.protein_residue_name,
                    c.protein_residue_id, c.protein_chain_id,
                    round(c.distance, 3), c.contact_type.value,
                    c.angle if c.angle is not None else "",
                    c.is_water_mediated, c.description,
                ])
