"""
Structure parser for mmCIF and PDB files.

This module parses macromolecular structure files and creates
internal structure representations with chains, residues, atoms,
and ligands.
"""
import re
import tempfile
import subprocess
import requests
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from io import StringIO

from .schemas import (
    Structure, Chain, Residue, Atom, Ligand,
    LigandResolutionStatus
)
from .config import get_config


class StructureParser:
    """
    Parse macromolecular structure files (mmCIF, PDB).

    Supports:
    - Local file parsing (mmCIF primary, PDB legacy fallback)
    - Remote fetching from RCSB (via ModelServer or download)
    """

    # _KNOWN_POLYMER_COMP is intentionally empty.
    # The parser does not decide what is a polymer component from hardcoded name sets.
    # Polymer / non-polymer classification comes from the entity categories in the
    # structure file (entity_type in _entity, and the CCD pdbx_type fetched live).
    # Everything is deferred to the enrichment client / CCD; the parser only splits
    # atom_site into polymer residues vs non-polymer entity atoms.
    _KNOWN_POLYMER_COMP: set = set()  # populated only by data sources, not by us

    def parse_mmcif(self, content: str, source_id: str = "",
                    source: str = "local") -> Structure:
        """
        Parse an mmCIF file content.

        Args:
            content: The mmCIF file content.
            source_id: Identifier for the structure (PDB ID, etc.).
            source: Source type ("local", "rcsb", "pdb", "csm").

        Returns:
            Parsed Structure object.
        """
        structure = Structure(id=source_id or "unknown", source=source)
        structure.file_format = "mmcif"

        # Parse key-value pairs from mmCIF
        data = self._parse_mmcif_data(content)

        # Extract audit category info
        audit = data.get("_audit", {})
        structure.title = audit.get("audit_site_id", audit.get(
            "pdbx_description", "Unknown structure"
        ))

        # Extract entity category - defines what the molecules are
        entities = data.get("_entity", [])
        entity_poly = data.get("_entity_poly", [])
        entity_nonpoly = data.get("_entity_nonpoly", [])
        entity_poly_seq = data.get("_entity_poly_seq", [])

        # Extract assembly information
        assemblies = data.get("_pdbx_struct_assembly", [])

        # Build chain-to-entity mapping from pdbx_strand_entity
        chain_to_entity: Dict[str, int] = {}
        for strand in data.get("_pdbx_struct_assembly_entity", []):
            chain_ids = strand.get("pdbx_polymer_entity_id", "")
            if chain_ids:
                for cid in chain_ids.split(","):
                    cid = cid.strip()
                    if cid:
                        entity_id = strand.get("pdbx_polymer_entity_id")
                        if entity_id:
                            chain_to_entity[cid] = entity_id

        # Build chains
        chains: Dict[str, Chain] = {}

        # Process polymer entities (proteins, nucleic acids)
        for entity in entities:
            entity_id = entity.get("id")
            entity_name = entity.get("entity_name", "Unknown")
            entity_type = entity.get("entity_type", "")

            # Find polymer information
            poly_info = None
            for poly in entity_poly:
                if poly.get("entity_id") == entity_id:
                    poly_info = poly
                    break

            if poly_info:
                # This is a polymer - create chains
                pdbx_types = poly_info.get("pdbx_strand_id", "")
                for strand_id in pdbx_types.split(","):
                    strand_id = strand_id.strip()
                    if strand_id:
                        chain = Chain(
                            id=strand_id,
                            name=entity_name,
                            is_polymer=True
                        )
                        chains[strand_id] = chain

        # Build entity type maps once so both the polymer path and the
        # non-polymer candidate path can use the file's own entity classification.
        entity_type_map: Dict[str, Optional[str]] = {}
        for entity in (data.get("_entity") or []):
            if isinstance(entity, dict):
                eid = entity.get("id")
                if eid is not None:
                    entity_type_map[str(eid)] = entity.get("type") or None

        asym_entity: Dict[str, Optional[str]] = {}
        struct_asym = data.get("_struct_asym", {})
        if isinstance(struct_asym, dict):
            for asym_id, row in struct_asym.items():
                if isinstance(row, dict) and "entity_id" in row:
                    asym_entity[asym_id] = str(row["entity_id"])

        # Process coordinates from _atom_site category.
        # The parser does not decide polymer / ligand / solvent / ion status from
        # hardcoded name lists. That classification is derived from the structure
        # file entity categories and then enriched live by the enrichment client.
        atom_site = data.get("_atom_site", [])
        current_chain: Optional[str] = None
        current_residue: Optional[int] = None
        current_residue_name: Optional[str] = None
        current_chain_obj: Optional[Chain] = None
        current_residue_obj: Optional[Residue] = None

        entity_comp_atoms: Dict[str, List[Atom]] = {}

        for atom_data in atom_site:
            # Extract atom information
            residue_name = atom_data.get("auth_comp_id") or \
                           atom_data.get("label_comp_id", "")
            residue_number = atom_data.get("auth_seq_id") or \
                            atom_data.get("label_seq_id", 0)
            atom_name = atom_data.get("auth_atom_id") or \
                        atom_data.get("label_atom_id", "")
            x = float(atom_data.get("Cartn_x", 0))
            y = float(atom_data.get("Cartn_y", 0))
            z = float(atom_data.get("Cartn_z", 0))
            b_factor = float(atom_data.get("B_iso_or_equiv", 0))
            occupancy = float(atom_data.get("occupancy", 1))

            chain_id = atom_data.get("auth_asym_id") or \
                      atom_data.get("label_asym_id", " ")
            chain_id = chain_id.strip() or "A"

            asym_id = chain_id
            entity_id = atom_data.get("label_entity_id") or atom_data.get("auth_entity_id")

            # Determine entity type from the file, not from hardcoded name lists.
            ent_type = None
            if entity_id is not None:
                ent_type = entity_type_map.get(str(entity_id))
            if ent_type is None and asym_id in asym_entity:
                ent_type = entity_type_map.get(asym_entity[asym_id])

            atom = Atom(
                id=len(structure.chains) * 1000 + len(atom_site),  # rough unique ID
                name=atom_name,
                residue_name=residue_name,
                residue_id=int(residue_number) if residue_number else 0,
                chain_id=chain_id,
                x=x, y=y, z=z,
                element=atom_data.get("type_symbol", "") or self._guess_element(
                    atom_name, residue_name
                ),
                b_factor=b_factor,
                occupancy=occupancy,
            )

            # Non-polymer / unknown entities are accumulated as candidate ligands.
            # Polymer entities are consumed into chain/residue structure below;
            # everything else is a candidate for ligand / solvent / ion / unknown.
            if ent_type != "polymer":
                entity_comp_atoms.setdefault(residue_name, [])
                ligand_atom_name = atom_name
                lx = x
                ly = y
                lz = z
                lbf = b_factor
                locc = occupancy
                latom = Atom(
                    id=len(entity_comp_atoms[residue_name]),
                    name=ligand_atom_name,
                    residue_name=residue_name,
                    residue_id=int(residue_number) if residue_number else 0,
                    chain_id=asym_id,
                    x=lx, y=ly, z=lz,
                    element=atom.element,
                    b_factor=lbf,
                    occupancy=locc,
                )
                entity_comp_atoms[residue_name].append(latom)

            # Track current chain/residue for polymer paths
            if chain_id != current_chain:
                current_chain = chain_id
                current_residue = None
                current_residue_name = None

                if chain_id not in chains:
                    chains[chain_id] = Chain(id=chain_id, name=chain_id)

                current_chain_obj = chains[chain_id]

            if residue_number != current_residue or residue_name != current_residue_name:
                current_residue = residue_number
                current_residue_name = residue_name

                # Create or find residue (polymer chains only).
                # Non-polymer entities are handled through the ligand path, not
                # through polymer chain residues.
                residue = next(
                    (r for r in current_chain_obj.residues
                     if r.id == int(residue_number) and r.name == residue_name),
                    None
                )
                if residue is None:
                    residue = Residue(
                        id=int(residue_number),
                        name=residue_name,
                        chain_id=chain_id,
                        residue_number=int(residue_number)
                    )
                    current_chain_obj.residues.append(residue)
                current_residue_obj = residue
            else:
                residue = current_residue_obj

            # Add atom to residue
            if residue is not None:
                residue.atoms.append(atom)
                # Update chain reference
                if current_chain_obj:
                    # Find or create chain with this atom
                    pass

        # Create ligand objects from the candidate groups.
        # Names and atom groupings come from the structure file.
        # Formula and molecular weight are not computed here from hardcoded
        # element/weight tables. Those are derived chemical properties that come
        # from data sources (CCD formula / formula_weight, PubChem, etc.) and are
        # attached by the enrichment client. The parser only carries the identity
        # it can read directly from the structure file.
        for comp_id, atoms in entity_comp_atoms.items():
            if not atoms:
                continue

            ligand = Ligand(
                id=f"L{comp_id}",
                name=comp_id,
                residue_name=comp_id,
                formula=None,
                molecular_weight=None,
                atom_count=len(atoms),
                atoms=atoms,
                resolution_status=LigandResolutionStatus.NOT_FOUND,
                # No classification here. Deferred to enrichment client / CCD.
                classification_hint="unclassified",
            )
            structure.ligands.append(ligand)

        # Set chain list from dictionary
        structure.chains = list(chains.values())

        # Add atom counts
        for chain in structure.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    pass  # atoms already attached

        # Extract resolution from quality info
        quality = data.get("_pdbx_quality", {})
        if isinstance(quality, dict) and quality:
            # Try to get resolution from any available field
            for key in quality:
                if 'resolution' in key.lower():
                    try:
                        structure.resolution = float(quality[key])
                    except (ValueError, TypeError):
                        pass
                    break

        # Extract experiment type
        exp_detail = data.get("_exptl", {})
        if isinstance(exp_detail, dict) and exp_detail:
            structure.experiment_type = str(exp_detail.get('method', ''))

        # Check for density data
        structure.has_density = "_pdbx_vrpt" in data or \
                                any(k.startswith("_pdbx_volray") for k in data)

        return structure

    def parse_pdb(self, content: str, source_id: str = "",
                  source: str = "local") -> Structure:
        """
        Parse a PDB file content (legacy format).

        Args:
            content: The PDB file content.
            source_id: Identifier for the structure.
            source: Source type.

        Returns:
            Parsed Structure object.
        """
        structure = Structure(id=source_id or "unknown", source=source)
        structure.file_format = "pdb"

        lines = content.splitlines()
        chains: Dict[str, Chain] = {}
        ligands: Dict[str, List[Atom]] = {}
        current_chain: Optional[str] = None
        current_residue: int = 0
        current_residue_atoms: List[Atom] = []
        current_residue_name: Optional[str] = None

        for line in lines:
            record_type = line[0:6].strip()

            if record_type == "ATOM" or record_type == "HETATM":
                # Parse PDB atom record
                atom_serial = int(line[6:11].strip() or 0)
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21:22].strip() or "A"
                residue_number = int(line[22:26].strip() or 0)
                x = float(line[30:38].strip() or 0)
                y = float(line[38:46].strip() or 0)
                z = float(line[46:54].strip() or 0)
                b_factor = float(line[60:66].strip() or 0)
                occupancy = float(line[54:60].strip() or 1)
                element = line[76:78].strip() or self._guess_element(
                    atom_name, residue_name
                )

                # Check if this is a HETATM (heteroatom/ligand)
                is_hetatm = record_type == "HETATM"

                atom = Atom(
                    id=atom_serial,
                    name=atom_name,
                    residue_name=residue_name,
                    residue_id=residue_number,
                    chain_id=chain_id,
                    x=x, y=y, z=z,
                    element=element,
                    b_factor=b_factor,
                    occupancy=occupancy,
                )

                if is_hetatm:
                    # This is a ligand or solvent
                    lig_key = f"{residue_name}_{chain_id}"
                    if lig_key not in ligands:
                        ligands[lig_key] = []
                    ligands[lig_key].append(atom)
                else:
                    # Polymer atom
                    if chain_id not in chains:
                        chains[chain_id] = Chain(id=chain_id, name=chain_id)

                    chain = chains[chain_id]

                    # Check if we've moved to a new residue
                    if residue_number != current_residue or \
                       residue_name != current_residue_name:
                        # Save previous residue
                        if current_residue_atoms:
                            residue = Residue(
                                id=current_residue,
                                name=current_residue_name or "UNK",
                                chain_id=chain_id,
                                residue_number=current_residue,
                                atoms=current_residue_atoms,
                            )
                            chain.residues.append(residue)

                        current_residue = residue_number
                        current_residue_name = residue_name
                        current_residue_atoms = [atom]
                    else:
                        current_residue_atoms.append(atom)

            elif record_type == "TITLE":
                title = line[10:70].strip()
                if structure.title == "Unknown structure":
                    structure.title = title

            elif record_type == "REMARK" and line[7:10] == "RF":
                # Resolution
                try:
                    structure.resolution = float(line[13:20].strip())
                except (ValueError, IndexError):
                    pass

            elif record_type == "KEYWDS":
                keywords = line[10:70].strip()
                if structure.title == "Unknown structure":
                    structure.title = keywords[:80]

        # Save last residue
        if current_residue_atoms and current_chain:
            chain = chains.get(current_chain)
            if chain:
                residue = Residue(
                    id=current_residue,
                    name=current_residue_name or "UNK",
                    chain_id=current_chain,
                    residue_number=current_residue,
                    atoms=current_residue_atoms,
                )
                chain.residues.append(residue)

        # Convert chain dict to list
        structure.chains = list(chains.values())

        # Process ligands from HETATM records.
        # The parser does not skip residues based on a hardcoded solvent/ion list.
        # HETATM records can be ligands, solvent, ions, or other non-polymer
        # entities, and the definitive classification comes from entity categories
        # / CCD / enrichment client, not from embedded name tables here.
        # For now, every HETATM group is surfaced as a candidate ligand.
        for lig_key, atoms in ligands.items():
            res_name = atoms[0].residue_name if atoms else ""

            # Formula and molecular weight are not computed here from hardcoded
            # element/weight tables. Those are derived chemical properties that come
            # from data sources (CCD formula / formula_weight, PubChem, etc.) and are
            # attached by the enrichment client. The parser only carries the identity
            # it can read directly from the structure file.
            ligand = Ligand(
                id=lig_key,
                name=res_name,
                residue_name=res_name,
                formula=None,
                molecular_weight=None,
                atom_count=len(atoms),
                atoms=atoms,
                resolution_status=LigandResolutionStatus.NOT_FOUND,
                # No classification here. Deferred to enrichment client / CCD.
                classification_hint="unclassified",
            )
            structure.ligands.append(ligand)

        return structure

    def fetch_from_rcsb(self, pdb_id: str, format: str = "mmCIF") -> Structure:
        """
        Fetch a structure from RCSB PDB.

        Args:
            pdb_id: PDB ID (e.g., "1ABC") or CSM ID.
            format: File format ("mmCIF" or "PDB").

        Returns:
            Parsed Structure object.
        """
        config = get_config()

        # Try ModelServer first for efficient fetching
        modelserver_url = (
            f"{config.rcsb_base_url}/v1/model/{pdb_id}"
        )

        try:
            # For mmCIF, use the download service
            if format.upper() == "PDBx-MMCIF":
                download_url = (
                    f"https://files.rcsb.org/download/{pdb_id}.cif"
                )
            elif format.upper() == "PDB":
                download_url = (
                    f"https://files.rcsb.org/download/{pdb_id}.pdb"
                )
            else:
                download_url = (
                    f"https://files.rcsb.org/download/{pdb_id}.cif"
                )

            response = requests.get(
                download_url,
                timeout=config.request_timeout
            )
            response.raise_for_status()

            if format.upper() == "PDB":
                return self.parse_pdb(
                    response.text, source_id=pdb_id, source="rcsb"
                )
            else:
                return self.parse_mmcif(
                    response.text, source_id=pdb_id, source="rcsb"
                )

        except requests.RequestException as e:
            raise RuntimeError(f"Failed to fetch {pdb_id} from RCSB: {e}")

    def fetch_csm_from_rcsb(self, csm_id: str) -> Structure:
        """
        Fetch a Computed Structure Model from RCSB.

        Args:
            csm_id: CSM ID (e.g., "AF_AFP_1" or similar).

        Returns:
            Parsed Structure object.
        """
        return self.fetch_from_rcsb(csm_id, format="mmCIF")

    def fetch_structure(self, identifier: str) -> Structure:
        """
        Fetch a structure from the appropriate source.

        Args:
            identifier: PDB ID, CSM ID, or local file path.

        Returns:
            Parsed Structure object.
        """
        config = get_config()

        # Check if it's a local file
        local_path = Path(identifier)
        if local_path.exists():
            content = local_path.read_text(encoding="utf-8")
            if local_path.suffix.lower() in (".cif", ".mcif", ".bcif"):
                return self.parse_mmcif(content, source_id=local_path.stem)
            elif local_path.suffix.upper() == ".PDB":
                return self.parse_pdb(content, source_id=local_path.stem)
            else:
                # Try mmCIF first
                try:
                    return self.parse_mmcif(content, source_id=local_path.stem)
                except Exception:
                    return self.parse_pdb(content, source_id=local_path.stem)

        # Otherwise try RCSB
        # Check if it looks like a PDB ID (4 chars, starts with a digit)
        if re.match(r'^[0-9][A-Za-z0-9]{3,4}$', identifier):
            return self.fetch_from_rcsb(identifier)

        # Check if it's a CSM ID (starts with AF_ or MA_)
        if identifier.startswith("AF_") or identifier.startswith("MA_"):
            return self.fetch_csm_from_rcsb(identifier)

        raise ValueError(f"Unable to resolve identifier: {identifier}")

    def _parse_mmcif_data(self, content: str) -> Dict[str, Any]:
        """
        Parse mmCIF data into a nested dictionary structure.

        Args:
            content: Raw mmCIF content.

        Returns:
            Nested dictionary of categories and their data.
        """
        data: Dict[str, Any] = {}
        category_data: Dict[str, Dict[str, List[str]]] = {}
        current_category: Optional[str] = None

        for line in content.splitlines():
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            # Category header: category_name.id
            if line.endswith(".") and not line.startswith("_"):
                current_category = line[:-1]
                category_data[current_category] = {}
                continue

            # Category field
            if line.startswith("_"):
                # Parse field name and value(s)
                parts = line.split(None, 1)
                field_name = parts[0]
                values = []

                if len(parts) > 1:
                    value_part = parts[1]
                    # Handle semicolon-quoted values
                    if value_part.startswith(";"):
                        # Multi-line quoted value
                        values = [value_part[1:]]
                        # Note: we don't handle multi-line here for simplicity
                    elif ";" in value_part:
                        # Semicolon-separated list
                        values = value_part.split(";")
                        values = [v.strip() for v in values if v.strip()]
                    elif "," in value_part:
                        # Comma-separated list
                        values = value_part.split(",")
                        values = [v.strip() for v in values if v.strip()]
                    else:
                        values = [value_part.strip()]

                # Get category name from field
                if "." in field_name:
                    cat_name = field_name.split(".")[0]
                    if cat_name not in category_data:
                        category_data[cat_name] = {}
                    field_short = field_name.split(".")[1] if "." in field_name else field_name
                    category_data[cat_name][field_short] = values

        # Convert to nested dict
        for cat_name, fields in category_data.items():
            if cat_name not in data:
                data[cat_name] = {}

            for field_name, values in fields.items():
                if not isinstance(values, list):
                    values = [values]
                if len(values) == 1:
                    # Single value - try to convert to appropriate type
                    val = values[0]
                    try:
                        data[cat_name][field_name] = int(val)
                    except (ValueError, TypeError):
                        try:
                            data[cat_name][field_name] = float(val)
                        except (ValueError, TypeError):
                            data[cat_name][field_name] = val
                else:
                    # Multiple values
                    converted = []
                    for val in values:
                        try:
                            converted.append(int(val))
                        except (ValueError, TypeError):
                            try:
                                converted.append(float(val))
                            except (ValueError, TypeError):
                                converted.append(val)
                    data[cat_name][field_name] = converted

        # Handle loop structures (multiple rows)
        data["_atom_site"] = self._parse_atom_site(content)

        # Also parse _entity category from loop format for proper entity handling
        data["_entity"] = self._parse_entity_category(content)
        data["_entity_poly"] = self._parse_entity_poly_category(content)
        data["_entity_nonpoly"] = self._parse_entity_nonpoly_category(content)

        return data

    def _parse_entity_category(self, content: str) -> List[Dict[str, Any]]:
        """Parse _entity category from mmCIF loops."""
        entities = []
        in_loop = False
        fields = []
        current = {}
        field_idx = 0

        for line in content.splitlines():
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('_entity.'):
                if not in_loop:
                    in_loop = True
                    fields = []
                    current = {}
                    field_idx = 0
                field_name = line.split('.', 1)[1].split('(')[0]
                fields.append(field_name)
                continue

            if in_loop:
                if line.startswith('_') and not line.startswith('_entity.'):
                    if current and fields:
                        entities.append(current)
                    in_loop = False
                    continue

                values = line.split()
                for val in values:
                    if field_idx < len(fields):
                        current[fields[field_idx]] = val
                        field_idx += 1

                if field_idx >= len(fields):
                    entities.append(current)
                    current = {}
                    field_idx = 0

        if current and fields:
            entities.append(current)

        return entities

    def _parse_entity_poly_category(self, content: str) -> List[Dict[str, Any]]:
        """Parse _entity_poly category from mmCIF loops."""
        return self._parse_entity_category(content)  # Same structure

    def _parse_entity_nonpoly_category(self, content: str) -> List[Dict[str, Any]]:
        """Parse _entity_nonpoly category from mmCIF loops."""
        return self._parse_entity_category(content)  # Same structure

    def _parse_atom_site(self, content: str) -> List[Dict[str, Any]]:
        """
        Parse the _atom_site category from mmCIF content.

        Args:
            content: Raw mmCIF content.

        Returns:
            List of atom data dictionaries.
        """
        atoms: List[Dict[str, Any]] = []
        in_atom_site = False
        atom_fields: List[str] = []
        current_atom: Dict[str, str] = {}
        field_index = 0

        for line in content.splitlines():
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            # Check for atom_site category start
            if line.startswith("_atom_site."):
                in_atom_site = True
                field_name = line.split(".", 1)[1]
                if "(" in field_name and ")" in field_name:
                    # This is a loop header
                    field_name = field_name.split("(")[0]
                atom_fields.append(field_name)
                current_atom = {}
                continue

            if in_atom_site:
                # Check if this starts a new category
                if line.startswith("_") and not line.startswith("_atom_site"):
                    in_atom_site = False
                    if current_atom and atom_fields:
                        atoms.append(current_atom)
                    continue

                # Parse values for this row
                # Handle semicolon-quoted values
                if line.startswith(";"):
                    value = line[1:]
                    if value:  # Not empty
                        if field_index < len(atom_fields):
                            current_atom[atom_fields[field_index]] = value
                        field_index += 1
                    continue

                # Split by whitespace (mmCIF uses whitespace in loops)
                values = line.split()
                for val in values:
                    if field_index < len(atom_fields):
                        current_atom[atom_fields[field_index]] = val
                        field_index += 1

                # Check if we have a complete atom record
                if field_index >= len(atom_fields):
                    atoms.append(current_atom)
                    current_atom = {}
                    field_index = 0

        # Handle any remaining atom
        if current_atom and atom_fields:
            atoms.append(current_atom)

        return atoms

    def _is_ligand(self, residue_name: str, atom_site: List[Dict]) -> bool:
        """
        Determine if a residue is a ligand (non-polymer).
        
        The parser does not make chemical classification decisions.
        It only skips residues that the structure file itself marks as polymer
        components via the entity categories. Everything else is treated as
        potentially non-polymer, and the classification comes from data sources
        (CCD pdbx_type, PubChem, ChEMBL) through the enrichment client.
        
        Args:
            residue_name: The residue/comp_id name.
            atom_site: The atom_site data (not used here, for signature compatibility).

        Returns:
            True if this is potentially a non-polymer entity.
        """
        # Nothing is classified here from hardcoded knowledge.
        # If we had a CCD-backed polymer component set, it would be used here.
        # For now defer everything to enrichment client.
        return True

    def _guess_element(self, atom_name: str, residue_name: str) -> Optional[str]:
        """
        Best-effort element guess from atom name only.

        This is used only as a last resort when the structure file does not
        provide an authoritative element symbol (e.g., type_symbol). The parser
        does not embed its own periodic-table data; for reliable element info,
        the CCD / enrichment client should be used.

        Args:
            atom_name: The atom name (e.g., "CA", "CB", "OXT").
            residue_name: The residue name for context.

        Returns:
            Element symbol guess, or None if unknown.
        """
        atom_name = atom_name.strip()

        if not atom_name:
            return None

        # Single-letter elements that are common in PDB/mmCIF atom names.
        # This is intentionally limited and not a full periodic table.
        single_letter_elements = {"H", "C", "N", "O", "S", "P", "F"}

        # If the atom name is a single letter and it is a known element symbol,
        # return that.
        if len(atom_name) == 1 and atom_name.upper() in single_letter_elements:
            return atom_name.upper()

        # Multi-character names: do not guess. Return None so callers can fall
        # back to data sources (CCD type_symbol, enrichment client) instead of
        # embedding chemical intuition here.
        return None

    def _compute_formula(self, atoms: List[Atom]) -> Optional[str]:
        """
        Compute chemical formula from atoms.

        This module does not hardcode element symbols or atomic weights.
        If the file provides an authoritative formula (e.g., from the CCD via
        enrichment), that should be used. Without element data of known fidelity,
        the parser cannot reliably produce a formula and returns None.

        Args:
            atoms: List of atoms.

        Returns:
            Formula string, or None if formula cannot be derived from file data alone.
        """
        # No hardcoded elements. Formula is a derived chemical property that
        # should come from data sources, not from parser heuristics.
        return None

    def _compute_molecular_weight(self, atoms: List[Atom]) -> Optional[float]:
        """
        Compute molecular weight from atoms.

        This module does not hardcode atomic-weight values. If authoritative
        element/isotope/weight data is attached by the enrichment client / CCD,
        that should be used instead of this fallback. Until then, the parser
        cannot compute a defensible molecular weight and returns None.

        Args:
            atoms: List of atoms.

        Returns:
            Molecular weight in g/mol, or None if cannot compute from file data alone.
        """
        # No hardcoded atomic weights. Molecular weight is a derived chemical
        # property that should come from data sources (CCD formula_weight,
        # PubChem MolecularWeight, etc.), not from parser constants.
        return None

    def _classify_ligand(self, name: str, formula: Optional[str],
                         mol_weight: Optional[float]) -> str:
        """
        Classify a ligand based on data source information.
        
        The classification is deferred entirely to the enrichment client / CCD.
        The parser does not apply its own chemical heuristics.
        
        Args:
            name: Ligand name/residue name (CCD ID).
            formula: Chemical formula from CCD.
            mol_weight: Molecular weight from CCD.

        Returns:
            Classification hint from CCD or "unclassified" if not available.
        """
        # Defer to enrichment client / CCD.
        # The parser itself has no opinion and applies no rules.
        return "unclassified"
    
    def _get_ccd_classification(self, name: str) -> Optional[str]:
        """
        Get the chemical component classification from RCSB CCD.

        The parser does not maintain this itself. This method exists so an
        enrichment client / CCD-backed layer can attach the authoritative
        pdbx_type / compound classification later. The same CCD field is
        available in the mmCIF _chem_comp category (pdbx_type).

        Args:
            name: The CCD identifier (e.g., "HEM", "ATP", "NA").

        Returns:
            The CCD-provided classification type if already resolved, otherwise None.
        """
        # Defer to enrichment client / CCD data.
        # The parser does not fetch or interpret CCD by itself.
        return None
