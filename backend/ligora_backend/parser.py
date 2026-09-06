"""
Structure parser for mmCIF and PDB files.

Parses macromolecular structure files into an internal representation of
chains, residues, atoms, and non-polymer (ligand) entities.

Classification policy (no local heuristics):
- Polymer vs non-polymer status comes from the file's own _entity category
  (entity.type == "polymer").
- The `_atom_site.group_PDB` token (ATOM/HETATM) is part of the file itself
  and is used to route records to the polymer or non-polymer path.
- Every non-polymer component is surfaced as a candidate ligand. Whether a
  candidate is solvent, an ion, or a small molecule is decided ONLY by the
  RCSB Chemical Component Dictionary (CCD) `chem_comp.type` field, fetched
  live and attached by the ligand resolver / enrichment layer.
- Element symbols come from the file (`type_symbol` / columns 77-78). When a
  file does not provide them, the CCD's authoritative atom list is consulted
  (live); nothing is guessed from atom names.
"""
import re
import requests
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

from .schemas import (
    Structure, Chain, Residue, Atom, Ligand,
    LigandResolutionStatus
)
from .config import get_config


class StructureParser:
    """
    Parse macromolecular structure files (mmCIF, PDB).

    Supports:
    - Local file parsing (mmCIF primary, PDB legacy)
    - Remote fetching from RCSB (files.rcsb.org download service)
    """

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

        data = self._parse_mmcif_data(content)

        # Title: use the file's own description fields when present.
        title = self._get_nested_str(data, "struct.title")
        if not title:
            entities = self._get_loop_rows(data, "_entity")
            for ent in entities:
                if ent.get("type") == "polymer":
                    desc = ent.get("pdbx_description") or ent.get("details")
                    if desc:
                        title = desc
                        break
        structure.title = title or "Unknown structure"

        # Resolution: from the file's own refine / em_3d_reconstruction /
        # rcsb_entry_info categories when present.
        structure.resolution = self._read_resolution(data)
        structure.experiment_type = self._read_experiment_type(data)

        # Entity type map from the file's own _entity table.
        entity_type_map: Dict[str, str] = {}
        for entity in self._get_loop_rows(data, "_entity"):
            eid = entity.get("id")
            etype = entity.get("type")
            if eid is not None and etype:
                entity_type_map[str(eid)] = str(etype)

        # Entity -> strand mapping from _entity_poly (polymer entities).
        polymer_strands: Dict[str, List[str]] = {}
        entity_names: Dict[str, str] = {}
        for entity in self._get_loop_rows(data, "_entity"):
            eid = str(entity.get("id", ""))
            entity_names[eid] = entity.get("entity_name") or entity.get(
                "pdbx_description") or eid
        for poly in self._get_loop_rows(data, "_entity_poly"):
            eid = str(poly.get("entity_id", ""))
            strand_list = poly.get("pdbx_strand_id") or ""
            polymer_strands[eid] = [
                s.strip() for s in strand_list.split(",") if s.strip()
            ]

        chains: Dict[str, Chain] = {}

        def _get_chain(chain_id: str, is_polymer: bool) -> Chain:
            if chain_id not in chains:
                chains[chain_id] = Chain(
                    id=chain_id, name=chain_id, is_polymer=is_polymer)
            elif is_polymer and not chains[chain_id].is_polymer:
                chains[chain_id].is_polymer = True
            return chains[chain_id]

        # Register polymer chains from _entity_poly so they exist even if
        # some residues have no coordinates.
        for eid, strands in polymer_strands.items():
            if entity_type_map.get(eid) == "polymer":
                for sid in strands:
                    _get_chain(sid, True)

        # Parse atom_site rows.
        atom_site = data.get("_atom_site", [])
        if isinstance(atom_site, dict):
            atom_site = [atom_site]

        ligand_atoms: Dict[Tuple[str, str], List[Atom]] = {}
        serial_counter = 0

        for row in atom_site:
            group = (row.get("group_PDB") or
                     row.get("group_pdb") or "ATOM").strip().upper()
            atom_name = (row.get("label_atom_id") or
                         row.get("auth_atom_id") or "").strip()
            comp_id = (row.get("label_comp_id") or
                       row.get("auth_comp_id") or "").strip()
            label_asym = (row.get("label_asym_id") or "").strip()
            auth_asym = (row.get("auth_asym_id") or "").strip()
            chain_id = auth_asym or label_asym or " "
            label_seq = row.get("label_seq_id")
            auth_seq = row.get("auth_seq_id")
            seq_token = auth_seq if auth_seq not in (None, "", ".") else label_seq
            residue_number = self._parse_int(seq_token)
            label_entity = str(row.get("label_entity_id", "")).strip()

            try:
                x = float(row.get("Cartn_x"))
                y = float(row.get("Cartn_y"))
                z = float(row.get("Cartn_z"))
            except (TypeError, ValueError):
                # Model / deuterated rows without coordinates are skipped;
                # we only carry atoms the file actually located in space.
                continue

            b_factor = self._parse_float(row.get("B_iso_or_equiv"), 0.0)
            occupancy = self._parse_float(row.get("occupancy"), 1.0)
            element = (row.get("type_symbol") or "").strip().upper() or None

            serial_counter += 1
            atom = Atom(
                id=serial_counter,
                name=atom_name,
                residue_name=comp_id,
                residue_id=residue_number,
                chain_id=chain_id,
                x=x, y=y, z=z,
                element=element,
                b_factor=b_factor,
                occupancy=occupancy,
            )

            # Route by the file's own tokens:
            # - ATOM records go to the polymer path.
            # - HETATM records for entities the file declares as polymer
            #   (e.g. modified residues) also go to the polymer path.
            # - Everything else is a non-polymer entity (ligand candidate).
            is_polymer_record = (
                group == "ATOM" or entity_type_map.get(label_entity) == "polymer"
            )

            if is_polymer_record:
                chain = _get_chain(chain_id, True)
                residue = next(
                    (r for r in chain.residues
                     if r.id == residue_number and r.name == comp_id),
                    None
                )
                if residue is None:
                    residue = Residue(
                        id=residue_number,
                        name=comp_id,
                        chain_id=chain_id,
                        residue_number=residue_number,
                    )
                    chain.residues.append(residue)
                residue.atoms.append(atom)
            else:
                # Non-polymer entity: group by (component, instance chain).
                key = (comp_id, chain_id)
                ligand_atoms.setdefault(key, []).append(atom)

        # Build candidate ligands from the file's own non-polymer groups.
        # Identity/classification values (formula, weight, SMILES, type) are
        # attached later by the resolver from the live CCD.
        for (comp_id, chain_id), atoms in ligand_atoms.items():
            ligand = Ligand(
                id=f"L{comp_id}_{chain_id}" if chain_id else f"L{comp_id}",
                name=comp_id,
                residue_name=comp_id,
                formula=None,
                molecular_weight=None,
                atom_count=len(atoms),
                atoms=atoms,
                resolution_status=LigandResolutionStatus.NOT_FOUND,
                classification_hint=None,
            )
            structure.ligands.append(ligand)

        structure.chains = list(chains.values())

        structure.has_density = "_pdbx_vrpt" in data or any(
            k.startswith("_pdbx_volray") for k in data)

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
        ligand_atoms: Dict[Tuple[str, str], List[Atom]] = {}
        # Residues are accumulated per (chain, resseq, icode, comp) so
        # insertion codes and multi-chain files are handled correctly.
        residues: Dict[Tuple[str, int, str, str], Residue] = {}

        for line in lines:
            record = line[0:6].strip()

            if record in ("ATOM", "HETATM"):
                altloc = line[16:17].strip()
                if altloc and altloc != "A":
                    # Keep only the first alternate location per atom; the
                    # file's altLoc column is authoritative, not a heuristic.
                    continue
                atom_name = line[12:16].strip()
                # PDB name field is left-justified for element symbols
                # starting at column 14 when the name begins with a digit
                # (e.g. '1HG2'); for 1-char names the element sits in
                # column 14 (index 13). When the name field contains an
                # interior space the atom name is the first token.
                if " " in atom_name:
                    atom_name = atom_name.split()[0] if atom_name.split() \
                        else atom_name
                residue_name = line[17:20].strip()
                chain_id = line[21:22].strip() or " "
                residue_number = self._parse_int(line[22:26])
                insertion_code = line[26:27].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                occupancy = self._parse_float(line[54:60], 1.0)
                b_factor = self._parse_float(line[60:66], 0.0)
                element = line[76:78].strip().upper() or None

                atom = Atom(
                    id=self._parse_int(line[6:11]),
                    name=atom_name,
                    residue_name=residue_name,
                    residue_id=residue_number,
                    chain_id=chain_id,
                    x=x, y=y, z=z,
                    element=element,
                    b_factor=b_factor,
                    occupancy=occupancy,
                )

                if record == "HETATM":
                    key = (residue_name, chain_id)
                    ligand_atoms.setdefault(key, []).append(atom)
                else:
                    chain = chains.get(chain_id)
                    if chain is None:
                        chain = Chain(id=chain_id, name=chain_id,
                                      is_polymer=True)
                        chains[chain_id] = chain
                    res_key = (chain_id, residue_number, insertion_code,
                               residue_name)
                    residue = residues.get(res_key)
                    if residue is None:
                        residue = Residue(
                            id=residue_number,
                            name=residue_name,
                            chain_id=chain_id,
                            residue_number=residue_number,
                        )
                        residues[res_key] = residue
                        chain.residues.append(residue)
                    residue.atoms.append(atom)

            elif record == "TITLE":
                piece = line[10:80].rstrip()
                if piece:
                    # TITLE records continue across lines; PDB convention is
                    # that continuation lines repeat no leading whitespace.
                    if structure.title == "Unknown structure":
                        structure.title = piece.strip()
                    else:
                        existing = structure.title
                        addition = piece.strip()
                        if addition and not existing.endswith(addition):
                            structure.title = existing + " " + addition

            elif record == "REMARK" and line[7:10] == "   2":
                # REMARK   2 RESOLUTION. ### ANGSTROMS.
                m = re.search(r"RESOLUTION\.\s+([\d.]+)\s+ANGSTROM", line)
                if m:
                    try:
                        structure.resolution = float(m.group(1))
                    except ValueError:
                        pass

            elif record == "EXPDTA":
                technique = line[10:80].strip()
                if technique:
                    structure.experiment_type = technique

        structure.chains = list(chains.values())

        for (comp_id, chain_id), atoms in ligand_atoms.items():
            ligand = Ligand(
                id=f"L{comp_id}_{chain_id}" if chain_id else f"L{comp_id}",
                name=comp_id,
                residue_name=comp_id,
                formula=None,
                molecular_weight=None,
                atom_count=len(atoms),
                atoms=atoms,
                resolution_status=LigandResolutionStatus.NOT_FOUND,
                classification_hint=None,
            )
            structure.ligands.append(ligand)

        return structure

    def fetch_from_rcsb(self, pdb_id: str, format: str = "mmCIF") -> Structure:
        """
        Fetch a structure from RCSB (files.rcsb.org download service).

        Args:
            pdb_id: PDB ID (e.g., "1ABC") or CSM ID.
            format: File format ("mmCIF" or "PDB").

        Returns:
            Parsed Structure object.
        """
        config = get_config()
        pdb_id = pdb_id.strip().upper()

        if format.upper() in ("PDB",):
            download_url = f"{config.rcsb_files_url}/download/{pdb_id}.pdb"
        else:
            download_url = f"{config.rcsb_files_url}/download/{pdb_id}.cif"

        try:
            response = requests.get(
                download_url, timeout=config.request_timeout)
            response.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"Failed to fetch {pdb_id} from RCSB: {e}")

        if format.upper() == "PDB":
            return self.parse_pdb(response.text, source_id=pdb_id,
                                  source="rcsb")
        return self.parse_mmcif(response.text, source_id=pdb_id,
                                source="rcsb")

    def fetch_csm_from_rcsb(self, csm_id: str) -> Structure:
        """Fetch a Computed Structure Model from RCSB."""
        return self.fetch_from_rcsb(csm_id, format="mmCIF")

    def fetch_structure(self, identifier: str) -> Structure:
        """
        Fetch a structure from the appropriate source.

        Args:
            identifier: PDB ID, CSM ID, or local file path.

        Returns:
            Parsed Structure object.
        """
        local_path = Path(identifier)
        if local_path.exists():
            content = local_path.read_text(encoding="utf-8")
            if local_path.suffix.lower() in (".cif", ".mcif"):
                return self.parse_mmcif(content, source_id=local_path.stem)
            if local_path.suffix.lower() == ".pdb":
                return self.parse_pdb(content, source_id=local_path.stem)
            # Unknown extension: sniff the file's own content.
            if self._looks_like_mmcif(content):
                return self.parse_mmcif(content, source_id=local_path.stem)
            return self.parse_pdb(content, source_id=local_path.stem)

        if re.match(r'^[0-9][A-Za-z0-9]{3}$', identifier.strip()):
            return self.fetch_from_rcsb(identifier.strip())

        if identifier.startswith(("AF_", "MA_")):
            return self.fetch_csm_from_rcsb(identifier)

        raise ValueError(f"Unable to resolve identifier: {identifier}")

    # ------------------------------------------------------------------
    # CCD-backed element completion (data source, not heuristic)
    # ------------------------------------------------------------------

    def complete_missing_elements(self, structure: Structure) -> int:
        """
        Fill in missing element symbols from the RCSB CCD.

        The CCD's `atoms` category is the authoritative source for each
        component's atom types. Only atoms whose component could be resolved
        in the CCD are updated; everything else keeps element=None and the
        rest of the app reports it as unknown.

        Returns:
            Number of atoms updated.
        """
        from .ligand import ccd_element_map

        updated = 0
        components: Dict[str, Dict[str, Optional[str]]] = {}

        all_atoms: List[Atom] = []
        for chain in structure.chains:
            for residue in chain.residues:
                all_atoms.extend(residue.atoms)
        for ligand in structure.ligands:
            all_atoms.extend(ligand.atoms)

        for atom in all_atoms:
            if atom.element:
                continue
            comp = atom.residue_name
            if comp not in components:
                components[comp] = ccd_element_map(comp)
            mapping = components[comp]
            if mapping and atom.name in mapping and mapping[atom.name]:
                atom.element = mapping[atom.name]
                updated += 1

        return updated

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _looks_like_mmcif(self, content: str) -> bool:
        """Detect mmCIF by the file's own data_ header token."""
        for line in content.splitlines()[:20]:
            s = line.strip()
            if s.startswith("data_"):
                return True
            if s.startswith(("HEADER", "ATOM  ", "HETATM", "CRYST1")):
                return False
        return False

    def _read_resolution(self, data: Dict[str, Any]) -> Optional[float]:
        """Read resolution from the file's own quality categories."""
        # refine.ls_d_res_high (X-ray) - loop or scalar category
        for row in self._get_loop_rows(data, "_refine"):
            val = row.get("ls_d_res_high")
            if val:
                try:
                    return float(val)
                except (ValueError, TypeError):
                    pass
        # em_3d_reconstruction.resolution (cryo-EM)
        for row in self._get_loop_rows(data, "_em_3d_reconstruction"):
            val = row.get("resolution")
            if val:
                try:
                    return float(val)
                except (ValueError, TypeError):
                    pass
        # rcsb_entry_info.resolution_combined (loop or scalar)
        for row in self._get_loop_rows(data, "_rcsb_entry_info"):
            v = row.get("resolution_combined")
            if v:
                first = str(v).split(",")[0].strip("[]' ")
                try:
                    return float(first)
                except ValueError:
                    pass
        return None

    def _read_experiment_type(self, data: Dict[str, Any]) -> Optional[str]:
        """Read the experimental method from the file's own _exptl table."""
        rows = self._get_loop_rows(data, "_exptl")
        for row in rows:
            method = row.get("method")
            if method:
                return str(method)
        return None

    @staticmethod
    def _parse_int(token) -> int:
        """Parse an integer token; unparseable tokens yield 0."""
        try:
            return int(str(token).strip())
        except (ValueError, TypeError):
            # Strip insertion codes like "100A"
            m = re.match(r"^(-?\d+)", str(token).strip())
            return int(m.group(1)) if m else 0

    @staticmethod
    def _parse_float(token, default: float) -> float:
        """Parse a float token; unparseable tokens yield the default."""
        try:
            return float(str(token).strip())
        except (ValueError, TypeError):
            return default

    def _parse_mmcif_data(self, content: str) -> Dict[str, Any]:
        """
        Parse mmCIF data into a nested dictionary structure.

        Single-value categories (``_category.field value``) are stored as
        nested dicts; ``loop_`` tables are stored as lists of row dicts.
        Values use the CIF conventions: ``?`` = unknown, ``.`` = null,
        quoting with single/double quotes, and ``;``-delimited multi-line
        text blocks.
        """
        data: Dict[str, Any] = {}
        lines = content.splitlines()
        n = len(lines)
        i = 0

        while i < n:
            line = lines[i].strip()
            i += 1

            if not line or line.startswith('#') or line.startswith('data_'):
                continue

            if line == 'loop_':
                fields: List[str] = []
                rows: List[Dict[str, str]] = []
                category: Optional[str] = None

                # Read column headers.
                while i < n:
                    h = lines[i].strip()
                    if not h or h.startswith('#'):
                        i += 1
                        continue
                    if h == 'loop_' or not h.startswith('_'):
                        break
                    col = h.split()[0]
                    if '.' in col:
                        cat, field = col.split('.', 1)
                        category = cat.lstrip('_')
                    else:
                        field = col
                        category = category or col.lstrip('_')
                    fields.append(field)
                    i += 1

                if not fields:
                    continue

                # Read data rows using the tokenizer.
                tokenizer = _CifLoopTokenizer(fields, rows)
                while i < n:
                    raw = lines[i]
                    s = raw.strip()
                    if s == 'loop_':
                        break
                    if s.startswith('_') and not tokenizer.open_quote:
                        break
                    if s.startswith('#') and not tokenizer.open_quote:
                        i += 1
                        continue
                    i += 1
                    tokenizer.feed(raw)
                data['_' + category] = rows
                continue

            if line.startswith('_'):
                parts = line.split(None, 1)
                field_full = parts[0]
                value_raw = parts[1].strip() if len(parts) > 1 else ''

                cat_part, _, field_part = field_full.partition('.')
                cat_name = cat_part.lstrip('_')
                if not field_part:
                    # "_field value" without a category.
                    cat_name, field_part = field_full.lstrip('_'), 'value'

                if value_raw.startswith(';'):
                    # Multi-line text block: value continues until a line
                    # that is exactly ';'.
                    buf = [value_raw[1:]]
                    closed = False
                    while i < n:
                        t = lines[i]
                        i += 1
                        if t.rstrip() == ';':
                            closed = True
                            break
                        buf.append(t)
                    value = '\n'.join(buf) if closed else '\n'.join(buf)
                    data.setdefault('_' + cat_name, {})[field_part] = \
                        self._parse_scalar(value.strip())
                    continue

                if not value_raw:
                    # Value on following lines (quoted continuation).
                    if i < n and lines[i].strip().startswith(';'):
                        i += 1
                        buf = []
                        while i < n:
                            t = lines[i]
                            i += 1
                            if t.rstrip() == ';':
                                break
                            buf.append(t)
                        data.setdefault('_' + cat_name, {})[field_part] = \
                            self._parse_scalar('\n'.join(buf).strip())
                        continue
                    data.setdefault('_' + cat_name, {})[field_part] = None
                    continue

                data.setdefault('_' + cat_name, {})[field_part] = \
                    self._parse_scalar(value_raw)
                continue

        return data

    def _tokenize_row(
        self,
        line: str,
        fields: List[str],
        rows: List[Dict[str, str]],
        pending: List[str],
    ) -> List[str]:
        """Legacy single-line tokenizer entry point (kept for tests).

        Feeds one physical line into a fresh tokenizer and returns the
        still-pending token list.
        """
        tok = _CifLoopTokenizer(fields, rows, pending=pending)
        tok.feed(line)
        return tok.pending_tokens

    def _parse_scalar(self, raw: str) -> Any:
        """Convert a single mmCIF scalar value to a Python value."""
        s = raw.strip()
        if s in ('', '.', '?'):
            return None
        if (s.startswith("'") and s.endswith("'")) or \
           (s.startswith('"') and s.endswith('"')):
            return s[1:-1]
        lowered = s.lower()
        if lowered in ('true', 'yes'):
            return True
        if lowered in ('false', 'no'):
            return False
        try:
            return int(s)
        except ValueError:
            pass
        try:
            return float(s)
        except ValueError:
            pass
        return s

    def _get_loop_rows(self, data: Dict[str, Any],
                       category: str) -> List[Dict[str, Any]]:
        """Return the rows of a mmCIF category (loop or single-value)."""
        rows = data.get(category)
        if isinstance(rows, list):
            return rows
        if isinstance(rows, dict) and rows:
            return [rows]
        return []

    def _get_nested_str(self, data: Dict[str, Any],
                        dotted_key: str) -> Optional[str]:
        """Resolve a dotted key like ``refine.ls_d_res_high``."""
        cur: Any = data
        for part in dotted_key.split('.'):
            if not isinstance(cur, dict):
                return None
            cur = cur.get(part)
            if cur is None:
                return None
        if isinstance(cur, str):
            return cur
        if isinstance(cur, (int, float)):
            return str(cur)
        return None


class _CifLoopTokenizer:
    """Stateful tokenizer for mmCIF loop data sections.

    Implements CIF 1.1 value quoting: a quote character delimits a value
    only when followed by whitespace (or end of line), so embedded
    apostrophes/quotes inside values are preserved. Values may span
    multiple physical lines. Logical rows are emitted once one token per
    column has accumulated.
    """

    def __init__(self, fields: List[str], rows: List[Dict[str, str]],
                 pending: Optional[List[str]] = None):
        self.fields = fields
        self.rows = rows
        self.tokens: List[str] = list(pending or [])
        self.open_quote: Optional[str] = None
        self._buf: List[str] = []

    def feed(self, line: str):
        """Feed one physical line into the tokenizer."""
        i = 0
        n = len(line)
        while i < n:
            ch = line[i]

            if self.open_quote is not None:
                # Inside a quoted value: a matching quote followed by
                # whitespace/EOL closes the value.
                if ch == self.open_quote and \
                        (i + 1 >= n or line[i + 1].isspace()):
                    self.tokens.append(''.join(self._buf))
                    self._buf = []
                    self.open_quote = None
                    i += 1
                    continue
                self._buf.append(ch)
                i += 1
                continue

            if ch.isspace():
                i += 1
                continue

            if ch in "'\"" and (i + 1 >= n or not line[i + 1].isspace()):
                # Opening quote of a quoted value.
                self.open_quote = ch
                i += 1
                continue

            # Bare token: read until whitespace or an opening quote.
            j = i
            while j < n and not line[j].isspace() and line[j] not in "'\"":
                j += 1
            self.tokens.append(line[i:j])
            i = j

        self._try_emit_row()

    def _try_emit_row(self):
        if self.open_quote is None and \
                len(self.tokens) >= len(self.fields):
            row: Dict[str, str] = {}
            for idx, field in enumerate(self.fields):
                row[field] = self.tokens[idx]
            self.rows.append(row)
            self.tokens = self.tokens[len(self.fields):]

    @property
    def pending_tokens(self) -> List[str]:
        return self.tokens
