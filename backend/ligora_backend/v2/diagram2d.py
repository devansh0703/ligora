"""
2D interaction diagram from real analysis data.

Layout inputs, all real:
- Ligand atom positions: RDKit 2D depiction of the resolved CCD chemistry
  when RDKit can build the molecule; otherwise a principal-component
  projection of the structure's actual 3D coordinates (pure coordinate
  math, no invented chemistry).
- Ligand bonds: the CCD's own _chem_comp_bond records (via the same
  Cheminformatics path the 2D editor uses).
- Residues shown: only residues with real PLIP contacts.
- Edge geometry: the real contact types and distances from PLIP.

The diagram is an SVG; residue nodes are placed on a circle around the
ligand (angular positions from each residue's real direction relative to
the ligand center — projection of the actual 3D geometry, not a layout
heuristic about chemistry).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from ..cheminformatics import Cheminformatics
from ..schemas import Contact, Ligand, Structure


def _escape(text: str) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def _pca_projection(points: List[Tuple[float, float, float]]
                    ) -> List[Tuple[float, float]]:
    """2D projection of 3D coordinates onto the two largest-variance axes."""
    import numpy as np
    arr = np.array(points, dtype=float)
    centered = arr - arr.mean(axis=0)
    # Singular vectors of the centered coordinates span the principal axes.
    _u, _s, vt = np.linalg.svd(centered, full_matrices=False)
    proj = centered @ vt[:2].T
    return [(float(x), float(y)) for x, y in proj]


def _rdkit_2d_positions(
    ligand: Ligand,
) -> Optional[Tuple[List[Tuple[float, float]], List[Dict[str, Any]]]]:
    """
    Real 2D depiction coordinates + CCD bond list via RDKit.

    Returns (positions_by_atom_id, bonds) or None when RDKit cannot build
    the molecule (e.g. unresolved chemistry). Bonds come from the CCD
    records inside Cheminformatics (no distance-based guessing).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None
    chem = Cheminformatics()
    atoms = ligand.atoms
    mol_block = chem.mol_block_from_atoms(atoms)
    if mol_block is None:
        return None
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False, sanitize=True)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()

    # Mol-block atom order: input atoms first (same order), then any
    # appended hydrogens — same mapping contract as the 2D editor.
    n_input = len(atoms)
    positions: List[Tuple[float, float]] = []
    for i, atom in enumerate(atoms):
        if i >= mol.GetNumAtoms():
            break
        p = conf.GetAtomPosition(i)
        positions.append((float(p.x), float(p.y)))
    if len(positions) < n_input:
        return None

    bonds: List[Dict[str, Any]] = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i >= n_input or j >= n_input:
            continue
        ai, aj = atoms[i], atoms[j]
        if ((ai.element or "").upper() == "H" or
                (aj.element or "").upper() == "H"):
            continue
        order = {1.0: 1, 1.5: 2, 2.0: 2, 3.0: 3}.get(
            bond.GetBondTypeAsDouble(), 1)
        bonds.append({"from": ai.id, "to": aj.id, "order": order})
    return positions, bonds


_CONTACT_COLORS = {
    "hydrogen_bond": "#2f80ed",
    "hydrophobic": "#f2994a",
    "pi_stacking": "#9b51e0",
    "salt_bridge": "#eb5757",
    "halogen_bond": "#56ccf2",
    "metal_coordination": "#b3782d",
    "water_mediated": "#6fcf97",
    "van_der_waals": "#bdbdbd",
    "unknown": "#828282",
}

_CONTACT_DASH = {
    "hydrogen_bond": "6 3",
    "salt_bridge": "6 3",
    "halogen_bond": "6 3",
    "water_mediated": "2 3",
}


def build_interaction_diagram(
    structure: Structure,
    ligand: Ligand,
    contacts: List[Contact],
    width: int = 900,
    height: int = 640,
) -> Dict[str, Any]:
    """
    Build the SVG interaction diagram from the structure, the ligand's
    real chemistry, and PLIP's real contacts.

    Returns {'svg': str, 'residues': [...], 'ligand_bounds': {...}} — the
    residue list doubles as the structured data behind the diagram.
    """
    if not ligand.atoms:
        raise ValueError("Selected ligand has no atoms")

    # --- Ligand 2D layout -------------------------------------------------
    rdkit_layout = _rdkit_2d_positions(ligand)
    if rdkit_layout is not None:
        raw_positions, bonds = rdkit_layout
        layout_source = "rdkit_2d_depiction"
    else:
        raw_positions = _pca_projection(
            [(a.x, a.y, a.z) for a in ligand.atoms])
        bonds = []
        layout_source = "pca_projection_of_3d_coordinates"

    positions = _normalize_layout(raw_positions, width, height,
                                  reserve=0.42)

    # --- Residues from real contacts --------------------------------------
    residue_ids: List[str] = []
    for c in contacts:
        key = f"{c.protein_chain_id}:{c.protein_residue_id}"
        if key not in residue_ids:
            residue_ids.append(key)

    lig_center = _ligand_center_2d(positions)
    residues = []
    for key in residue_ids:
        chain, _, resid = key.partition(":")
        first = next(c for c in contacts
                     if f"{c.protein_chain_id}:{c.protein_residue_id}" == key)
        name = first.protein_residue_name
        # Angular position from the residue's real 3D direction relative
        # to the ligand, projected onto the ligand's principal plane.
        theta = _residue_angle(ligand, structure, chain, resid)
        rx, ry = _ring_position(lig_center, theta, width, height)
        residues.append({
            "key": key,
            "chain": chain,
            "residue_id": int(resid),
            "name": name,
            "x": round(rx, 1),
            "y": round(ry, 1),
        })
    residue_pos = {r["key"]: (r["x"], r["y"]) for r in residues}

    # --- Contact edges -----------------------------------------------------
    edges = []
    for c in contacts:
        key = f"{c.protein_chain_id}:{c.protein_residue_id}"
        rx, ry = residue_pos.get(key, (0, 0))
        latt = next((a for a in ligand.atoms
                     if a.name == c.ligand_atom), None)
        lx, ly = positions.get(latt.id, lig_center) if latt else lig_center
        edges.append({
            "x1": round(lx, 1), "y1": round(ly, 1),
            "x2": rx, "y2": ry,
            "color": _CONTACT_COLORS.get(c.contact_type.value,
                                         _CONTACT_COLORS["unknown"]),
            "dash": _CONTACT_DASH.get(c.contact_type.value),
            "type": c.contact_type.value,
            "distance": c.distance,
            "ligand_atom": c.ligand_atom,
            "protein_atom": c.protein_atom,
        })

    # --- SVG ---------------------------------------------------------------
    parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}" '
        f'font-family="ui-sans-serif, system-ui, sans-serif">',
        f'<rect width="{width}" height="{height}" fill="transparent"/>',
    ]
    # Bonds first (under atoms)
    for b in bonds:
        (x1, y1) = positions.get(b["from"], (0, 0))
        (x2, y2) = positions.get(b["to"], (0, 0))
        for offset in range(b["order"]):
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / length, dx / length
            shift = (offset - (b["order"] - 1) / 2) * 2.6
            parts.append(
                f'<line x1="{x1 + nx*shift:.1f}" y1="{y1 + ny*shift:.1f}" '
                f'x2="{x2 + nx*shift:.1f}" y2="{y2 + ny*shift:.1f}" '
                f'stroke="#7a828a" stroke-width="1.6"/>')
    # Contact edges
    for e in edges:
        dash = f' stroke-dasharray="{e["dash"]}"' if e["dash"] else ""
        parts.append(
            f'<line x1="{e["x1"]}" y1="{e["y1"]}" x2="{e["x2"]}" '
            f'y2="{e["y2"]}" stroke="{e["color"]}" stroke-width="1.8" '
            f'stroke-opacity="0.85"{dash}/>')
        mx, my = (e["x1"] + e["x2"]) / 2, (e["y1"] + e["y2"]) / 2
        parts.append(
            f'<text x="{mx:.1f}" y="{my - 3:.1f}" font-size="9" '
            f'fill="{e["color"]}" text-anchor="middle">'
            f'{_escape(e["type"].replace("_", " "))} '
            f'{e["distance"]:.2f}</text>')
    # Ligand atoms
    atom_elements = {"N": "#3050c8", "O": "#e04030", "S": "#c8a020",
                     "P": "#ff8000", "FE": "#b07a30", "MG": "#40a050",
                     "ZN": "#7a60a0", "CA": "#808090"}
    for atom in ligand.atoms:
        x, y = positions.get(atom.id, (0, 0))
        elem = (atom.element or "C").upper()
        color = atom_elements.get(elem, "#3a3f45")
        r = 5.2 if elem != "H" else 3.0
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{color}" '
            f'stroke="#1c2024" stroke-width="0.8"/>')
        if elem not in ("C", "H"):
            parts.append(
                f'<text x="{x:.1f}" y="{y + 3:.1f}" font-size="8" '
                f'fill="#fff" text-anchor="middle">{_escape(elem)}</text>')
    # Residue nodes (on top)
    for r in residues:
        parts.append(
            f'<g class="diagram-residue" data-chain="{_escape(r["chain"])}" '
            f'data-resi="{r["residue_id"]}">'
            f'<circle cx="{r["x"]}" cy="{r["y"]}" r="17" fill="#22262b" '
            f'stroke="#566070" stroke-width="1.4"/>'
            f'<text x="{r["x"]}" y="{r["y"] + 4}" font-size="11" '
            f'fill="#e8eaed" text-anchor="middle">'
            f'{_escape(r["name"][:3])}{r["residue_id"]}</text></g>')
    # Legend from the contact types actually present
    present = sorted({e["type"] for e in edges})
    lx, ly = 12, height - 12 - 14 * len(present)
    for i, t in enumerate(present):
        color = _CONTACT_COLORS.get(t, _CONTACT_COLORS["unknown"])
        parts.append(
            f'<line x1="{lx}" y1="{ly + i*14}" x2="{lx + 22}" '
            f'y2="{ly + i*14}" stroke="{color}" stroke-width="2"/>'
            f'<text x="{lx + 28}" y="{ly + i*14 + 3}" font-size="10" '
            f'fill="#aab2bc">{_escape(t.replace("_", " "))}</text>')
    parts.append("</svg>")

    return {
        "svg": "".join(parts),
        "layout_source": layout_source,
        "residues": residues,
        "contact_count": len(contacts),
        "residue_count": len(residues),
        "ligand_name": ligand.residue_name,
    }


def _normalize_layout(
    raw: List[Tuple[float, float]],
    width: int,
    height: int,
    reserve: float,
) -> Dict[int, Tuple[float, float]]:
    """Scale raw layout units into the central diagram area."""
    xs = [p[0] for p in raw]
    ys = [p[1] for p in raw]
    min_x, max_x, min_y, max_y = min(xs), max(xs), min(ys), max(ys)
    span_x = (max_x - min_x) or 1.0
    span_y = (max_y - min_y) or 1.0
    cx = width / 2
    cy = height / 2
    avail_w = width * (1 - 2 * reserve)
    avail_h = height * (1 - 2 * reserve)
    scale = min(avail_w / span_x, avail_h / span_y)
    return {
        i + 1: (cx + (x - (min_x + max_x) / 2) * scale,
                cy + (y - (min_y + max_y) / 2) * scale)
        for i, (x, y) in enumerate(raw)
    }


def _ligand_center_2d(positions: Dict[int, Tuple[float, float]]
                      ) -> Tuple[float, float]:
    xs = [p[0] for p in positions.values()]
    ys = [p[1] for p in positions.values()]
    return (sum(xs) / len(xs), sum(ys) / len(ys))


def _residue_angle(ligand: Ligand, structure: Structure, chain_id: str,
                   resid: str) -> float:
    """
    Screen angle of a residue around the ligand, from real 3D geometry.

    The ligand's principal axes (SVD of its centered coordinates) define
    a plane; the residue center is projected onto that plane and its
    angle measured. Pure coordinate math — reproducible and honest.
    """
    import numpy as np
    lig_pts = np.array([[a.x, a.y, a.z] for a in ligand.atoms],
                       dtype=float)
    lig_center = lig_pts.mean(axis=0)
    centered = lig_pts - lig_center
    _u, _s, vt = np.linalg.svd(centered, full_matrices=False)
    axis_u, axis_v = vt[0], vt[1]  # two largest-variance directions

    res_pts = []
    for chain in structure.chains:
        if chain.id != chain_id:
            continue
        for residue in chain.residues:
            if str(residue.id) == resid:
                res_pts = [(a.x, a.y, a.z) for a in residue.atoms]
                break
    if not res_pts:
        return 0.0
    res_center = np.array(res_pts).mean(axis=0)
    delta = res_center - lig_center
    u = float(delta @ axis_u)
    v = float(delta @ axis_v)
    return math.atan2(v, u)


def _ring_position(center, theta: float, width: int,
                   height: int) -> Tuple[float, float]:
    """Place a residue node on an ellipse around the ligand at angle theta."""
    rx = width * 0.36
    ry = height * 0.36
    return (center[0] + rx * math.cos(theta),
            center[1] + ry * math.sin(theta))
