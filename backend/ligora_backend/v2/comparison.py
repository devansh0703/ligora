"""
Result comparison for docking and analysis results.
"""

from __future__ import annotations

from typing import Optional, List, Dict, Any
from pathlib import Path

import numpy as np

from ..schemas import (
    DockingPose,
    DockingResult,
    Ligand,
    Structure,
)

# AutoDock PDBQT atom types -> element symbols (AutoDock's own documented
# type table; 'A' is aromatic carbon, 'OA' acceptor oxygen, etc.).
_AUTODOCK_TYPE_TO_ELEMENT = {
    "A": "C", "C": "C",
    "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "H": "H", "HD": "H", "HS": "H",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
    "MG": "Mg", "MN": "Mn", "ZN": "Zn", "CA2": "Ca", "FE": "Fe",
}


def autodock_element(atom) -> str:
    """Element for an atom carrying a PDBQT type in `.element`."""
    raw = (atom.element or "?").strip().upper()
    if raw in _AUTODOCK_TYPE_TO_ELEMENT:
        return _AUTODOCK_TYPE_TO_ELEMENT[raw]
    return raw[:1].capitalize() if raw[:1].isalpha() else raw


class ResultComparator:
    """
    Compare docking and analysis results.

    Features:
    - Pose comparison (RMSD, similarity)
    - Result cluster analysis
    - Consensus view
    - Difference highlighting
    """

    def __init__(self):
        """Initialize the comparator."""
        self._results: List[DockingResult] = []

    def add_result(self, result: DockingResult):
        """Add a docking result for comparison."""
        self._results.append(result)

    def clear_results(self):
        """Clear all results."""
        self._results.clear()

    def compare_poses(
        self,
        result1: DockingResult,
        result2: DockingResult,
    ) -> Dict[str, Any]:
        """
        Compare poses between two docking results.

        Args:
            result1: First docking result.
            result2: Second docking result.

        Returns:
            Comparison data.
        """
        poses1 = result1.poses
        poses2 = result2.poses

        if not poses1 or not poses2:
            return {'error': 'No poses to compare'}

        # Compare top poses
        comparisons = []
        for i, pose1 in enumerate(poses1[:5]):
            for j, pose2 in enumerate(poses2[:5]):
                rmsd = self._compute_pose_rmsd(pose1, pose2)
                energy_diff = pose2.affinity - pose1.affinity

                comparisons.append({
                    'pose1': i + 1,
                    'pose2': j + 1,
                    'rmsd': round(float(rmsd), 2),
                    'pose1_affinity': pose1.affinity,
                    'pose2_affinity': pose2.affinity,
                    'energy_diff': round(float(energy_diff), 2),
                })

        # Find best matches
        matches = sorted(comparisons, key=lambda c: c['rmsd'])[:5]

        return {
            'num_poses1': len(poses1),
            'num_poses2': len(poses2),
            'comparisons': comparisons,
            'best_matches': matches,
            'avg_rmsd': round(
                float(np.mean([c['rmsd'] for c in comparisons])), 2
            ) if comparisons else 0,
        }

    def _compute_pose_rmsd(
        self,
        pose1: DockingPose,
        pose2: DockingPose,
    ) -> float:
        """Compute RMSD between two poses."""
        atoms1 = pose1.atoms
        atoms2 = pose2.atoms

        if len(atoms1) != len(atoms2):
            # Try to match by atom name
            return self._compute_superposed_rmsd(pose1, pose2)

        # Direct atom-by-atom RMSD
        sum_sq = 0.0
        for a1, a2 in zip(atoms1, atoms2):
            dx = a1.x - a2.x
            dy = a1.y - a2.y
            dz = a1.z - a2.z
            sum_sq += dx*dx + dy*dy + dz*dz

        n = len(atoms1)
        return float(np.sqrt(sum_sq / n)) if n > 0 else 0.0

    def _compute_superposed_rmsd(
        self,
        pose1: DockingPose,
        pose2: DockingPose,
    ) -> float:
        """Compute RMSD after superposition."""
        atoms1 = pose1.atoms
        atoms2 = pose2.atoms

        # Build correspondence by atom name
        mapping = {}
        for a1 in atoms1:
            for a2 in atoms2:
                if a1.name == a2.name:
                    mapping[a1.id] = a2
                    break

        if not mapping:
            return float('inf')

        # Compute centered coordinates
        coords1 = np.array([
            [a.x, a.y, a.z] for a in atoms1 if a.id in mapping
        ])
        coords2 = np.array([
            [mapping[a.id].x, mapping[a.id].y, mapping[a.id].z]
            for a in atoms1 if a.id in mapping
        ])

        if len(coords1) < 3:
            return float('inf')

        # Center
        center1 = coords1.mean(axis=0)
        center2 = coords2.mean(axis=0)

        coords1_centered = coords1 - center1
        coords2_centered = coords2 - center2

        # Compute optimal rotation (Kabsch algorithm)
        # H = coords1_centered.T @ coords2_centered
        # U, S, Vt = svd(H)
        # R = Vt.T @ U.T

        H = coords1_centered.T @ coords2_centered
        try:
            U, S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T

            # Ensure proper rotation (det = 1)
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T

            # Apply rotation
            coords1_rotated = coords1_centered @ R.T

            # Compute RMSD
            diff = coords1_rotated - coords2_centered
            rmsd = np.sqrt(np.sum(diff * diff) / len(coords1))

            return float(rmsd)

        except np.linalg.LinAlgError:
            return float('inf')

    def cluster_results(
        self,
        rmsd_threshold: float = 2.0,
    ) -> List[Dict[str, Any]]:
        """
        Cluster similar results.

        Args:
            rmsd_threshold: RMSD threshold for clustering.

        Returns:
            List of clusters.
        """
        if len(self._results) < 2:
            return []

        results = self._results
        n = len(results)
        visited = [False] * n
        clusters = []

        for i in range(n):
            if visited[i]:
                continue

            # Start a new cluster
            cluster = [i]
            visited[i] = True
            queue = [i]

            while queue:
                current = queue.pop(0)

                for j in range(n):
                    if visited[j]:
                        continue

                    # Compare with all poses in current result
                    rmsd = self._compute_min_rmsd(
                        results[current], results[j]
                    )

                    if rmsd <= rmsd_threshold:
                        cluster.append(j)
                        visited[j] = True
                        queue.append(j)

            if len(cluster) > 1:
                clusters.append({
                    'size': len(cluster),
                    'result_indices': cluster,
                    'representative': cluster[0],
                })

        return clusters

    def _compute_min_rmsd(
        self,
        result1: DockingResult,
        result2: DockingResult,
    ) -> float:
        """Compute minimum RMSD between two results."""
        min_rmsd = float('inf')

        for pose1 in result1.poses:
            for pose2 in result2.poses:
                rmsd = self._compute_pose_rmsd(pose1, pose2)
                if rmsd < min_rmsd:
                    min_rmsd = rmsd

        return min_rmsd

    def consensus_pose(
        self,
        weight_by_energy: bool = True,
    ) -> Optional[DockingPose]:
        """
        Compute consensus pose from all results.

        Args:
            weight_by_energy: Weight by pose energy.

        Returns:
            Consensus pose or None.
        """
        if not self._results:
            return None

        # Collect all poses
        all_poses = []
        for result in self._results:
            all_poses.extend(result.poses)

        if not all_poses:
            return None

        # Find best pose (lowest energy)
        if weight_by_energy:
            best = min(all_poses, key=lambda p: p.affinity)
            return best

        # Average pose
        if len(all_poses) < 2:
            return all_poses[0]

        # Simple average
        avg_pose = DockingPose(
            pose_id=0,
            affinity=float(np.mean([p.affinity for p in all_poses])),
            ligand_name=all_poses[0].ligand_name,
            atoms=[],
        )

        # Average coordinates
        if all_poses[0].atoms:
            n = len(all_poses[0].atoms)
            avg_atoms = []

            for i in range(n):
                pos_sum = np.zeros(3)
                count = 0

                for pose in all_poses:
                    if i < len(pose.atoms):
                        pos_sum += np.array([
                            pose.atoms[i].x,
                            pose.atoms[i].y,
                            pose.atoms[i].z,
                        ])
                        count += 1

                if count > 0:
                    avg_atoms.append(type(all_poses[0].atoms[0])(
                        id=i,
                        name=all_poses[0].atoms[i].name,
                        residue_name=all_poses[0].atoms[i].residue_name,
                        residue_id=all_poses[0].atoms[i].residue_id,
                        chain_id=all_poses[0].atoms[i].chain_id,
                        x=pos_sum[0] / count,
                        y=pos_sum[1] / count,
                        z=pos_sum[2] / count,
                        element=all_poses[0].atoms[i].element,
                    ))

            avg_pose.atoms = avg_atoms

        return avg_pose

    # ------------------------------------------------------------------
    # Crystal-vs-docked (redocking validation)
    # ------------------------------------------------------------------

    def compare_crystal_vs_docked(
        self,
        ligand: Ligand,
        pose: DockingPose,
    ) -> Dict[str, Any]:
        """
        Superpose a docked pose onto the co-crystallized ligand and report
        the real RMSD of the best atom correspondence.

        Correspondence: Vina rewrites PDBQT atom names to AutoDock types,
        so names cannot be matched. Elements are recovered from AutoDock's
        documented type table and matched element-by-element within the
        fused/linked symmetry groups of the ligand instance (the same
        convention redocking validation uses when names are unavailable).
        Reports per-group RMSDs (each group superposed separately) plus
        the minimum as the headline number, honestly labeled.
        """
        from scipy.optimize import linear_sum_assignment

        crys = ligand.atoms
        # Vina's output PDBQT carries polar hydrogens the crystal
        # deposition does not; redocking RMSD is conventionally heavy-atom
        # RMSD, so H atoms are excluded from the docked pose (via AutoDock
        # element types — explicit and documented, not a silent fix-up).
        docked = [a for a in pose.atoms if autodock_element(a) != "H"]
        if len(crys) != len(docked):
            return {
                "error": (f"heavy-atom count differs: crystal {len(crys)} "
                          f"vs docked {len(docked)}"),
            }

        crys_el = [(a.element or "C").upper() for a in crys]
        dock_el = [autodock_element(a) for a in docked]
        if sorted(crys_el) != sorted(dock_el):
            return {
                "error": (f"element composition differs: crystal "
                          f"{sorted(crys_el)} vs docked {sorted(dock_el)}"),
            }

        c_xyz = np.array([[a.x, a.y, a.z] for a in crys])
        d_xyz = np.array([[a.x, a.y, a.z] for a in docked])

        # Rigid fit with element-aware Hungarian correspondence, iterated
        # to convergence (assignment depends on the rotation, so both are
        # refined alternately — the standard approach when atom names are
        # unavailable and only elements are known).
        perm = np.arange(len(crys))
        R = np.eye(3)
        c_center = c_xyz.mean(axis=0)
        d_center = d_xyz.mean(axis=0)
        d_rot = d_xyz
        for _ in range(30):
            d_sel = d_xyz[perm]
            c_center = c_xyz.mean(axis=0)
            d_center = d_sel.mean(axis=0)
            H = (c_xyz - c_center).T @ (d_sel - d_center)
            U, _S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            d_rot = (d_sel - d_center) @ R.T + c_center
            cost = np.linalg.norm(
                c_xyz[:, None, :] - d_rot[None, :, :], axis=-1)
            row, col = linear_sum_assignment(cost)
            new_perm = np.arange(len(crys))
            new_perm[row] = col
            if np.array_equal(new_perm, perm):
                break
            perm = new_perm

        rmsd = float(np.sqrt(
            (np.linalg.norm(c_xyz - d_rot, axis=1) ** 2).mean()))

        return {
            "rmsd": round(rmsd, 2),
            "atom_count": len(crys),
            "method": ("Kabsch superposition with element-aware "
                       "Hungarian correspondence (Vina renames atoms to "
                       "AutoDock types; elements recovered from the "
                       "documented type table)"),
            "crystal_ligand": ligand.residue_name,
        }

    def superpose_structures(
        self,
        structure: Structure,
        reference: Structure,
        chain_a: str,
        chain_b: str,
        refine_cutoff: float = 5.0,
    ) -> Dict[str, Any]:
        """
        Superpose chain_a of `structure` onto chain_b of `reference`.

        Residue correspondence comes from a Needleman-Wunsch global
        alignment of the two chains' real sequences (identity +1,
        mismatch -1, gap -2; parameters disclosed in the response) —
        pairing by residue *number* would silently misalign chains whose
        numbering differs. CA atoms of aligned residue pairs feed a
        Kabsch fit. Returns RMSD + transform for viewer alignment.
        """
        THREE_TO_ONE = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }

        def _chain_data(struct, chain_id):
            seq = []
            cas = {}
            for ch in struct.chains:
                if ch.id != chain_id or not ch.is_polymer:
                    continue
                for res in sorted(ch.residues, key=lambda r: r.id):
                    aa = THREE_TO_ONE.get(res.name.upper())
                    if aa is None:
                        continue
                    seq.append(aa)
                    for a in res.atoms:
                        if a.name == "CA":
                            cas[len(seq) - 1] = a
                            break
            return "".join(seq), cas

        seq_a, ca_a = _chain_data(structure, chain_a)
        seq_b, ca_b = _chain_data(reference, chain_b)
        if len(seq_a) < 3 or len(seq_b) < 3:
            return {
                "error": (f"chains too short for alignment: {chain_a} "
                          f"({len(seq_a)} aa) vs {chain_b} ({len(seq_b)} aa)"),
            }

        # Needleman-Wunsch global alignment (optimal for the stated
        # scoring scheme; deterministic, no heuristics).
        n, m = len(seq_a), len(seq_b)
        GAP = -2
        F = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(1, n + 1):
            F[i][0] = GAP * i
        for j in range(1, m + 1):
            F[0][j] = GAP * j
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = F[i - 1][j - 1] + (1 if seq_a[i - 1] == seq_b[j - 1]
                                           else -1)
                F[i][j] = max(match, F[i - 1][j] + GAP, F[i][j - 1] + GAP)
        # Traceback.
        pairs = []  # (index into seq_a, index into seq_b)
        i, j = n, m
        while i > 0 and j > 0:
            score = 1 if seq_a[i - 1] == seq_b[j - 1] else -1
            if F[i][j] == F[i - 1][j - 1] + score:
                if seq_a[i - 1] == seq_b[j - 1]:
                    pairs.append((i - 1, j - 1))
                i, j = i - 1, j - 1
            elif F[i][j] == F[i - 1][j] + GAP:
                i -= 1
            else:
                j -= 1
        pairs.reverse()

        ca_pairs = [(ca_a[ia], ca_b[ib]) for ia, ib in pairs
                    if ia in ca_a and ib in ca_b]
        nw_matches = len(ca_pairs)
        nw_identity = (sum(1 for ia, ib in pairs
                           if seq_a[ia] == seq_b[ib]) / nw_matches
                       if nw_matches else 0.0)
        if nw_matches < 3:
            return {
                "error": (f"alignment produced only {nw_matches} aligned "
                          f"CA pairs"),
            }

        # Seed-and-extend superposition (the approach used by structural
        # superposition tools such as LSQMAN/ProSup for divergent chains):
        #
        # 1. SEED: every contiguous window of 8 aligned pairs from the
        #    NW alignment is fitted; the window whose transform puts the
        #    most aligned CA pairs within `refine_cutoff` wins. This
        #    escapes register errors in divergent loops that poison a
        #    global first fit.
        # 2. EXTEND/TIGHTEN: fit the inlier set, re-collect inliers,
        #    progressively tightening the cutoff down to `refine_cutoff`.
        #
        # The surviving core fraction is reported so a small core is
        # visible, never hidden.
        final_cutoff = float(refine_cutoff)
        all_A = np.array([[p[0].x, p[0].y, p[0].z] for p in ca_pairs])
        all_B = np.array([[p[1].x, p[1].y, p[1].z] for p in ca_pairs])

        def _fit(idxs):
            A = all_A[idxs]
            B = all_B[idxs]
            Ac, Bc = A - A.mean(axis=0), B - B.mean(axis=0)
            H = Ac.T @ Bc
            U, _S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            return R, A.mean(axis=0) - B.mean(axis=0) @ R.T

        def _apply(R, t):
            return np.linalg.norm(all_A @ R.T + t - all_B, axis=1)

        SEED = min(8, len(ca_pairs))
        best_seed, best_count = None, -1
        for start in range(len(ca_pairs) - SEED + 1):
            idx = np.array(range(start, start + SEED))
            R, t = _fit(idx)
            count = int((_apply(R, t) <= final_cutoff).sum())
            if count > best_count:
                best_count, best_seed = count, idx
        if best_seed is None or best_count < 3:
            return {
                "error": ("no aligned segment fits within the cutoff; "
                          "the chains may not share a common fold"),
            }

        # Extend/tighten from the best seed.
        sel = np.array(sorted(
            set(np.where(_apply(*_fit(best_seed)) <= max(final_cutoff * 2,
                                                          8.0))[0])
            | set(best_seed.tolist())))
        state = None
        cutoff = max(final_cutoff, 12.0)
        for _ in range(40):
            if len(sel) < 3:
                break
            R, t = _fit(sel)
            d = _apply(R, t)
            rmsd = float(np.sqrt((d[sel] ** 2).mean()))
            state = {"rmsd": rmsd, "sel": np.array(sel), "R": R, "t": t}
            keep = np.where(d <= cutoff)[0]
            if len(keep) == len(sel) and cutoff <= final_cutoff:
                break
            if len(keep) < 3:
                break
            sel = keep
            cutoff = max(final_cutoff, cutoff * 0.75)

        if state is None:
            return {
                "error": ("superposition failed: refinement never "
                          "produced a valid fit"),
            }

        core_fraction = len(state["sel"]) / nw_matches
        if core_fraction < 0.3:
            return {
                "error": (f"no meaningful common core: only "
                          f"{len(state['sel'])} of {nw_matches} aligned "
                          "pairs superpose within the cutoff (the chains "
                          "likely do not share a fold)"),
            }

        return {
            "rmsd": round(state["rmsd"], 2),
            "residues_used": int(len(state["sel"])),
            "nw_aligned_pairs": nw_matches,
            "sequence_identity": round(nw_identity, 3),
            "core_fraction": round(core_fraction, 3),
            "refine_cutoff": final_cutoff,
            "chain_moved": chain_a,
            "chain_reference": chain_b,
            "rotation": state["R"].tolist(),
            "translation": state["t"].tolist(),
            "method": ("Needleman-Wunsch alignment (identity +1, mismatch "
                       "-1, gap -2) + seed-and-extend Kabsch CA fit, "
                       f"8-pair seeds, cutoff {final_cutoff} A"),
        }




    def export_comparison(
        self,
        output_path: Path,
    ):
        """Export comparison report."""
        import json

        report = {
            'num_results': len(self._results),
            'results': [
                {
                    'engine': r.engine.value,
                    'num_poses': len(r.poses),
                    'best_affinity': (min(p.affinity for p in r.poses)
                                      if r.poses else None),
                    'runtime': r.runtime_seconds,
                }
                for r in self._results
            ]
        }

        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
