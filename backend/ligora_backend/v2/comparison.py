"""
Result comparison for docking and analysis results.
"""

from typing import Optional, List, Dict, Any
from pathlib import Path

import numpy as np

from ..schemas import (
    DockingPose,
    DockingResult,
)


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
                    'best_affinity': min(p.affinity for p in r.poses) if r.poses else None,
                    'runtime': r.runtime_seconds,
                }
                for r in self._results
            ]
        }

        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
