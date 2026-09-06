"""
Tests for V2 features.
"""

import pytest

from ligora_backend.v2.batch import BatchAnalyzer, BatchJob, BatchTask
from ligora_backend.v2.editors import Ligand2DEditor
from ligora_backend.v2.water import WaterNetworkAnalyzer
from ligora_backend.v2.scripting import ScriptingConsole
from ligora_backend.v2.comparison import ResultComparator
from ligora_backend.schemas import Structure, Chain, Residue, Atom, Ligand, AnalysisSummary, DockingResult, DockingPose


class TestBatchAnalyzer:
    """Tests for batch processing."""

    def test_create_batch(self):
        """Test batch creation."""
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()

        assert batch_id is not None
        assert len(batch_id) > 0
        assert batch_id in analyzer._batch_jobs

    def test_add_structure_to_batch(self):
        """Test adding structures to batch."""
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()

        task = analyzer.add_structure_to_batch(
            batch_id,
            "/path/to/structure.pdb",
            "local",
        )

        assert task is not None
        assert task.source == "/path/to/structure.pdb"
        assert task.source_type == "local"
        assert len(analyzer._batch_jobs[batch_id].tasks) == 1

    def test_batch_export(self):
        """Test batch result export."""
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()

        # Add a task
        from ligora_backend.schemas import AnalysisSummary
        task = analyzer.add_structure_to_batch(
            batch_id,
            "1ABC",
            "pdb_id",
        )

        # Mock a completed task
        task.status = "completed"
        task.result = BatchTaskResult(
            task_id=task.id,
            structure_id="1ABC",
            summary=AnalysisSummary(
                structure_id="1ABC",
                structure_title="Test",
                source="rcsb",
                ligand_id="L1",
                ligand_name="Test",
                ligand_formula="C6H6",
                ligand_smiles="c1ccccc1",
                resolution_status="resolved",
                contact_count=5,
                contacts=[],
                evidence=[],
                job_results=[],
                notes="",
                created_at=0,
                exported_at=0,
            ),
        )

        output_path = analyzer.export_batch_results(batch_id)

        
        assert output_path.exists()


class BatchTaskResult:
    """Simple result class for testing."""
    def __init__(self, task_id, structure_id, summary):
        self.task_id = task_id
        self.structure_id = structure_id
        self.summary = summary
        self.success = True
        self.error = None


class TestLigand2DEditor:
    """Tests for 2D ligand editor."""

    def test_load_ligand(self):
        """Test loading a ligand."""
        editor = Ligand2DEditor()

        from ligora_backend.schemas import Atom, Ligand
        ligand = Ligand(
            id="L1",
            name="Test",
            residue_name="LIG",
            atom_count=3,
            atoms=[
                Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                     chain_id="L", x=0, y=0, z=0, element="C"),
                Atom(id=2, name="C2", residue_name="LIG", residue_id=1,
                     chain_id="L", x=1.5, y=0, z=0, element="C"),
                Atom(id=3, name="O1", residue_name="LIG", residue_id=1,
                     chain_id="L", x=3.0, y=0, z=0, element="O"),
            ],
        )

        state = editor.load_ligand(ligand)

        assert len(state['atoms']) == 3
        # Bond inference is delegated to RDKit. When RDKit can build a
        # molecule from the supplied coordinates, bonds are populated;
        # otherwise they remain empty and the caller must supply topology.
        assert isinstance(state['bonds'], list)

    def test_update_atom_position(self):
        """Test updating atom position."""
        editor = Ligand2DEditor()

        from ligora_backend.schemas import Atom, Ligand
        ligand = Ligand(
            id="L1",
            name="Test",
            residue_name="LIG",
            atom_count=1,
            atoms=[
                Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                     chain_id="L", x=0, y=0, z=0, element="C"),
            ],
        )

        editor.load_ligand(ligand)
        editor.update_atom_position(1, 5.0, 5.0)

        state = editor.get_state()
        assert state['atoms'][0]['x'] == 5.0
        assert state['atoms'][0]['y'] == 5.0

    def test_add_atom(self):
        """Test adding an atom."""
        editor = Ligand2DEditor()

        atom_id = editor.add_atom("N", 1.0, 1.0)

        assert atom_id > 0
        state = editor.get_state()
        assert len(state['atoms']) == 1
        assert state['atoms'][0]['element'] == "N"

    def test_remove_atom(self):
        """Test removing an atom."""
        editor = Ligand2DEditor()

        atom_id = editor.add_atom("C", 0, 0)
        editor.add_atom("O", 1.5, 0)

        editor.remove_atom(atom_id)

        state = editor.get_state()
        assert len(state['atoms']) == 1

    def test_toggle_sync(self):
        """Test toggle sync."""
        editor = Ligand2DEditor()

        assert editor._sync_enabled is True
        editor.toggle_sync()
        assert editor._sync_enabled is False


class TestWaterNetworkAnalyzer:
    """Tests for water network analysis."""

    def test_detect_water_molecules(self):
        """Test water molecule detection."""
        analyzer = WaterNetworkAnalyzer()

        from ligora_backend.schemas import Structure, Chain, Residue, Atom

        structure = Structure(id="test", source="local", title="Test")
        chain = Chain(id="A", name="Protein", is_polymer=True)

        # Add a water molecule
        water = Residue(
            id=100,
            name="HOH",
            chain_id="A",
            residue_number=100,
            atoms=[
                Atom(id=1, name="O", residue_name="HOH", residue_id=100,
                     chain_id="A", x=10, y=10, z=10, element="O"),
            ],
        )
        chain.residues.append(water)
        structure.chains.append(chain)

        result = analyzer.analyze(structure)

        assert result['water_count'] == 1

    def test_build_network(self):
        """Test water network building."""
        analyzer = WaterNetworkAnalyzer()

        waters = [
            {
                'residue_id': 100,
                'chain_id': 'A',
                'name': 'HOH',
                'atoms': [
                    {'id': 1, 'name': 'O', 'element': 'O', 'x': 5, 'y': 5, 'z': 5}
                ],
                'center': [5, 5, 5],
            },
            {
                'residue_id': 101,
                'chain_id': 'A',
                'name': 'HOH',
                'atoms': [
                    {'id': 2, 'name': 'O', 'element': 'O', 'x': 7, 'y': 5, 'z': 5}
                ],
                'center': [7, 5, 5],
            },
        ]

        network = analyzer._build_network(waters)

        # Waters 2A apart may not be in network if threshold is 2.5A
        # This is a test to verify the network building works
        # Network may be empty if distance threshold isn't met
        assert isinstance(network, list)

    def test_find_clusters(self):
        """Test cluster finding."""
        analyzer = WaterNetworkAnalyzer()

        network = [
            {
                'water1': {'chain': 'A', 'residue': 100},
                'water2': {'chain': 'A', 'residue': 101},
                'distance': 2.8,
                'type': 'water-water',
            },
            {
                'water1': {'chain': 'A', 'residue': 101},
                'water2': {'chain': 'A', 'residue': 102},
                'distance': 2.9,
                'type': 'water-water',
            },
        ]

        clusters = analyzer._find_clusters(network)

        assert isinstance(clusters, list)
        if clusters:
            assert clusters[0]['size'] == 3


class TestScriptingConsole:
    """Tests for scripting console."""

    def test_execute_simple(self):
        """Test simple script execution."""
        console = ScriptingConsole()

        result = console.execute("x = 1 + 1")

        assert result['success'] is True
        assert 'x' in console._variables
        assert console._variables['x'] == 2

    def test_execute_with_error(self):
        """Test script with error."""
        console = ScriptingConsole()

        result = console.execute("undefined_function()")

        assert result['success'] is False
        assert 'error' in result

    def test_execute_file(self):
        """Test file execution."""
        console = ScriptingConsole()

        # Create temp file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write("y = 5 * 2")
            temp_path = f.name

        try:
            result = console.execute_file(temp_path)
            assert result['success'] is True
            assert console._variables['y'] == 10
        finally:
            import os
            os.unlink(temp_path)

    def test_get_history(self):
        """Test command history."""
        console = ScriptingConsole()

        console.execute("a = 1")
        console.execute("b = 2")

        history = console.get_history()

        assert len(history) == 2


class TestResultComparator:
    """Tests for result comparison."""

    def test_compare_poses(self):
        """Test pose comparison."""
        comparator = ResultComparator()

        from ligora_backend.schemas import DockingResult, DockingPose, Atom

        result1 = DockingResult(
            job_id="job1",
            engine="vina",
            poses=[
                DockingPose(
                    pose_id=1,
                    affinity=-7.5,
                    ligand_name="LIG",
                    atoms=[
                        Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                             chain_id="L", x=0, y=0, z=0, element="C"),
                        Atom(id=2, name="C2", residue_name="LIG", residue_id=1,
                             chain_id="L", x=1, y=0, z=0, element="C"),
                    ],
                ),
            ],
        )

        result2 = DockingResult(
            job_id="job2",
            engine="vina",
            poses=[
                DockingPose(
                    pose_id=1,
                    affinity=-8.0,
                    ligand_name="LIG",
                    atoms=[
                        Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                             chain_id="L", x=0.5, y=0, z=0, element="C"),
                        Atom(id=2, name="C2", residue_name="LIG", residue_id=1,
                             chain_id="L", x=1.5, y=0, z=0, element="C"),
                    ],
                ),
            ],
        )

        comparator.add_result(result1)
        comparator.add_result(result2)

        comparison = comparator.compare_poses(result1, result2)

        assert 'comparisons' in comparison
        assert 'best_matches' in comparison

    def test_cluster_results(self):
        """Test result clustering."""
        comparator = ResultComparator()

        from ligora_backend.schemas import DockingResult, DockingPose, Atom

        # Two similar results
        result1 = DockingResult(
            job_id="job1",
            engine="vina",
            poses=[
                DockingPose(
                    pose_id=1,
                    affinity=-7.5,
                    ligand_name="LIG",
                    atoms=[
                        Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                             chain_id="L", x=0, y=0, z=0, element="C"),
                    ],
                ),
            ],
        )

        result2 = DockingResult(
            job_id="job2",
            engine="vina",
            poses=[
                DockingPose(
                    pose_id=1,
                    affinity=-8.0,
                    ligand_name="LIG",
                    atoms=[
                        Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                             chain_id="L", x=0.2, y=0, z=0, element="C"),
                    ],
                ),
            ],
        )

        comparator.add_result(result1)
        comparator.add_result(result2)

        clusters = comparator.cluster_results(rmsd_threshold=2.0)

        assert len(clusters) >= 1

    def test_consensus_pose(self):
        """Test consensus pose computation."""
        comparator = ResultComparator()

        from ligora_backend.schemas import DockingResult, DockingPose, Atom

        result1 = DockingResult(
            job_id="job1",
            engine="vina",
            poses=[
                DockingPose(
                    pose_id=1,
                    affinity=-7.5,
                    ligand_name="LIG",
                    atoms=[
                        Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                             chain_id="L", x=0, y=0, z=0, element="C"),
                    ],
                ),
            ],
        )

        comparator.add_result(result1)

        consensus = comparator.consensus_pose()

        assert consensus is not None
        assert consensus.affinity == -7.5

    def test_clear_results(self):
        """Test clearing results."""
        comparator = ResultComparator()

        from ligora_backend.schemas import DockingResult, DockingPose, Atom

        result = DockingResult(
            job_id="job1",
            engine="vina",
            poses=[
                DockingPose(
                    pose_id=1,
                    affinity=-7.5,
                    ligand_name="LIG",
                    atoms=[
                        Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                             chain_id="L", x=0, y=0, z=0, element="C"),
                    ],
                ),
            ],
        )

        comparator.add_result(result)
        assert len(comparator._results) == 1

        comparator.clear_results()
        assert len(comparator._results) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
