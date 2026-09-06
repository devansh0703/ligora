"""
Tests for V2 features: batch, editors, water, scripting, comparison.

Pure-unit where possible; network-dependent behavior is exercised in the
server-level live tests, not here.
"""

import pytest

from ligora_backend.v2.batch import BatchAnalyzer
from ligora_backend.v2.editors import Ligand2DEditor
from ligora_backend.v2.water import WaterNetworkAnalyzer
from ligora_backend.v2.scripting import ScriptingConsole
from ligora_backend.v2.comparison import ResultComparator
from ligora_backend.schemas import (
    Structure, Chain, Residue, Atom, Ligand,
    DockingResult, DockingPose,
)


def make_editor_ligand() -> Ligand:
    return Ligand(
        id="L1", name="Test", residue_name="LIG", atom_count=3,
        atoms=[
            Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=0, y=0, z=0, element="C"),
            Atom(id=2, name="C2", residue_name="LIG", residue_id=1,
                 chain_id="L", x=1.5, y=0, z=0, element="C"),
            Atom(id=3, name="O1", residue_name="LIG", residue_id=1,
                 chain_id="L", x=3.0, y=0, z=0, element="O"),
        ],
    )


class TestBatchAnalyzer:
    """Tests for batch processing (offline paths only)."""

    def test_create_batch(self):
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()
        assert batch_id in analyzer._batch_jobs

    def test_add_structure_to_batch(self):
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()
        task = analyzer.add_structure_to_batch(
            batch_id, "1ABC", "pdb_id")
        assert task is not None
        assert task.source == "1ABC"
        assert task.source_type == "pdb_id"
        assert len(analyzer._batch_jobs[batch_id].tasks) == 1

    def test_add_to_unknown_batch_returns_none(self):
        analyzer = BatchAnalyzer()
        assert analyzer.add_structure_to_batch(
            "nope", "1ABC", "pdb_id") is None

    def test_run_empty_batch(self):
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()
        results = analyzer.run_batch(batch_id)
        assert results == []

    def test_run_batch_reports_load_failure_honestly(self):
        analyzer = BatchAnalyzer()
        batch_id = analyzer.create_batch()
        task = analyzer.add_structure_to_batch(
            batch_id, "/nonexistent/structure.cif", "local")
        results = analyzer.run_batch(batch_id)
        # The failure is reported honestly: as a failed BatchResult AND on
        # the task itself, so it appears in summaries and exports.
        assert len(results) == 1
        assert results[0].success is False
        assert results[0].error == "Failed to load structure"
        assert task.status == "failed"
        assert task.error == "Failed to load structure"


class TestLigand2DEditor:
    """Tests for the 2D ligand editor."""

    def test_load_ligand(self):
        editor = Ligand2DEditor()
        state = editor.load_ligand(make_editor_ligand())
        assert len(state['atoms']) == 3
        assert isinstance(state['bonds'], list)

    def test_update_atom_position(self):
        editor = Ligand2DEditor()
        ligand = Ligand(
            id="L1", name="Test", residue_name="LIG", atom_count=1,
            atoms=[Atom(id=1, name="C1", residue_name="LIG", residue_id=1,
                        chain_id="L", x=0, y=0, z=0, element="C")],
        )
        editor.load_ligand(ligand)
        editor.update_atom_position(1, 5.0, 5.0)
        state = editor.get_state()
        assert state['atoms'][0]['x'] == 5.0
        assert state['atoms'][0]['y'] == 5.0

    def test_add_atom(self):
        editor = Ligand2DEditor()
        atom_id = editor.add_atom("N", 1.0, 1.0)
        assert atom_id > 0
        state = editor.get_state()
        assert state['atoms'][0]['element'] == "N"

    def test_remove_atom(self):
        editor = Ligand2DEditor()
        first = editor.add_atom("C", 0, 0)
        editor.add_atom("O", 1.5, 0)
        editor.remove_atom(first)
        assert len(editor.get_state()['atoms']) == 1

    def test_add_and_remove_bond(self):
        editor = Ligand2DEditor()
        a = editor.add_atom("C", 0, 0)
        b = editor.add_atom("O", 1.5, 0)
        editor.add_bond(a, b, 2)
        assert {'from': a, 'to': b, 'order': 2} in \
            editor.get_state()['bonds']
        editor.remove_bond(a, b)
        assert editor.get_state()['bonds'] == []

    def test_bond_order_lookup_no_heuristics(self):
        """_estimate_bond_order reads the inferred topology only."""
        editor = Ligand2DEditor()
        a = editor.add_atom("C", 0, 0)
        b = editor.add_atom("O", 1.5, 0)
        editor.add_bond(a, b, 2)
        c = editor.add_atom("N", 9.0, 9.0)
        atoms = {x['id']: x for x in editor.get_state()['atoms']}
        assert editor._estimate_bond_order(atoms[a], atoms[b], 1.5) == 2
        assert editor._estimate_bond_order(atoms[b], atoms[a], 1.5) == 2
        # Unbonded pair returns the neutral single order.
        assert editor._estimate_bond_order(atoms[a], atoms[c], 12.0) == 1

    def test_toggle_sync(self):
        editor = Ligand2DEditor()
        assert editor._sync_enabled is True
        editor.toggle_sync()
        assert editor._sync_enabled is False


class TestWaterNetworkAnalyzer:
    """Tests for water network analysis."""

    def test_detect_water_molecules_on_chains(self):
        analyzer = WaterNetworkAnalyzer(water_names=['HOH'])
        structure = Structure(id="test", source="local", title="Test")
        chain = Chain(id="A", name="Protein", is_polymer=True)
        chain.residues.append(Residue(
            id=100, name="HOH", chain_id="A", residue_number=100,
            atoms=[Atom(id=1, name="O", residue_name="HOH", residue_id=100,
                        chain_id="A", x=10, y=10, z=10, element="O")],
        ))
        structure.chains.append(chain)
        result = analyzer.analyze(structure)
        assert result['water_count'] == 1

    def test_detect_water_in_ligand_instances(self):
        """Waters parsed as non-polymer instances are found and split
        per residue."""
        analyzer = WaterNetworkAnalyzer(water_names=['HOH'])
        structure = Structure(id="test", source="local")
        ligand = Ligand(id="LHOH_A", name="HOH", residue_name="HOH",
                        atom_count=2)
        ligand.atoms = [
            Atom(id=1, name="O", residue_name="HOH", residue_id=501,
                 chain_id="A", x=0, y=0, z=0, element="O"),
            Atom(id=2, name="O", residue_name="HOH", residue_id=502,
                 chain_id="A", x=5, y=0, z=0, element="O"),
        ]
        structure.ligands.append(ligand)
        result = analyzer.analyze(structure)
        assert result['water_count'] == 2

    def test_network_respects_cutoff(self):
        analyzer = WaterNetworkAnalyzer(water_names=['HOH'],
                                        network_cutoff=3.0,
                                        contact_cutoff=4.0)
        waters = [
            {'residue_id': 100, 'chain_id': 'A', 'name': 'HOH',
             'atoms': [{'id': 1, 'name': 'O', 'element': 'O',
                        'x': 0, 'y': 0, 'z': 0}],
             'center': [0, 0, 0]},
            {'residue_id': 101, 'chain_id': 'A', 'name': 'HOH',
             'atoms': [{'id': 2, 'name': 'O', 'element': 'O',
                        'x': 2, 'y': 0, 'z': 0}],
             'center': [2, 0, 0]},
            {'residue_id': 102, 'chain_id': 'A', 'name': 'HOH',
             'atoms': [{'id': 3, 'name': 'O', 'element': 'O',
                        'x': 50, 'y': 0, 'z': 0}],
             'center': [50, 0, 0]},
        ]
        network = analyzer._build_network(waters)
        # Only the 2 A pair is within the 3 A cutoff.
        assert len(network) == 1
        assert network[0]['distance'] == 2.0

    def test_clusters_from_network(self):
        analyzer = WaterNetworkAnalyzer(water_names=['HOH'])
        network = [
            {'water1': {'chain': 'A', 'residue': 100},
             'water2': {'chain': 'A', 'residue': 101},
             'distance': 2.8, 'type': 'water-water'},
            {'water1': {'chain': 'A', 'residue': 101},
             'water2': {'chain': 'A', 'residue': 102},
             'distance': 2.9, 'type': 'water-water'},
        ]
        clusters = analyzer._find_clusters(network)
        assert len(clusters) == 1
        assert clusters[0]['size'] == 3

    def test_conservation_not_fabricated(self):
        analyzer = WaterNetworkAnalyzer(water_names=['HOH'])
        assert analyzer._analyze_conservation([]) == {}


class TestScriptingConsole:
    """Tests for the scripting console."""

    def test_execute_simple(self):
        console = ScriptingConsole()
        result = console.execute("x = 1 + 1")
        assert result['success'] is True
        assert console._variables['x'] == 2

    def test_execute_with_error(self):
        console = ScriptingConsole()
        result = console.execute("undefined_function_xyz()")
        assert result['success'] is False
        assert 'error' in result

    def test_execute_with_context(self):
        console = ScriptingConsole()
        result = console.execute(
            "n = len(structure.ligands)",
            context={'structure': Structure(id="s", source="local")})
        assert result['success'] is True
        assert console._variables['n'] == 0

    def test_execute_file(self):
        import tempfile
        import os
        console = ScriptingConsole()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py',
                                         delete=False) as f:
            f.write("y = 5 * 2")
            temp_path = f.name
        try:
            result = console.execute_file(temp_path)
            assert result['success'] is True
            assert console._variables['y'] == 10
        finally:
            os.unlink(temp_path)

    def test_history(self):
        console = ScriptingConsole()
        console.execute("a = 1")
        console.execute("b = 2")
        assert len(console.get_history()) == 2
        console.clear_history()
        assert console.get_history() == []


class TestResultComparator:
    """Tests for result comparison."""

    @staticmethod
    def _result(job_id, affinity, positions):
        return DockingResult(
            job_id=job_id,
            engine="vina",
            poses=[DockingPose(
                pose_id=1,
                affinity=affinity,
                ligand_name="LIG",
                atoms=[Atom(id=i + 1, name=f"C{i + 1}",
                            residue_name="LIG", residue_id=1,
                            chain_id="L", x=x, y=0, z=0, element="C")
                       for i, x in enumerate(positions)],
            )],
        )

    def test_compare_poses(self):
        comparator = ResultComparator()
        r1 = self._result("job1", -7.5, [0, 1])
        r2 = self._result("job2", -8.0, [0.5, 1.5])
        comparator.add_result(r1)
        comparator.add_result(r2)
        comparison = comparator.compare_poses(r1, r2)
        assert 'comparisons' in comparison

    def test_cluster_results(self):
        comparator = ResultComparator()
        comparator.add_result(self._result("job1", -7.5, [0]))
        comparator.add_result(self._result("job2", -8.0, [0.2]))
        clusters = comparator.cluster_results(rmsd_threshold=2.0)
        assert len(clusters) >= 1

    def test_consensus_pose(self):
        comparator = ResultComparator()
        comparator.add_result(self._result("job1", -7.5, [0]))
        consensus = comparator.consensus_pose()
        assert consensus is not None
        assert consensus.affinity == -7.5

    def test_clear_results(self):
        comparator = ResultComparator()
        comparator.add_result(self._result("job1", -7.5, [0]))
        comparator.clear_results()
        assert len(comparator._results) == 0


if __name__ == "__main__":
    pytest.main([__file__])
