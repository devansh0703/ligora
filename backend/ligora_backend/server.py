"""
Backend server for Ligora - handles IPC communication with the Tauri frontend.

Protocol: newline-delimited JSON over stdin/stdout.
- Requests:  {"command_id": "...", "type": "<command>", "payload": {...},
               "session_id": "..."}
- Responses: {"command_id": "...", "success": true/false, "data": {...},
               "error": "..."}
- Events (unsolicited): {"event": "job", "data": {...}}

Every command runs the real implementation. When an external engine or data
source is unavailable the command fails with a clear error - results are
never simulated.
"""

import sys
import io
import json
import contextlib
import threading
import time
import traceback
from pathlib import Path
from typing import Optional, Dict, Any, Callable, List

from .schemas import (
    CommandResponse,
    Structure,
    Ligand,
    Contact,
    Job,
    EngineType,
    DockingResult,
)
from .config import get_config, Config
from .workspace import WorkspaceManager, Session
from .parser import StructureParser
from .ligand import LigandResolver
from .contacts import ContactAnalyzer
from .cheminformatics import Cheminformatics
from .enrichment import EnrichmentClient
from .engines import EngineRegistry, VinaAdapter
from .search import SearchIndex
from .jobs import JobManager
from .export import ArtifactExporter
from .v2.water import WaterNetworkAnalyzer
from .v2.editors import Ligand2DEditor
from .v2.batch import BatchAnalyzer
from .v2.comparison import ResultComparator
from .v2.scripting import ScriptingConsole

import numpy as np


class BackendServer:
    """
    Backend server for Ligora.

    Handles communication with the Tauri frontend via stdin/stdout
    newline-delimited JSON. Routes commands to real implementations.
    """

    def __init__(self, config: Optional[Config] = None):
        self.config = config or get_config()
        self.workspace_manager = WorkspaceManager(self.config.workspace_dir)
        self.parser = StructureParser()
        self.ligand_resolver = LigandResolver()
        self.contact_analyzer = ContactAnalyzer()
        self.cheminformatics = Cheminformatics()
        self.enrichment_client = EnrichmentClient()
        self.engine_registry = EngineRegistry()
        self.job_manager = JobManager()
        self.exporter = ArtifactExporter(self.workspace_manager)
        self.batch_analyzer = BatchAnalyzer(
            workspace_manager=self.workspace_manager,
            parser=self.parser,
            ligand_resolver=self.ligand_resolver,
            contact_analyzer=self.contact_analyzer,
            enrichment_client=self.enrichment_client,
        )

        # Session storage
        self._sessions: Dict[str, Session] = {}
        self._current_session_id: Optional[str] = None

        # 2D editors per session
        self._editors: Dict[str, Ligand2DEditor] = {}

        self._handlers: Dict[str, Callable] = {
            'open_local_file': self._handle_open_local_file,
            'open_pdb_id': self._handle_open_pdb_id,
            'select_ligand': self._handle_select_ligand,
            'run_contact_analysis': self._handle_run_contact_analysis,
            'run_docking': self._handle_run_docking,
            'get_job': self._handle_get_job,
            'run_geometry_cleanup': self._handle_run_geometry_cleanup,
            'set_docking_box': self._handle_set_docking_box,
            'cancel_job': self._handle_cancel_job,
            'measure_distance': self._handle_measure_distance,
            'measure_angle': self._handle_measure_angle,
            'save_artifact': self._handle_save_artifact,
            'load_artifact': self._handle_load_artifact,
            'get_status': self._handle_get_status,
            'export_ligand_sdf': self._handle_export_ligand_sdf,
            'export_contacts_csv': self._handle_export_contacts_csv,
            'get_structure_file': self._handle_get_structure_file,
            'set_notes': self._handle_set_notes,
            # v2 features
            'analyze_water_network': self._handle_analyze_water_network,
            'get_ligand_2d': self._handle_get_ligand_2d,
            'editor_add_atom': self._handle_editor_add_atom,
            'editor_remove_atom': self._handle_editor_remove_atom,
            'editor_add_bond': self._handle_editor_add_bond,
            'editor_remove_bond': self._handle_editor_remove_bond,
            'editor_update_position': self._handle_editor_update_position,
            'editor_export_sdf': self._handle_editor_export_sdf,
            'batch_add': self._handle_batch_add,
            'batch_run': self._handle_batch_run,
            'batch_export': self._handle_batch_export,
            'run_script': self._handle_run_script,
            'list_jobs': self._handle_list_jobs,
            'compare_results': self._handle_compare_results,
            'search': self._handle_search,
            'search_index_refresh': self._handle_search_index_refresh,
        }

        self._running = False
        self._thread: Optional[threading.Thread] = None

    # ------------------------------------------------------------------
    # IPC loop
    # ------------------------------------------------------------------

    def start(self):
        """Start reading commands from stdin (blocking)."""
        self._running = True
        self._run()

    def stop(self):
        """Stop the server loop."""
        self._running = False

    def _run(self):
        """Main server loop: one JSON command per line on stdin."""
        while self._running:
            try:
                line = sys.stdin.readline()
                if not line:
                    break
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                except json.JSONDecodeError as e:
                    self._send_response(CommandResponse(
                        command_id='', success=False,
                        error=f'Invalid JSON: {e}'))
                    continue
                if isinstance(data, list):
                    for item in data:
                        self._process_command(item)
                else:
                    self._process_command(data)
            except Exception as e:
                self._send_response(CommandResponse(
                    command_id='', success=False, error=str(e)))

    def _process_command(self, data: Dict[str, Any]):
        command_type = data.get('type')
        payload = data.get('payload', {}) or {}
        session_id = data.get('session_id', '')
        command_id = data.get('command_id', '')

        handler = self._handlers.get(command_type)
        if handler:
            try:
                result = handler(payload, session_id)
                self._send_response(CommandResponse(
                    command_id=command_id, success=True, data=result))
            except Exception as e:
                traceback.print_exc(file=sys.stderr)
                self._send_response(CommandResponse(
                    command_id=command_id, success=False, error=str(e)))
        else:
            self._send_response(CommandResponse(
                command_id=command_id, success=False,
                error=f'Unknown command: {command_type}'))

    def _send_response(self, response: CommandResponse):
        message = {
            'command_id': response.command_id,
            'success': response.success,
            'data': response.data,
            'error': response.error,
        }
        print(json.dumps(message, default=str))
        sys.stdout.flush()

    def _send_event(self, event: str, data: Dict[str, Any]):
        message = {'event': event, 'data': data}
        print(json.dumps(message, default=str))
        sys.stdout.flush()

    # ------------------------------------------------------------------
    # Session helpers
    # ------------------------------------------------------------------

    def _require_session(self, session_id: str) -> Session:
        session = self._get_current_session(session_id)
        if not session:
            raise ValueError(
                f"No session: {session_id!r}. Open a structure first.")
        return session

    def _require_structure(self, session_id: str) -> Session:
        session = self._require_session(session_id)
        if not session.structure:
            raise ValueError("No structure loaded in session")
        return session

    def _create_session(self) -> Session:
        session = self.workspace_manager.create_session()
        self._sessions[session.id] = session
        self._current_session_id = session.id
        return session

    def _get_current_session(self, session_id: str) -> Optional[Session]:
        if session_id and session_id in self._sessions:
            return self._sessions[session_id]
        if self._current_session_id and self._current_session_id in self._sessions:
            return self._sessions[self._current_session_id]
        return None

    def _find_ligand(self, session: Session,
                     ligand_id: str) -> Optional[Ligand]:
        for lig in session.structure.ligands:
            if ligand_id in (lig.id, lig.name, lig.residue_name):
                return lig
        return None

    # ------------------------------------------------------------------
    # Structure loading
    # ------------------------------------------------------------------

    def _parse_structure_content(self, content: str, source_id: str,
                                 source: str, path: Path) -> Structure:
        suffix = path.suffix.lower()
        if suffix in ('.pdb', '.ent'):
            return self.parser.parse_pdb(content, source_id=source_id,
                                         source=source)
        # Default to mmCIF for .cif/.mcif and unknown suffixes.
        return self.parser.parse_mmcif(content, source_id=source_id,
                                       source=source)

    def _load_and_resolve(self, session: Session, content: str,
                          source_id: str, source: str, path: Path) -> Structure:
        structure = self._parse_structure_content(
            content, source_id, source, path)
        session.structure = structure
        session.file_path = str(path)
        for ligand in structure.ligands:
            self.ligand_resolver.resolve_ligand(ligand, structure.id)
        return structure

    def _handle_open_local_file(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Open a structure file from disk."""
        file_path = payload.get('file_path')
        if not file_path:
            raise ValueError("file_path is required")

        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        session = (self._sessions[session_id]
                   if session_id in self._sessions
                   else self._create_session())

        target = self.workspace_manager.copy_structure_to_workspace(
            session, path, path.name)
        content = path.read_text(encoding='utf-8', errors='replace')
        structure = self._load_and_resolve(
            session, content, path.stem, 'local', target)

        return {
            'structure': self._structure_to_dict(structure),
            'session_id': session.id,
            'workspace_path': str(session.get_workspace()),
        }

    def _handle_open_pdb_id(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Fetch a structure from RCSB, save it in the workspace, parse it."""
        pdb_id = (payload.get('pdb_id') or '').strip().upper()
        if not pdb_id:
            raise ValueError("pdb_id is required")

        session = (self._sessions[session_id]
                   if session_id in self._sessions
                   else self._create_session())

        # Download the authoritative mmCIF from files.rcsb.org.
        import requests
        url = f"{self.config.rcsb_files_url}/download/{pdb_id}.cif"
        try:
            response = requests.get(url, timeout=self.config.request_timeout)
            response.raise_for_status()
            content = response.text
        except requests.RequestException as e:
            raise RuntimeError(f"Failed to download {pdb_id} from RCSB: {e}")

        workspace_path = session.get_workspace() / f'{pdb_id}.cif'
        workspace_path.write_text(content, encoding='utf-8')

        structure = self._load_and_resolve(
            session, content, pdb_id, 'rcsb', workspace_path)

        return {
            'structure': self._structure_to_dict(structure),
            'session_id': session.id,
            'workspace_path': str(session.get_workspace()),
        }

    # ------------------------------------------------------------------
    # Ligand selection + contact analysis
    # ------------------------------------------------------------------

    def _handle_select_ligand(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        ligand_id = payload.get('ligand_id')
        if not ligand_id:
            raise ValueError("ligand_id is required")

        session = self._require_structure(session_id)
        ligand = self._find_ligand(session, ligand_id)
        if not ligand:
            raise ValueError(f"Ligand not found: {ligand_id}")

        # An honest warning (not a refusal): water or an ion can be selected
        # for inspection, but it is not a drug-like ligand. The classification
        # comes from the CCD pdbx_type attached during resolution.
        warning = None
        if ligand.classification_hint in ('HETAS', 'HETAI'):
            warning = (
                f"{ligand.residue_name} is classified as "
                f"{'solvent' if ligand.classification_hint == 'HETAS' else 'ion'} "
                f"by the RCSB CCD, not a drug-like ligand")

        session.selected_ligand_id = ligand.id
        enriched = self.ligand_resolver.resolve_ligand(
            ligand, session.structure.id)

        return {
            'ligand': self._ligand_to_dict(enriched),
            'warning': warning,
            'session_id': session.id,
        }

    def _handle_run_contact_analysis(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_structure(session_id)
        if not session.selected_ligand_id:
            raise ValueError("No ligand selected")

        ligand = self._find_ligand(session, session.selected_ligand_id)
        if not ligand:
            raise ValueError(
                f"Selected ligand not found: {session.selected_ligand_id}")

        # Real PLIP analysis. When PLIP is missing, plip_available is False
        # and is surfaced to the UI instead of an empty-contact success.
        contacts, plip_available = self.contact_analyzer.analyze_contacts(
            session.structure, ligand)
        pocket = self.contact_analyzer.compute_binding_pocket(
            session.structure, ligand)

        enriched = self.ligand_resolver.resolve_ligand(
            ligand, session.structure.id)

        evidence = self.enrichment_client.get_all_evidence(
            {
                'name': enriched.name,
                'residue_name': enriched.residue_name,
                'pubchem_cid': enriched.pubchem_cid,
                'chembl_id': enriched.chembl_id,
                'pdbbind_affinity': enriched.pdbbind_affinity,
                'classification_hint': enriched.classification_hint,
            },
            pdb_id=session.structure.id,
        )
        # Evidence chain recorded by the resolver itself (CCD/PubChem/ChEMBL).
        evidence.extend([
            {'source': e.get('source'), 'field': 'resolution',
             'value': e.get('fields'), 'url': e.get('url')}
            for e in self.ligand_resolver.get_last_evidence()
        ])

        self.contact_analyzer.export_contacts_csv(
            contacts, session.get_contacts_path())

        return {
            'contacts': [self._contact_to_dict(c) for c in contacts],
            'contact_count': len(contacts),
            'plip_available': plip_available,
            'pocket': pocket,
            'ligand': self._ligand_to_dict(enriched),
            'evidence': evidence,
        }

    # ------------------------------------------------------------------
    # Docking (real Vina via JobManager)
    # ------------------------------------------------------------------

    def _handle_run_docking(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_structure(session_id)
        if not session.selected_ligand_id:
            raise ValueError("No ligand selected")

        ligand = self._find_ligand(session, session.selected_ligand_id)
        if not ligand:
            raise ValueError(
                f"Selected ligand not found: {session.selected_ligand_id}")

        vina: VinaAdapter = self.engine_registry.get_adapter(EngineType.VINA)
        if not vina.is_available():
            raise RuntimeError(
                "AutoDock Vina is not installed. Install it (e.g. "
                "'apt install autodock-vina') to run docking - docking "
                "results are never simulated.")

        # Preparation can be slow (Open Babel) - do it before the job so
        # errors surface synchronously to the user.
        output_dir = session.get_job_output_dir(f"vina_{int(time.time())}")
        box_center = payload.get('box_center')
        box_size = payload.get('box_size')
        session.docking_box = {
            'center': box_center, 'size': box_size} if box_center else {}
        input_files = vina.prepare(
            session.structure, ligand, output_dir,
            box_center=box_center, box_size=box_size)

        job = self.job_manager.create_job(
            EngineType.VINA,
            {k: v for k, v in payload.items() if k != 'session_id'},
            session.id,
        )
        session.add_job(job)

        def run_fn(job: Job, params: Dict[str, Any]):
            result = vina.run(input_files, params)
            self.job_manager.update_progress(
                job.id, 1.0, f"{len(result['poses'])} poses, "
                f"{result['runtime_seconds']:.1f}s")
            return result

        def on_event(job_id: str, event: str, job: Job):
            self._send_event('job', {
                'job_id': job.id,
                'status': job.status.value,
                'progress': job.progress,
                'error': job.error,
                'session_id': session.id,
            })

        self.job_manager.register_callback(job.id, on_event)
        self.job_manager.start_job(job, run_fn)

        return {
            'job_id': job.id,
            'status': 'running',
            'message': 'Docking started',
        }

    def _handle_get_job(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        job_id = payload.get('job_id')
        if not job_id:
            raise ValueError("job_id is required")
        job = self.job_manager.get_job(job_id)
        if not job:
            raise ValueError(f"Job not found: {job_id}")

        data: Dict[str, Any] = {
            'job_id': job.id,
            'type': job.type.value,
            'status': job.status.value,
            'progress': job.progress,
            'error': job.error,
        }
        if isinstance(job.result, dict):
            result = dict(job.result)
            result['poses'] = [
                self._pose_to_dict(p) for p in result.get('poses', [])]
            data['result'] = result
        return data

    # ------------------------------------------------------------------
    # Geometry cleanup (real Open Babel minimization)
    # ------------------------------------------------------------------

    def _handle_run_geometry_cleanup(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_structure(session_id)
        if not session.selected_ligand_id:
            raise ValueError("No ligand selected")

        ligand = self._find_ligand(session, session.selected_ligand_id)
        if not ligand:
            raise ValueError(
                f"Selected ligand not found: {session.selected_ligand_id}")

        result = self.cheminformatics.cleanup_geometry(ligand.atoms)
        if result.get('status') != 'ok':
            raise RuntimeError(result.get(
                'message', 'Geometry cleanup failed'))

        atoms = result.pop('atoms')
        return {
            'status': 'ok',
            'force_field': result.get('force_field'),
            'steps': result.get('steps'),
            'rmsd': result.get('rmsd'),
            'energy_kcal_mol': result.get('energy_kcal_mol'),
            'atoms': [self._atom_to_dict(a) for a in atoms],
        }

    # ------------------------------------------------------------------
    # Docking box
    # ------------------------------------------------------------------

    def _handle_set_docking_box(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        required = ('center_x', 'center_y', 'center_z',
                    'size_x', 'size_y', 'size_z')
        missing = [k for k in required if payload.get(k) is None]
        if missing:
            raise ValueError(f"Missing box parameters: {missing}")
        try:
            box = {
                'center': [float(payload['center_x']),
                           float(payload['center_y']),
                           float(payload['center_z'])],
                'size': [float(payload['size_x']), float(payload['size_y']),
                         float(payload['size_z'])],
            }
        except (TypeError, ValueError):
            raise ValueError("Box parameters must be numbers")
        if any(s <= 0 for s in box['size']):
            raise ValueError("Box sizes must be positive")
        session.docking_box = box
        return {'status': 'ok', 'box': box}

    # ------------------------------------------------------------------
    # Jobs
    # ------------------------------------------------------------------

    def _handle_cancel_job(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        job_id = payload.get('job_id')
        if not job_id:
            raise ValueError("job_id is required")
        cancelled = self.job_manager.cancel_job(job_id)
        return {'cancelled': cancelled}

    # ------------------------------------------------------------------
    # Measurements (real coordinates from the loaded structure)
    # ------------------------------------------------------------------

    def _find_atom(self, structure: Structure, atom_ref) -> Optional[Any]:
        """Find an atom by identifier.

        Accepts {"chain": "A", "residue_id": 12, "residue_name": "GLU",
        "name": "OE1"} objects, "chain:resid:name" strings, and
        "chain:resid:resname:name" strings.
        """
        if isinstance(atom_ref, dict):
            chain_id = str(atom_ref.get('chain', ''))
            res_id = atom_ref.get('residue_id')
            atom_name = atom_ref.get('name')
            for chain in structure.chains:
                if chain.id != chain_id:
                    continue
                for residue in chain.residues:
                    if res_id is not None and \
                            residue.id != res_id and \
                            str(residue.id) != str(res_id):
                        continue
                    for atom in residue.atoms:
                        if atom_name is None or atom.name == atom_name:
                            return atom
            # Ligand instances are not chain residues; search them too so
            # measurements can reference ligand atoms. The instance's own
            # atoms carry the residue id/chain of that instance.
            for ligand in structure.ligands:
                if not ligand.atoms:
                    continue
                if chain_id and \
                        (ligand.atoms[0].chain_id or '') != chain_id:
                    continue
                if res_id is not None:
                    inst_res = ligand.atoms[0].residue_id
                    if inst_res != res_id and str(inst_res) != str(res_id):
                        continue
                if atom_name is None:
                    return ligand.atoms[0]
                for atom in ligand.atoms:
                    if atom.name == atom_name:
                        return atom
            return None

        atom_id = str(atom_ref)
        parts = atom_id.split(':')
        if len(parts) >= 3:
            chain_id, res_id_s, atom_name = parts[0], parts[-2], parts[-1]
            res_name = parts[2] if len(parts) >= 4 else None
            try:
                res_id = int(res_id_s)
            except ValueError:
                return None
            for chain in structure.chains:
                if chain.id != chain_id:
                    continue
                for residue in chain.residues:
                    if residue.id != res_id:
                        continue
                    if res_name and residue.name != res_name:
                        continue
                    for atom in residue.atoms:
                        if atom.name == atom_name:
                            return atom
        return None

    def _handle_measure_distance(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        atom1 = payload.get('atom1')
        atom2 = payload.get('atom2')
        if not atom1 or not atom2:
            raise ValueError("atom1 and atom2 are required")

        session = self._require_structure(session_id)
        a1 = self._find_atom(session.structure, atom1)
        a2 = self._find_atom(session.structure, atom2)
        if not a1 or not a2:
            raise ValueError("Atoms not found")

        pos1 = np.array([a1.x, a1.y, a1.z])
        pos2 = np.array([a2.x, a2.y, a2.z])
        return {
            'distance': round(float(np.linalg.norm(pos1 - pos2)), 4),
            'atom1': self._atom_to_dict(a1),
            'atom2': self._atom_to_dict(a2),
        }

    def _handle_measure_angle(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        atom1, atom2, atom3 = (payload.get('atom1'), payload.get('atom2'),
                               payload.get('atom3'))
        if not all([atom1, atom2, atom3]):
            raise ValueError("atom1, atom2, and atom3 are required")

        session = self._require_structure(session_id)
        a1 = self._find_atom(session.structure, atom1)
        a2 = self._find_atom(session.structure, atom2)
        a3 = self._find_atom(session.structure, atom3)
        if not all([a1, a2, a3]):
            raise ValueError("Atoms not found")

        v1 = np.array([a1.x, a1.y, a1.z]) - np.array([a2.x, a2.y, a2.z])
        v2 = np.array([a3.x, a3.y, a3.z]) - np.array([a2.x, a2.y, a2.z])
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 == 0 or n2 == 0:
            raise ValueError("Degenerate angle: coincident atoms")
        cos = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
        return {
            'angle': round(float(np.degrees(np.arccos(cos))), 3),
            'atom1': self._atom_to_dict(a1),
            'atom2': self._atom_to_dict(a2),
            'atom3': self._atom_to_dict(a3),
        }

    # ------------------------------------------------------------------
    # Artifacts
    # ------------------------------------------------------------------

    def _handle_save_artifact(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        if not session.structure:
            raise ValueError("No structure loaded in session")

        ligand = (self._find_ligand(session, session.selected_ligand_id)
                  if session.selected_ligand_id else None)
        if not ligand:
            raise ValueError("No ligand selected; nothing to save")

        contacts = []
        contacts_csv = session.get_contacts_path()
        if contacts_csv.exists():
            # Re-read the contacts exported during analysis.
            import csv as _csv
            with open(contacts_csv, newline='', encoding='utf-8') as f:
                contacts = list(_csv.DictReader(f))
            for i, c in enumerate(contacts):
                c['id'] = i

        evidence = self.enrichment_client.get_all_evidence(
            {
                'name': ligand.name,
                'residue_name': ligand.residue_name,
                'pubchem_cid': ligand.pubchem_cid,
                'chembl_id': ligand.chembl_id,
                'pdbbind_affinity': ligand.pdbbind_affinity,
                'classification_hint': ligand.classification_hint,
            },
            pdb_id=session.structure.id,
        )

        # Export the ligand SDF with real perceived bonds next to the summary.
        ligand_sdf_path = session.get_workspace() / \
            f"ligand_{ligand.residue_name}.sdf"
        try:
            self.exporter.export_ligand_sdf(ligand, ligand_sdf_path)
        except (ValueError, OSError):
            ligand_sdf_path = None

        summary = self.exporter.export_analysis_summary(
            structure=session.structure,
            ligand=ligand,
            contacts=contacts,
            evidence=[e.__dict__ if hasattr(e, '__dict__') else e
                      for e in evidence],
            job_results=[
                {'id': j.id, 'type': j.type.value, 'status': j.status.value}
                for j in session.jobs.values()
            ],
            notes=session.notes,
        )

        out_dir = session.get_workspace()
        summary_path = out_dir / 'analysis_summary.json'
        self.exporter.save_analysis_summary(
            session.id, summary, summary_path)
        if ligand_sdf_path is not None:
            summary.ligand_export_path = str(ligand_sdf_path)
            self.exporter.save_analysis_summary(
                session.id, summary, summary_path)
        contacts_path = out_dir / 'contacts.csv'
        if contacts_path.exists():
            summary.contact_export_path = str(contacts_path)

        return {
            'path': str(summary_path),
            'ligand_sdf': (str(ligand_sdf_path)
                           if ligand_sdf_path else None),
            'contacts_csv': (str(contacts_path)
                             if contacts_path.exists() else None),
            'session_id': session.id,
        }

    def _handle_load_artifact(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        path = payload.get('path')
        if not path:
            raise ValueError("path is required")
        artifact_path = Path(path)
        if not artifact_path.exists():
            raise FileNotFoundError(f"Artifact not found: {path}")

        try:
            with open(artifact_path, encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Artifact is not valid JSON: {e}")

        required = ('structure_id', 'ligand_name', 'contacts')
        missing = [k for k in required if k not in data]
        if missing:
            raise ValueError(
                f"Artifact missing fields: {missing}; not a Ligora "
                f"analysis summary")

        return {'summary': data}

    # ------------------------------------------------------------------
    # Exports
    # ------------------------------------------------------------------

    def _handle_export_ligand_sdf(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_structure(session_id)
        if not session.selected_ligand_id:
            raise ValueError("No ligand selected")
        ligand = self._find_ligand(session, session.selected_ligand_id)
        if not ligand:
            raise ValueError(
                f"Selected ligand not found: {session.selected_ligand_id}")

        out_path = session.get_workspace() / \
            f"ligand_{ligand.residue_name}.sdf"
        self.exporter.export_ligand_sdf(ligand, out_path)
        return {'path': str(out_path), 'format': 'sdf'}

    def _handle_export_contacts_csv(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        contacts_path = session.get_contacts_path()
        if not contacts_path.exists():
            raise ValueError(
                "No contacts exported yet; run contact analysis first")
        return {'path': str(contacts_path), 'format': 'csv'}

    # ------------------------------------------------------------------
    # Structure file access (for the 3D viewer)
    # ------------------------------------------------------------------

    def _handle_get_structure_file(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Return the loaded structure file content for the viewer."""
        session = self._require_structure(session_id)
        path = session.get_workspace() / f"{session.structure.id}.cif"
        if not path.exists():
            # Local files are copied under their original name.
            for candidate in session.get_workspace().glob('*'):
                if candidate.suffix.lower() in ('.cif', '.mcif', '.pdb',
                                                '.ent'):
                    path = candidate
                    break
        if not path.exists():
            raise FileNotFoundError(
                "Structure file not found in workspace")
        content = path.read_text(encoding='utf-8', errors='replace')
        fmt = ('pdb' if path.suffix.lower() in ('.pdb', '.ent')
               else 'cif')
        return {
            'path': str(path),
            'format': fmt,
            'content': content,
        }

    def _handle_set_notes(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Store the session notes (included in saved artifacts)."""
        session = self._require_session(session_id)
        notes = payload.get('notes')
        if notes is None:
            raise ValueError("notes is required")
        session.notes = str(notes)
        return {'status': 'ok'}

    # ------------------------------------------------------------------
    # Status
    # ------------------------------------------------------------------

    def _handle_get_status(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._get_current_session(session_id)

        job_status = self.job_manager.get_status_summary()
        engine_status = self.engine_registry.get_engine_status()
        data_sources_status = self.enrichment_client.health_check()

        session_payload = None
        if session:
            session_payload = {
                'id': session.id,
                'structure_id': (session.structure.id
                                 if session.structure else None),
                'selected_ligand_id': session.selected_ligand_id,
                'notes': session.notes,
            }

        vina = self.engine_registry.get_adapter(EngineType.VINA)
        engines_out = {name: {'available': info.get('available', False)}
                       for name, info in engine_status.items()}
        engines_out['plip'] = {
            'available': self.contact_analyzer.is_plip_available()}
        engines_out['openbabel'] = {
            'available': self.cheminformatics.is_obabel_available()}
        engines_out['rdkit'] = {'available': True}

        return {
            'session': session_payload,
            'jobs': job_status,
            'engines': engines_out,
            'data_sources': data_sources_status,
            'vina_preparation_available': vina.preparation_available(),
        }

    # ------------------------------------------------------------------
    # v2: water network analysis
    # ------------------------------------------------------------------

    def _handle_analyze_water_network(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_structure(session_id)

        # Water names come from the CCD classification of the structure's own
        # components (pdbx_type == HETAS = solvent), not a local list.
        water_names = set()
        for ligand in session.structure.ligands:
            if ligand.classification_hint == 'HETAS':
                water_names.add(ligand.residue_name)
        # Waters parsed from the structure but not resolved as ligand
        # instances: residue-level scan for component names whose atoms are
        # a single oxygen named 'O' (the structure file's own data).
        for chain in session.structure.chains:
            if chain.is_polymer:
                continue
            for residue in chain.residues:
                atoms = residue.atoms
                if (len(atoms) == 1 and atoms[0].name == 'O' and
                        (atoms[0].element or '').upper() == 'O'):
                    water_names.add(residue.name)

        analyzer = WaterNetworkAnalyzer(
            water_names=sorted(water_names),
            network_cutoff=payload.get(
                'network_cutoff', self.config.water_network_cutoff),
            contact_cutoff=payload.get(
                'contact_cutoff', self.config.water_contact_cutoff),
        )
        return analyzer.analyze(session.structure)

    # ------------------------------------------------------------------
    # v2: 2D ligand editor
    # ------------------------------------------------------------------

    def _get_editor(self, session: Session) -> Ligand2DEditor:
        editor = self._editors.get(session.id)
        if editor is None:
            editor = Ligand2DEditor()
            self._editors[session.id] = editor
        return editor

    def _require_editor_ligand(self, session: Session) -> Ligand:
        if not session.selected_ligand_id:
            raise ValueError("No ligand selected")
        ligand = self._find_ligand(session, session.selected_ligand_id)
        if not ligand:
            raise ValueError(
                f"Selected ligand not found: {session.selected_ligand_id}")
        return ligand

    def _handle_get_ligand_2d(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_structure(session_id)
        ligand = self._require_editor_ligand(session)
        editor = self._get_editor(session)
        state = editor.load_ligand(ligand)
        # Attach the resolved SMILES so the UI can render a proper 2D layout
        # when the ligand chemistry is known.
        state['smiles'] = ligand.smiles
        state['residue_name'] = ligand.residue_name
        return state

    def _handle_editor_add_atom(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        self._require_editor_ligand(session)
        editor = self._get_editor(session)
        element = payload.get('element')
        if not element or not isinstance(element, str) or \
                len(element) > 2 or not element.isalpha():
            raise ValueError("element must be a 1-2 letter symbol")
        try:
            atom_id = editor.add_atom(
                element.capitalize(), float(payload['x']), float(payload['y']))
        except KeyError as e:
            raise ValueError(f"Missing field: {e}")
        return {'atom_id': atom_id, 'state': editor.get_state()}

    def _handle_editor_remove_atom(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        editor = self._get_editor(session)
        editor.remove_atom(int(payload['atom_id']))
        return {'state': editor.get_state()}

    def _handle_editor_add_bond(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        editor = self._get_editor(session)
        order = int(payload.get('order', 1))
        if order not in (1, 2, 3):
            raise ValueError("order must be 1, 2 or 3")
        editor.add_bond(int(payload['from_atom']), int(payload['to_atom']),
                        order)
        return {'state': editor.get_state()}

    def _handle_editor_remove_bond(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        editor = self._get_editor(session)
        editor.remove_bond(int(payload['from_atom']), int(payload['to_atom']))
        return {'state': editor.get_state()}

    def _handle_editor_update_position(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        editor = self._get_editor(session)
        state = editor.update_atom_position(
            int(payload['atom_id']), float(payload['x']), float(payload['y']))
        return {'state': state}

    def _handle_editor_export_sdf(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        session = self._require_session(session_id)
        editor = self._get_editor(session)
        try:
            sdf = editor.export_sdf()
        except ValueError as e:
            raise ValueError(
                f"Edited molecule is not chemically valid: {e}")
        out_path = session.get_workspace() / 'edited_ligand.sdf'
        out_path.write_text(sdf, encoding='utf-8')
        return {'path': str(out_path)}

    # ------------------------------------------------------------------
    # v2: batch analysis
    # ------------------------------------------------------------------

    def _handle_batch_add(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        sources = payload.get('sources') or []
        if isinstance(sources, str):
            sources = [sources]
        if not sources:
            raise ValueError("sources is required")
        batch_id = self.batch_analyzer.create_batch()
        added = []
        for source in sources:
            task = self.batch_analyzer.add_structure_to_batch(
                batch_id, source,
                source_type=payload.get('source_type', 'pdb_id'))
            if task:
                added.append({'task_id': task.id, 'source': source})
        if not added:
            raise RuntimeError("No tasks could be added to the batch")
        return {'batch_id': batch_id, 'tasks': added}

    def _handle_batch_run(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        batch_id = payload.get('batch_id')
        if not batch_id:
            raise ValueError("batch_id is required")

        progress_rows: List[Dict[str, Any]] = []

        def on_progress(progress: float, task):
            progress_rows.append({
                'task_id': task.id, 'progress': progress,
                'status': task.status,
            })

        results = self.batch_analyzer.run_batch(
            batch_id, on_progress=on_progress)
        summary = {
            'batch_id': batch_id,
            'total': len(results),
            'succeeded': sum(1 for r in results if r.success),
            'failed': sum(1 for r in results if not r.success),
            'results': [
                {
                    'structure_id': r.structure_id,
                    'success': r.success,
                    'error': r.error,
                    'contact_count': (r.summary.contact_count
                                      if r.summary else None),
                    'ligand_name': (r.summary.ligand_name
                                    if r.summary else None),
                }
                for r in results
            ],
        }
        self._send_event('batch_progress', {
            'batch_id': batch_id, 'tasks': progress_rows})
        return summary

    def _handle_batch_export(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        batch_id = payload.get('batch_id')
        if not batch_id:
            raise ValueError("batch_id is required")
        out_dir = (Path(payload['output_dir'])
                   if payload.get('output_dir') else None)
        path = self.batch_analyzer.export_batch_results(
            batch_id, output_dir=out_dir)
        return {'path': str(path)}

    # ------------------------------------------------------------------
    # v2: scripting console
    # ------------------------------------------------------------------

    def _handle_run_script(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        code = payload.get('code')
        if not code:
            raise ValueError("code is required")
        session = self._get_current_session(session_id)
        console = self._scripting_console()
        # Redirect the script's stdout (print etc.) into a buffer: user
        # output must never leak into the IPC stream on this process's
        # real stdout, and the console should display what was printed.
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            console.execute(code, context=self._scripting_context(session))
        return {
            'history': console.get_history()[-1:],
            'variables': list(console.get_variables().keys()),
            'output': buffer.getvalue(),
        }

    def _scripting_console(self):
        if not hasattr(self, '_script_console'):
            self._script_console = ScriptingConsole()
        return self._script_console

    # ------------------------------------------------------------------
    # v2: result comparison
    # ------------------------------------------------------------------

    def _handle_list_jobs(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """List all jobs with their results (for comparison selection)."""
        jobs = []
        for job in self.job_manager.get_all_jobs():
            entry = {
                'job_id': job.id,
                'type': job.type.value,
                'status': job.status.value,
                'progress': job.progress,
                'error': job.error,
            }
            if isinstance(job.result, dict):
                entry['pose_count'] = len(job.result.get('poses', []))
                entry['best_affinity'] = min(
                    (p.affinity for p in job.result.get('poses', [])
                     if p.affinity == p.affinity),  # skip NaN
                    default=None)
            jobs.append(entry)
        return {'jobs': jobs}

    # ------------------------------------------------------------------
    # BM25 search over live-fetched documents
    # ------------------------------------------------------------------

    def _search_index(self) -> SearchIndex:
        if not hasattr(self, '_search_idx'):
            self._search_idx = SearchIndex()
        return self._search_idx

    def _handle_search(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        query = (payload.get('query') or '').strip()
        if not query:
            raise ValueError("query is required")
        index = self._search_index()
        result = index.search(
            query,
            scope=payload.get('scope') or None,
            limit=int(payload.get('limit', 20)),
        )
        # Auto-refresh hint: if nothing is indexed yet, tell the UI.
        result['needs_refresh'] = result.get('documents_indexed', 0) == 0
        return result

    def _handle_search_index_refresh(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        index = self._search_index()
        # Seed from the current session's own real data first.
        session = self._get_current_session(session_id)
        pdb_ids = list(payload.get('pdb_ids') or [])
        chem_names = list(payload.get('chem_names') or [])
        if session and session.structure:
            pdb_ids.append(session.structure.id)
            for lig in session.structure.ligands:
                chem_names.append(lig.residue_name)
        if not pdb_ids and not chem_names:
            raise ValueError(
                "Nothing to index: open a structure or pass pdb_ids/"
                "chem_names to fetch from live sources")
        return index.refresh(pdb_ids=pdb_ids, chem_names=chem_names)

    def _handle_compare_results(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Compare two completed docking jobs (pose RMSD / affinity)."""
        job_a = payload.get('job_id_a')
        job_b = payload.get('job_id_b')
        if not job_a or not job_b:
            raise ValueError("job_id_a and job_id_b are required")

        def _result(job_id: str):
            job = self.job_manager.get_job(job_id)
            if not job:
                raise ValueError(f"Job not found: {job_id}")
            if job.status.value != 'completed' or \
                    not isinstance(job.result, dict):
                raise ValueError(
                    f"Job {job_id} has no completed docking result")
            poses = job.result.get('poses')
            if not poses:
                raise ValueError(f"Job {job_id} has no poses to compare")
            return DockingResult(
                job_id=job.id, engine=job.type, poses=poses)

        comparator = ResultComparator()
        comparison = comparator.compare_poses(_result(job_a), _result(job_b))
        return {
            'comparison': comparison,
            'best_matches': comparison.get('best_matches', []),
            'avg_rmsd': comparison.get('avg_rmsd'),
        }

    @staticmethod
    def _scripting_context(session: Optional[Session]) -> Dict[str, Any]:
        """Variables exposed to user scripts (real objects, no wrappers)."""
        if session is None:
            return {}
        return {
            'session_id': session.id,
            'structure': session.structure,
            'selected_ligand_id': session.selected_ligand_id,
        }

    # ------------------------------------------------------------------
    # Serialization helpers
    # ------------------------------------------------------------------

    def _structure_to_dict(self, structure: Structure) -> Dict[str, Any]:
        return {
            'id': structure.id,
            'title': structure.title,
            'source': structure.source,
            'file_format': structure.file_format,
            'resolution': structure.resolution,
            'experiment_type': structure.experiment_type,
            'chains': [
                {
                    'id': c.id,
                    'name': c.name,
                    'is_polymer': c.is_polymer,
                    'residue_count': len(c.residues),
                    'atom_count': sum(len(r.atoms) for r in c.residues),
                }
                for c in structure.chains
            ],
            'ligands': [self._ligand_to_dict(lig) for lig in structure.ligands],
        }

    def _ligand_to_dict(self, ligand: Ligand) -> Dict[str, Any]:
        return {
            'id': ligand.id,
            'name': ligand.name,
            'residue_name': ligand.residue_name,
            'formula': ligand.formula,
            'molecular_weight': ligand.molecular_weight,
            'atom_count': ligand.atom_count or len(ligand.atoms),
            'smiles': ligand.smiles,
            'inchi_key': ligand.inchi_key,
            'iupac_name': getattr(ligand, 'iupac_name', None),
            'pubchem_cid': ligand.pubchem_cid,
            'pubchem_name': ligand.pubchem_name,
            'chembl_id': ligand.chembl_id,
            'pdbbind_affinity': ligand.pdbbind_affinity,
            'classification_hint': ligand.classification_hint,
            'resolution_status': (
                ligand.resolution_status.value
                if hasattr(ligand.resolution_status, 'value')
                else str(ligand.resolution_status)),
            'has_2d_structure': ligand.has_2d_structure,
            'atoms': [self._atom_to_dict(a) for a in ligand.atoms[:500]],
        }

    def _contact_to_dict(self, contact: Contact) -> Dict[str, Any]:
        return {
            'id': contact.id,
            'ligand_atom': contact.ligand_atom,
            'ligand_residue_name': contact.ligand_residue_name,
            'ligand_residue_id': contact.ligand_residue_id,
            'ligand_chain_id': contact.ligand_chain_id,
            'protein_residue_name': contact.protein_residue_name,
            'protein_residue_id': contact.protein_residue_id,
            'protein_chain_id': contact.protein_chain_id,
            'protein_atom': contact.protein_atom,
            'distance': contact.distance,
            'contact_type': (
                contact.contact_type.value
                if hasattr(contact.contact_type, 'value')
                else str(contact.contact_type)),
            'angle': contact.angle,
            'is_water_mediated': contact.is_water_mediated,
            'description': contact.description,
        }

    def _pose_to_dict(self, pose) -> Dict[str, Any]:
        return {
            'pose_id': pose.pose_id,
            'affinity': pose.affinity,
            'ligand_name': pose.ligand_name,
            'atoms': [self._atom_to_dict(a) for a in pose.atoms],
        }

    def _atom_to_dict(self, atom: Any) -> Dict[str, Any]:
        return {
            'id': getattr(atom, 'id', 0),
            'name': getattr(atom, 'name', ''),
            'residue_name': getattr(atom, 'residue_name', ''),
            'residue_id': getattr(atom, 'residue_id', 0),
            'chain_id': getattr(atom, 'chain_id', ''),
            'x': getattr(atom, 'x', 0),
            'y': getattr(atom, 'y', 0),
            'z': getattr(atom, 'z', 0),
            'element': getattr(atom, 'element', None),
            'b_factor': getattr(atom, 'b_factor', 0),
        }


def start_server(config: Optional[Config] = None):
    """Start the backend server (blocks on stdin)."""
    server = BackendServer(config)
    server.start()
    return server


if __name__ == '__main__':
    print("Ligora backend ready", file=sys.stderr)
    sys.stderr.flush()
    try:
        start_server()
    except KeyboardInterrupt:
        pass
