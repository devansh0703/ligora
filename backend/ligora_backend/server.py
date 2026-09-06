"""
Backend server for Ligora - handles IPC communication with the Tauri frontend.

Provides:
- JSON-RPC style command handling
- WebSocket or stdin/stdout communication
- Command routing to appropriate handlers
- Session and job management
"""

import sys
import json
import threading
import time
import os
from pathlib import Path
from typing import Optional, Dict, Any, Callable, List
from datetime import datetime, timezone

from .schemas import (
    Command,
    CommandResponse,
    Structure,
    Ligand,
    Contact,
    Job,
    AnalysisSummary,
    CommandType,
    JobStatus,
    EngineType,
)
from .config import get_config, Config
from .workspace import WorkspaceManager, Session
from .parser import StructureParser
from .ligand import LigandResolver
from .contacts import ContactAnalyzer
from .cheminformatics import Cheminformatics
from .enrichment import EnrichmentClient
from .engines import EngineRegistry, VinaAdapter
from .jobs import JobManager
from .export import ArtifactExporter

import numpy as np


class BackendServer:
    """
    Backend server for Ligora.

    Handles communication with the Tauri frontend via stdin/stdout
    or a socket connection. Routes commands to appropriate handlers.
    """

    def __init__(self, config: Optional[Config] = None):
        """Initialize the backend server."""
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

        # Session storage
        self._sessions: Dict[str, Session] = {}
        self._current_session_id: Optional[str] = None

        # Command handlers
        self._handlers: Dict[str, Callable] = {
            'open_local_file': self._handle_open_local_file,
            'open_pdb_id': self._handle_open_pdb_id,
            'select_ligand': self._handle_select_ligand,
            'run_contact_analysis': self._handle_run_contact_analysis,
            'run_docking': self._handle_run_docking,
            'run_geometry_cleanup': self._handle_run_geometry_cleanup,
            'set_docking_box': self._handle_set_docking_box,
            'cancel_job': self._handle_cancel_job,
            'measure_distance': self._handle_measure_distance,
            'measure_angle': self._handle_measure_angle,
            'save_artifact': self._handle_save_artifact,
            'load_artifact': self._handle_load_artifact,
            'get_status': self._handle_get_status,
            'export_scene_image': self._handle_export_scene_image,
        }

        self._running = False
        self._thread: Optional[threading.Thread] = None

    def start(self):
        """Start the backend server."""
        self._running = True
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self):
        """Stop the backend server."""
        self._running = False
        if self._thread:
            self._thread.join(timeout=5)

    def _run(self):
        """Main server loop."""
        # Read commands from stdin
        buffer = ""

        while self._running:
            try:
                # Read line from stdin
                line = sys.stdin.readline()
                if not line:
                    break

                buffer += line.strip()

                # Process complete JSON commands
                if buffer.startswith('{') and buffer.endswith('}'):
                    try:
                        data = json.loads(buffer)
                        self._process_command(data)
                        buffer = ""
                    except json.JSONDecodeError:
                        # Incomplete JSON, continue reading
                        pass
                elif '\n' in buffer:
                    # Multiple lines, try to parse each
                    parts = buffer.split('\n')
                    for part in parts:
                        part = part.strip()
                        if part.startswith('{') and part.endswith('}'):
                            try:
                                data = json.loads(part)
                                self._process_command(data)
                            except json.JSONDecodeError:
                                pass
                    buffer = ""

            except Exception as e:
                self._send_response({
                    'success': False,
                    'error': str(e),
                })

    def _process_command(self, data: Dict[str, Any]):
        """Process a command from the frontend."""
        command_type = data.get('type')
        payload = data.get('payload', {})
        session_id = data.get('session_id', '')
        command_id = data.get('command_id', '')

        handler = self._handlers.get(command_type)
        if handler:
            try:
                result = handler(payload, session_id)
                self._send_response({
                    'command_id': command_id,
                    'success': True,
                    'data': result,
                })
            except Exception as e:
                self._send_response({
                    'command_id': command_id,
                    'success': False,
                    'error': str(e),
                })
        else:
            self._send_response({
                'command_id': command_id,
                'success': False,
                'error': f'Unknown command: {command_type}',
            })

    def _send_response(self, response: Dict[str, Any]):
        """Send a response to the frontend."""
        print(json.dumps(response))
        sys.stdout.flush()

    # Command handlers

    def _handle_open_local_file(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle opening a local file."""
        file_path = payload.get('file_path')
        if not file_path:
            raise ValueError("file_path is required")

        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        # Create or get session
        if not session_id or session_id not in self._sessions:
            session = self.workspace_manager.create_session()
            self._sessions[session.id] = session
            self._current_session_id = session.id
        else:
            session = self._sessions[session_id]
            self._current_session_id = session_id

        # Copy file to workspace
        target_name = path.name
        target = self.workspace_manager.copy_structure_to_workspace(
            session, path, target_name
        )

        # Parse the structure
        content = path.read_text(encoding='utf-8')

        if path.suffix.lower() in ('.cif', '.mcif', '.bcif'):
            structure = self.parser.parse_mmcif(
                content, source_id=path.stem, source='local'
            )
        elif path.suffix.upper() == '.PDB':
            structure = self.parser.parse_pdb(
                content, source_id=path.stem, source='local'
            )
        else:
            # Try mmCIF first
            try:
                structure = self.parser.parse_mmcif(
                    content, source_id=path.stem, source='local'
                )
            except Exception:
                structure = self.parser.parse_pdb(
                    content, source_id=path.stem, source='local'
                )

        # Save structure to workspace
        session.structure = structure
        session.file_path = str(target)

        # Resolve ligands
        for ligand in structure.ligands:
            self.ligand_resolver.resolve_ligand(ligand, structure.id)

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
        """Handle opening a structure from PDB ID."""
        pdb_id = payload.get('pdb_id')
        if not pdb_id:
            raise ValueError("pdb_id is required")

        # Create or get session
        if not session_id or session_id not in self._sessions:
            session = self.workspace_manager.create_session()
            self._sessions[session.id] = session
            self._current_session_id = session.id
        else:
            session = self._sessions[session_id]
            self._current_session_id = session_id

        # Fetch from RCSB
        structure = self.parser.fetch_from_rcsb(pdb_id, format='mmCIF')
        structure.id = pdb_id
        session.structure = structure

        # Save to workspace
        content = self.parser.parse_mmcif(
            Path(f"/dev/null").read_text(),  # dummy
            source_id=pdb_id,
            source='local'
        )
        # Use the parser to get content string
        try:
            import requests
            response = requests.get(
                f"https://files.rcsb.org/download/{pdb_id}.cif",
                timeout=30,
            )
            response.raise_for_status()
            content = response.text
            workspace_path = session.get_workspace() / f'{pdb_id}.cif'
            workspace_path.write_text(content)

            # Re-parse from the saved content
            structure = self.parser.parse_mmcif(
                content, source_id=pdb_id, source='rcsb'
            )
        except Exception as e:
            raise RuntimeError(f"Failed to download structure: {e}")

        # Resolve ligands
        for ligand in structure.ligands:
            self.ligand_resolver.resolve_ligand(ligand, pdb_id)

        return {
            'structure': self._structure_to_dict(structure),
            'session_id': session.id,
            'workspace_path': str(session.get_workspace()),
        }

    def _handle_select_ligand(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle selecting a ligand."""
        ligand_id = payload.get('ligand_id')
        if not ligand_id:
            raise ValueError("ligand_id is required")

        session = self._get_current_session(session_id)
        if not session or not session.structure:
            raise ValueError("No structure loaded")

        # Find the ligand
        ligand = None
        for l in session.structure.ligands:
            if l.id == ligand_id or l.name == ligand_id or l.residue_name == ligand_id:
                ligand = l
                break

        if not ligand:
            raise ValueError(f"Ligand not found: {ligand_id}")

        session.selected_ligand_id = ligand_id

        return {
            'ligand': self._ligand_to_dict(ligand),
            'session_id': session.id,
        }

    def _handle_run_contact_analysis(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle running contact analysis."""
        session = self._get_current_session(session_id)
        if not session or not session.structure:
            raise ValueError("No structure loaded")

        if not session.selected_ligand_id:
            raise ValueError("No ligand selected")

        # Find selected ligand
        ligand = None
        for l in session.structure.ligands:
            if l.id == session.selected_ligand_id or l.name == session.selected_ligand_id:
                ligand = l
                break

        if not ligand:
            raise ValueError(f"Selected ligand not found: {session.selected_ligand_id}")

        # Run contact analysis
        contacts = self.contact_analyzer.analyze_contacts(
            session.structure, ligand
        )

        # Compute binding pocket
        pocket = self.contact_analyzer.compute_binding_pocket(
            session.structure, ligand
        )

        # Run enrichment
        enriched_ligand = self.ligand_resolver.resolve_ligand(
            ligand, session.structure.id
        )

        # Get evidence
        evidence = self.enrichment_client.get_all_evidence(
            {
                'name': enriched_ligand.name,
                'residue_name': enriched_ligand.residue_name,
                'pubchem_cid': enriched_ligand.pubchem_cid,
                'chembl_id': enriched_ligand.chembl_id,
                'pdbbind_affinity': enriched_ligand.pdbbind_affinity,
            },
            pdb_id=session.structure.id,
        )

        # Assign contact IDs
        for i, c in enumerate(contacts):
            c.id = i

        # Export contacts CSV
        contacts_path = session.get_contacts_path()
        self.contact_analyzer.export_contacts_csv(contacts, contacts_path)

        return {
            'contacts': [
                self._contact_to_dict(c) for c in contacts
            ],
            'pocket': pocket,
            'ligand': self._ligand_to_dict(enriched_ligand),
            'evidence': [
                {'source': e.source, 'field': e.field, 'value': e.value}
                for e in evidence
            ],
            'contact_count': len(contacts),
        }

    def _handle_run_docking(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle running docking."""
        session = self._get_current_session(session_id)
        if not session or not session.structure:
            raise ValueError("No structure loaded")

        if not session.selected_ligand_id:
            raise ValueError("No ligand selected")

        # Find selected ligand
        ligand = None
        for l in session.structure.ligands:
            if l.id == session.selected_ligand_id:
                ligand = l
                break

        if not ligand:
            raise ValueError("Selected ligand not found")

        # Get Vina adapter
        vina = self.engine_registry.get_adapter(EngineType.VINA)
        if not vina.is_available():
            raise RuntimeError("Vina not available. Please install AutoDock Vina.")

        # Prepare docking
        output_dir = session.get_job_output_dir('vina_' + str(int(time.time())))
        input_files = vina.prepare(session.structure, ligand, output_dir)

        # Run docking
        job = self.job_manager.create_job(
            EngineType.VINA,
            payload,
            session_id,
        )

        def run_docking(job, params):
            # This would be the actual docking call
            # For now, return a simulated result
            return {
                'poses': [],
                'runtime_seconds': 0,
            }

        self.job_manager.start_job(job, run_docking)

        return {
            'job_id': job.id,
            'status': 'running',
            'message': 'Docking started',
        }

    def _handle_run_geometry_cleanup(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle running geometry cleanup."""
        # Placeholder - would use OpenBabel or similar
        return {
            'status': 'not_implemented',
            'message': 'Geometry cleanup requires OpenBabel or RDKit',
        }

    def _handle_set_docking_box(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle setting docking box."""
        # Store box parameters in session
        session = self._get_current_session(session_id)
        if session:
            if not hasattr(session, 'docking_box'):
                session.docking_box = {}
            session.docking_box.update(payload)
        return {'status': 'ok'}

    def _handle_cancel_job(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle cancelling a job."""
        job_id = payload.get('job_id')
        if not job_id:
            raise ValueError("job_id is required")

        cancelled = self.job_manager.cancel_job(job_id)
        return {'cancelled': cancelled}

    def _handle_measure_distance(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle measuring distance between atoms."""
        atom1 = payload.get('atom1')
        atom2 = payload.get('atom2')

        if not atom1 or not atom2:
            raise ValueError("atom1 and atom2 are required")

        # Find atoms in structure
        session = self._get_current_session(session_id)
        if not session or not session.structure:
            raise ValueError("No structure loaded")

        atom1_obj = self._find_atom(session.structure, atom1)
        atom2_obj = self._find_atom(session.structure, atom2)

        if not atom1_obj or not atom2_obj:
            raise ValueError("Atoms not found")

        # Calculate distance
        pos1 = np.array([atom1_obj.x, atom1_obj.y, atom1_obj.z])
        pos2 = np.array([atom2_obj.x, atom2_obj.y, atom2_obj.z])
        distance = float(np.linalg.norm(pos1 - pos2))

        return {
            'distance': distance,
            'atom1': self._atom_to_dict(atom1_obj),
            'atom2': self._atom_to_dict(atom2_obj),
        }

    def _handle_measure_angle(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle measuring angle between three atoms."""
        atom1 = payload.get('atom1')
        atom2 = payload.get('atom2')
        atom3 = payload.get('atom3')

        if not all([atom1, atom2, atom3]):
            raise ValueError("atom1, atom2, and atom3 are required")

        session = self._get_current_session(session_id)
        if not session or not session.structure:
            raise ValueError("No structure loaded")

        a1 = self._find_atom(session.structure, atom1)
        a2 = self._find_atom(session.structure, atom2)
        a3 = self._find_atom(session.structure, atom3)

        if not all([a1, a2, a3]):
            raise ValueError("Atoms not found")

        # Calculate angle
        v1 = np.array([a1.x, a1.y, a1.z]) - np.array([a2.x, a2.y, a2.z])
        v2 = np.array([a3.x, a3.y, a3.z]) - np.array([a2.x, a2.y, a2.z])

        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)

        angle = float(np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)))
        angle_deg = float(np.degrees(angle))

        return {
            'angle': angle_deg,
            'atom1': self._atom_to_dict(a1),
            'atom2': self._atom_to_dict(a2),
            'atom3': self._atom_to_dict(a3),
        }

    def _handle_save_artifact(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle saving an analysis artifact."""
        session = self._get_current_session(session_id)
        if not session:
            raise ValueError("No session")

        # Create analysis summary
        summary = self.exporter.export_analysis_summary(
            structure=session.structure or Structure(id='', source=''),
            ligand=session.structure and next(
                (l for l in session.structure.ligands
                 if l.id == session.selected_ligand_id), None
            ) or Ligand(id='', name=''),
            contacts=[],
            evidence=[],
            job_results=[],
            notes=session.notes,
        )

        # Save to workspace
        path = self.exporter.save_analysis_summary(session_id, summary)

        return {
            'path': str(path),
            'session_id': session.id,
        }

    def _handle_load_artifact(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle loading an analysis artifact."""
        path = payload.get('path')
        if not path:
            raise ValueError("path is required")

        artifact_path = Path(path)
        if not artifact_path.exists():
            raise FileNotFoundError(f"Artifact not found: {path}")

        # Load summary
        summary = self.exporter.load_analysis_summary(artifact_path)

        # Create new session
        session = self.workspace_manager.create_session()
        self._sessions[session.id] = session
        self._current_session_id = session.id

        return {
            'session_id': session.id,
            'summary': {
                'structure_id': summary.structure_id,
                'structure_title': summary.structure_title,
                'ligand_id': summary.ligand_id,
                'ligand_name': summary.ligand_name,
                'contact_count': summary.contact_count,
                'notes': summary.notes,
            },
        }

    def _handle_get_status(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle getting status."""
        session = self._get_current_session(session_id)

        job_status = self.job_manager.get_status_summary()
        engine_status = self.engine_registry.get_engine_status()

        return {
            'session': {
                'id': session.id if session else None,
                'structure_id': session.structure.id if session and session.structure else None,
                'selected_ligand_id': session.selected_ligand_id,
                'notes': session.notes,
            } if session else None,
            'jobs': job_status,
            'engines': engine_status,
            'data_sources': self.enrichment_client.health_check(),
        }

    def _handle_export_scene_image(
        self,
        payload: Dict[str, Any],
        session_id: str,
    ) -> Dict[str, Any]:
        """Handle exporting scene image."""
        # Placeholder - would use renderer
        return {
            'status': 'not_implemented',
            'message': 'Scene export requires 3D renderer integration',
        }

    # Utility methods

    def _get_current_session(self, session_id: str) -> Optional[Session]:
        """Get the current session."""
        if session_id and session_id in self._sessions:
            return self._sessions[session_id]
        if self._current_session_id and self._current_session_id in self._sessions:
            return self._sessions[self._current_session_id]
        return None

    def _structure_to_dict(self, structure: Structure) -> Dict[str, Any]:
        """Convert a structure to a dictionary."""
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
            'ligands': [
                self._ligand_to_dict(l) for l in structure.ligands
            ],
        }

    def _ligand_to_dict(self, ligand: Ligand) -> Dict[str, Any]:
        """Convert a ligand to a dictionary."""
        return {
            'id': ligand.id,
            'name': ligand.name,
            'residue_name': ligand.residue_name,
            'formula': ligand.formula,
            'molecular_weight': ligand.molecular_weight,
            'atom_count': ligand.atom_count,
            'smiles': ligand.smiles,
            'inchi_key': ligand.inchi_key,
            'pubchem_cid': ligand.pubchem_cid,
            'pubchem_name': ligand.pubchem_name,
            'chembl_id': ligand.chembl_id,
            'pdbbind_affinity': ligand.pdbbind_affinity,
            'classification_hint': ligand.classification_hint,
            'resolution_status': ligand.resolution_status.value if hasattr(ligand.resolution_status, 'value') else str(ligand.resolution_status),
            'has_2d_structure': ligand.has_2d_structure,
            'atoms': [
                self._atom_to_dict(a) for a in ligand.atoms[:500]  # Limit for performance
            ],
        }

    def _contact_to_dict(self, contact: Contact) -> Dict[str, Any]:
        """Convert a contact to a dictionary."""
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
            'contact_type': contact.contact_type.value if hasattr(contact.contact_type, 'value') else str(contact.contact_type),
            'angle': contact.angle,
            'is_water_mediated': contact.is_water_mediated,
            'description': contact.description,
        }

    def _atom_to_dict(self, atom: Any) -> Dict[str, Any]:
        """Convert an atom to a dictionary."""
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

    def _find_atom(self, structure: Structure, atom_id: str) -> Optional[Any]:
        """Find an atom by identifier."""
        # Try to parse as "chain:resid:resname:atomname"
        parts = atom_id.split(':')
        if len(parts) >= 4:
            chain_id, res_id, res_name, atom_name = parts[0], parts[1], parts[2], parts[3]
            try:
                res_id = int(res_id)
            except ValueError:
                pass

            for chain in structure.chains:
                if chain.id != chain_id:
                    continue
                for residue in chain.residues:
                    if residue.id != res_id and str(residue.id) != str(res_id):
                        continue
                    if residue.name != res_name and residue.residue_name != res_name:
                        continue
                    for atom in residue.atoms:
                        if atom.name == atom_name:
                            return atom

        # Try to parse as "chain:resid:name" (PDB format)
        if len(parts) == 3:
            chain_id, res_id, atom_name = parts
            try:
                res_id = int(res_id)
            except ValueError:
                pass

            for chain in structure.chains:
                if chain.id != chain_id:
                    continue
                for residue in chain.residues:
                    if residue.id != res_id and str(residue.id) != str(res_id):
                        continue
                    for atom in residue.atoms:
                        if atom.name == atom_name:
                            return atom

        # Search all atoms
        for chain in structure.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    if (atom.name == atom_id or
                        f"{chain.id}:{atom.residue_id}:{atom.residue_name}:{atom.name}" == atom_id):
                        return atom

        return None


def start_server(config: Optional[Config] = None):
    """Start the backend server."""
    server = BackendServer(config)
    server.start()
    return server


if __name__ == '__main__':
    # Run the server
    server = BackendServer()
    print("Ligora backend server started", file=sys.stderr)
    sys.stderr.flush()
    server.start()

    # Keep main thread alive
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        server.stop()
