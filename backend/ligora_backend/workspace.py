"""
Workspace and session management for Ligora.

Each analysis session gets its own workspace directory with:
- structure files
- ligand files
- engine input/output files
- analysis artifacts (CSV, SDF, images, JSON summaries)
- job logs
"""
import uuid
import json
import shutil
from pathlib import Path
from datetime import datetime, timezone
from typing import Optional, Dict, Any, List
from dataclasses import dataclass, field

from .schemas import (
    Structure,
    Job,
)


@dataclass
class Session:
    """
    Analysis session state.

    A session owns:
    - The loaded structure
    - Selected ligand
    - Jobs (running/completed)
    - Notes and annotations
    - Workspace path for all files
    """
    id: str
    structure: Optional[Structure] = None
    selected_ligand_id: Optional[str] = None
    jobs: Dict[str, Job] = field(default_factory=dict)
    active_job_id: Optional[str] = None
    notes: str = ""
    docking_box: Optional[Dict[str, Any]] = None
    created_at: str = ""
    workspace_path: str = ""

    def __post_init__(self):
        if not self.id:
            self.id = str(uuid.uuid4())
        if not self.created_at:
            self.created_at = datetime.now(timezone.utc).isoformat()

    def get_workspace(self) -> Path:
        """Get the workspace directory path."""
        return Path(self.workspace_path)

    def get_structure_path(self) -> Path:
        """Get path to the structure file."""
        return self.get_workspace() / "structure.mmcif"

    def get_ligand_path(self, ligand_id: str) -> Path:
        """Get path to a ligand SDF file."""
        return self.get_workspace() / f"ligand_{ligand_id}.sdf"

    def get_contacts_path(self) -> Path:
        """Get path to the contacts CSV export."""
        return self.get_workspace() / "contacts.csv"

    def get_summary_path(self) -> Path:
        """Get path to the analysis summary JSON."""
        return self.get_workspace() / "analysis_summary.json"

    def get_scene_image_path(self) -> Path:
        """Get path to the scene image."""
        return self.get_workspace() / "scene.png"

    def get_job_log_path(self, job_id: str) -> Path:
        """Get path to a job's log file."""
        return self.get_workspace() / f"job_{job_id}.log"

    def get_job_output_dir(self, job_id: str) -> Path:
        """Get the output directory for a job."""
        return self.get_workspace() / f"job_{job_id}" / "output"

    def add_job(self, job: Job):
        """Add a job to the session."""
        self.jobs[job.id] = job
        if self.active_job_id is None:
            self.active_job_id = job.id

    def complete_job(self, job: Job):
        """Mark a job as complete in the session."""
        if job.id in self.jobs:
            self.jobs[job.id] = job
        if self.active_job_id == job.id:
            self.active_job_id = None

    def add_note(self, note: str):
        """Append a note to the session."""
        if self.notes:
            self.notes += "\n" + note
        else:
            self.notes = note


class WorkspaceManager:
    """
    Manages workspaces for analysis sessions.

    Creates, manages, and cleans up workspace directories.
    """

    def __init__(self, base_dir: Path):
        """
        Initialize the workspace manager.

        Args:
            base_dir: Base directory for all workspaces.
        """
        self.base_dir = base_dir
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self._sessions: Dict[str, Session] = {}

    def create_session(self) -> Session:
        """
        Create a new analysis session with its own workspace.

        Returns:
            A new Session instance with a dedicated workspace directory.
        """
        session_id = str(uuid.uuid4())
        workspace = self.base_dir / session_id
        workspace.mkdir(parents=True, exist_ok=True)

        session = Session(
            id=session_id,
            workspace_path=str(workspace),
        )
        self._sessions[session_id] = session
        return session

    def get_session(self, session_id: str) -> Optional[Session]:
        """
        Get an existing session by ID.

        Args:
            session_id: The session ID.

        Returns:
            The Session if found, None otherwise.
        """
        return self._sessions.get(session_id)

    def remove_session(self, session_id: str):
        """
        Remove a session and clean up its workspace.

        Args:
            session_id: The session ID to remove.
        """
        if session_id in self._sessions:
            session = self._sessions[session_id]
            workspace = Path(session.workspace_path)
            if workspace.exists():
                shutil.rmtree(workspace)
            del self._sessions[session_id]

    def copy_structure_to_workspace(
        self,
        session: Session,
        source_path: Path,
        target_name: str = "structure.mmcif"
    ) -> Path:
        """
        Copy a structure file into the session's workspace.

        Args:
            session: The session to copy into.
            source_path: Path to the source file.
            target_name: Name for the target file.

        Returns:
            Path to the copied file.
        """
        target = session.get_workspace() / target_name
        shutil.copy2(source_path, target)
        return target

    def save_json(self, session: Session, filename: str, data: Any):
        """
        Save JSON data to the session's workspace.

        Args:
            session: The session.
            filename: Name of the file.
            data: Data to serialize.
        """
        target = session.get_workspace() / filename
        with open(target, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, default=str)

    def load_json(self, session: Session, filename: str) -> Any:
        """
        Load JSON data from the session's workspace.

        Args:
            session: The session.
            filename: Name of the file.

        Returns:
            The deserialized data.
        """
        target = session.get_workspace() / filename
        with open(target, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_text(self, session: Session, filename: str, text: str):
        """
        Save text to the session's workspace.

        Args:
            session: The session.
            filename: Name of the file.
            text: Text content.
        """
        target = session.get_workspace() / filename
        with open(target, "w", encoding="utf-8") as f:
            f.write(text)

    def load_text(self, session: Session, filename: str) -> str:
        """
        Load text from the session's workspace.

        Args:
            session: The session.
            filename: Name of the file.

        Returns:
            The file contents.
        """
        target = session.get_workspace() / filename
        with open(target, "r", encoding="utf-8") as f:
            return f.read()

    def list_files(self, session: Session, pattern: str = "*") -> List[Path]:
        """
        List files in the session's workspace matching a pattern.

        Args:
            session: The session.
            pattern: Glob pattern to match.

        Returns:
            List of matching file paths.
        """
        return list(session.get_workspace().glob(pattern))

    def cleanup_session(self, session: Session):
        """
        Clean up a session's workspace (keep structure, remove temp files).

        Args:
            session: The session to clean up.
        """
        workspace = session.get_workspace()
        # Remove temporary job output directories
        import shutil
        for job_dir in workspace.glob('job_*'):
            if job_dir.is_dir():
                shutil.rmtree(job_dir)
        # Remove temporary files
        for temp_file in workspace.glob('*.tmp'):
            temp_file.unlink()
        # Remove PDB files created for PLIP (if any remain)
        for pdb_file in workspace.glob('*_for_plip.pdb'):
            pdb_file.unlink()
