"""
Reproducible-analysis manifest.

A "reproduce this analysis" record: a JSON document of the real inputs,
parameters, engine outputs, and artifact paths behind everything a session
produced. Every entry cites what actually happened — file hashes, job
parameters and results, contact counts, notes — so a reader (or reviewer)
can re-run the workflow. Nothing is reconstructed from memory or
approximated; a value that was never recorded is simply not in the
manifest.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


def _sha256(path: Path) -> Optional[str]:
    """Hash a real artifact file; absent files are simply not hashed."""
    if not path.exists() or not path.is_file():
        return None
    h = hashlib.sha256()
    try:
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 16), b""):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


def _artifact_entry(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "size_bytes": path.stat().st_size,
    }


def build_repro_manifest(
    session,  # workspace.Session
    last_contacts_count: Optional[int] = None,
    last_contacts_ligand: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Assemble the manifest from real session state.

    `session` is a workspace.Session. Everything included is read from the
    session, its workspace files, and its recorded jobs.
    """
    structure = session.structure
    workspace = session.get_workspace()

    # The structure file: the parser may record it on the structure or
    # the loader on the session; otherwise fall back to the workspace
    # convention. Only files that actually exist become artifact entries.
    structure_file = None
    if structure and structure.file_path:
        structure_file = Path(structure.file_path)
    elif session.file_path:
        structure_file = Path(session.file_path)
    else:
        structure_file = workspace / "structure.mmcif"

    artifacts: Dict[str, Optional[Dict[str, Any]]] = {}
    for name, path in {
        "structure": structure_file,
        "contacts_csv": workspace / "contacts.csv",
        "analysis_summary": workspace / "analysis_summary.json",
        "scene_image": workspace / "scene.png",
        "struct_conn": workspace / "struct_conn.cif",
        "interaction_diagram": workspace / "interaction_diagram.svg",
    }.items():
        entry = _artifact_entry(Path(path))
        if entry:
            artifacts[name] = entry

    jobs: List[Dict[str, Any]] = []
    for job in session.jobs.values():
        jobs.append({
            "job_id": job.id,
            "type": job.type.value,
            "status": job.status.value,
            "parameters": job.parameters,
            "started_at": job.started_at,
            "finished_at": job.finished_at,
            "error": job.error,
            # Results are included verbatim when the engine recorded a
            # JSON-serializable dict (docking poses, energies, ...).
            "result": job.result
            if isinstance(job.result, dict) else None,
        })

    selected_ligand = None
    if session.selected_ligand_id and structure:
        for ligand in structure.ligands:
            if ligand.id == session.selected_ligand_id:
                selected_ligand = {
                    "id": ligand.id,
                    "residue_name": ligand.residue_name,
                    "name": ligand.name,
                    "formula": ligand.formula,
                    "smiles": ligand.smiles,
                }
                break

    manifest: Dict[str, Any] = {
        "manifest_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "session": {
            "id": session.id,
            "created_at": session.created_at,
            "workspace_path": str(workspace),
        },
        "structure": {
            "id": structure.id if structure else None,
            "source": structure.source if structure else None,
            "title": structure.title if structure else None,
            "file_format": structure.file_format if structure else None,
            "resolution": structure.resolution if structure else None,
            "experiment_type": structure.experiment_type if structure
            else None,
        },
        "selected_ligand": selected_ligand,
        "docking_box": session.docking_box or None,
        "contact_analysis": {
            "performed": last_contacts_count is not None,
            "contact_count": last_contacts_count,
            "ligand_id": last_contacts_ligand,
        },
        "jobs": jobs,
        "artifacts": artifacts,
        "notes": session.notes or "",
    }

    # Engine/config parameter snapshot: the actual settings the backend
    # would use for a re-run (all environment-sourced configuration).
    from ..config import get_config
    cfg = get_config()
    manifest["reproduction_parameters"] = {
        "water_network_cutoff": cfg.water_network_cutoff,
        "water_contact_cutoff": cfg.water_contact_cutoff,
        "docking_exhaustiveness": cfg.default_docking_exhaustiveness,
        "docking_num_modes": cfg.default_docking_num_modes,
        "docking_energy_range": cfg.default_docking_energy_range,
        "md_force_field": cfg.md_force_field,
        "md_water_model": cfg.md_water_model,
        "geometry_cleanup_force_field": cfg.geometry_cleanup_force_field,
        "geometry_cleanup_steps": cfg.geometry_cleanup_steps,
    }

    return manifest


def write_repro_manifest(
    session,  # workspace.Session
    last_contacts_count: Optional[int] = None,
    last_contacts_ligand: Optional[str] = None,
    output_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Build the manifest and write it into the session workspace."""
    manifest = build_repro_manifest(
        session, last_contacts_count, last_contacts_ligand)
    path = Path(output_path) if output_path else \
        session.get_workspace() / "repro_manifest.json"
    path.write_text(
        json.dumps(manifest, indent=2, default=str) + "\n",
        encoding="utf-8")
    return {
        "manifest": manifest,
        "path": str(path),
        "sha256": _sha256(path),
    }
