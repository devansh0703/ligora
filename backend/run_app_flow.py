#!/usr/bin/env python3
"""Local repro runner for the Ligora backend app flow.

Usage:
  PYTHONPATH=backend python3 run_app_flow.py
"""
from pathlib import Path
import sys
import json

sys.path.insert(0, str(Path(__file__).resolve().parent / "backend"))

from ligora_backend.server import BackendServer
from ligora_backend.workspace import WorkspaceManager

session = WorkspaceManager(BackendServer().config.workspace_dir).create_session()
server = BackendServer()
open_result = server._handle_open_pdb_id({"pdb_id": "1hcl"}, session.id)
session_id = open_result["session_id"]
structure = server._sessions[session_id].structure
lig = structure.ligands[0]
server._sessions[session_id].selected_ligand_id = lig.id

select_result = server._handle_select_ligand({"ligand_id": lig.id}, session_id)
analysis_result = server._handle_run_contact_analysis({}, session_id)
status_result = server._handle_get_status({}, session_id)

print(json.dumps({
    "structure": {
        "id": open_result["structure"]["id"],
        "title": open_result["structure"]["title"],
        "chains": [c["id"] for c in open_result["structure"]["chains"]],
        "ligands": [lig["residue_name"] for lig in open_result["structure"]["ligands"]],
    },
    "select": {
        "ligand_name": select_result["ligand"]["name"],
        "ligand_residue_name": select_result["ligand"]["residue_name"],
        "classification_hint": select_result["ligand"]["classification_hint"],
        "resolution_status": select_result["ligand"]["resolution_status"],
        "pubchem_cid": select_result["ligand"]["pubchem_cid"],
        "iupac_name": select_result["ligand"].get("iupac_name"),
    },
    "analysis": {
        "contact_count": analysis_result["contact_count"],
        "pocket_residue_count": analysis_result["pocket"]["pocket_residue_count"],
        "ligand_name": analysis_result["ligand"]["name"],
        "ligand_smiles": analysis_result["ligand"]["smiles"],
        "ligand_formula": analysis_result["ligand"]["formula"],
        "evidence_count": len(analysis_result["evidence"]),
        "ligand_resolved": analysis_result.get("ligand_resolved"),
    },
    "status": {
        "engines": {k: v.get("available") for k, v in status_result["engines"].items()},
        "data_sources": status_result["data_sources"],
    },
}, default=str))
