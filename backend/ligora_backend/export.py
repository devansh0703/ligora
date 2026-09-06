"""
Artifact exporter for saving and loading analysis results.

Provides:
- Export of analysis summaries as JSON
- Export of contacts as CSV
- Export of ligands as SDF/MOL
- Export of full analysis artifacts
"""

import json
import csv
import shutil
from pathlib import Path
from datetime import datetime, timezone
from typing import Optional, Dict, Any, List

import numpy as np

from .schemas import (
    Structure,
    Ligand,
    Contact,
    EvidenceItem,
    AnalysisSummary,
    DockingResult,
    GeometryCleanupResult,
    Atom,
)


class ArtifactExporter:
    """Exports analysis artifacts in various formats."""

    def __init__(self, workspace_manager=None):
        self._workspace_manager = workspace_manager

    def export_analysis_summary(
        self,
        structure: Structure,
        ligand: Ligand,
        contacts: List[Contact],
        evidence: List[Dict[str, Any]],
        job_results: List[Dict[str, Any]],
        notes: str,
        **kwargs
    ) -> AnalysisSummary:
        """Create an analysis summary object."""
        return AnalysisSummary(
            structure_id=structure.id,
            structure_title=structure.title,
            source=structure.source,
            ligand_id=ligand.id,
            ligand_name=ligand.name,
            ligand_formula=ligand.formula,
            ligand_smiles=ligand.smiles,
            resolution_status=str(ligand.resolution_status),
            contact_count=len(contacts),
            contacts=[{
                'id': c.id,
                'ligand_atom': c.ligand_atom,
                'ligand_residue': c.ligand_residue_name,
                'ligand_residue_id': c.ligand_residue_id,
                'ligand_chain': c.ligand_chain_id,
                'protein_atom': c.protein_atom,
                'protein_residue': c.protein_residue_name,
                'protein_residue_id': c.protein_residue_id,
                'protein_chain': c.protein_chain_id,
                'distance': c.distance,
                'contact_type': str(c.contact_type),
                'angle': c.angle,
                'water_mediated': c.is_water_mediated,
                'description': c.description,
            } for c in contacts],
            evidence=evidence,
            job_results=job_results,
            notes=notes,
            created_at=datetime.now(timezone.utc).isoformat(),
            exported_at=datetime.now(timezone.utc).isoformat(),
        )

    def save_analysis_summary(
        self,
        session_id: str,
        summary: AnalysisSummary,
        output_path: Optional[Path] = None,
    ) -> Path:
        """Save analysis summary to JSON file."""
        if output_path is None:
            output_path = Path.cwd() / 'analysis_summary.json'

        data = {
            'structure_id': summary.structure_id,
            'structure_title': summary.structure_title,
            'source': summary.source,
            'ligand_id': summary.ligand_id,
            'ligand_name': summary.ligand_name,
            'ligand_formula': summary.ligand_formula,
            'ligand_smiles': summary.ligand_smiles,
            'resolution_status': summary.resolution_status,
            'contact_count': summary.contact_count,
            'contacts': summary.contacts,
            'evidence': summary.evidence,
            'job_results': summary.job_results,
            'notes': summary.notes,
            'created_at': summary.created_at,
            'exported_at': summary.exported_at,
        }

        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)

        return output_path

    def export_contacts_csv(
        self,
        contacts: List[Contact],
        output_path: Path,
    ) -> Path:
        """Export contacts to CSV."""
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'ID', 'Ligand_Atom', 'Ligand_Residue', 'Ligand_Residue_ID',
                'Ligand_Chain', 'Protein_Atom', 'Protein_Residue',
                'Protein_Residue_ID', 'Protein_Chain', 'Distance',
                'Contact_Type', 'Angle', 'Water_Mediated', 'Description'
            ])
            for c in contacts:
                writer.writerow([
                    c.id, c.ligand_atom, c.ligand_residue_name,
                    c.ligand_residue_id, c.ligand_chain_id,
                    c.protein_atom, c.protein_residue_name,
                    c.protein_residue_id, c.protein_chain_id,
                    round(c.distance, 3), str(c.contact_type),
                    c.angle if c.angle else '', c.is_water_mediated, c.description
                ])
        return output_path

    def export_ligand_sdf(
        self,
        ligand: Ligand,
        output_path: Path,
    ) -> Path:
        """Export ligand to SDF format."""
        lines = []
        mol_name = ligand.name or ligand.residue_name or 'LIG'
        lines.append(mol_name[:80].ljust(80))
        lines.append(''.ljust(80))
        n_atoms = len(ligand.atoms)
        n_bonds = max(0, n_atoms - 1)
        lines.append(f'{n_atoms:>3} {n_bonds:>3}'.ljust(16))
        lines.append(''.ljust(80))
        lines.append(''.ljust(80))
        lines.append(''.ljust(80))

        for i, atom in enumerate(ligand.atoms):
            elem = (atom.element or 'C').strip()
            atomic_num = self._atomic_number(elem)
            x, y, z = atom.x, atom.y, atom.z
            lines.append(
                f'{atomic_num:>3} {elem:<2} {x:>10.4f} {y:>10.4f} {z:>10.4f}'
                f'    0.0000           0'
            )

        for i in range(1, n_atoms):
            lines.append(f'{i:>3} {i+1:>3} {1:>3}')

        lines.append('$$$$')

        with open(output_path, 'w') as f:
            f.write('\n'.join(lines))

        return output_path

    def _atomic_number(self, symbol: str) -> int:
        """Get atomic number from element symbol."""
        table = {
            'H': 1, 'HE': 2, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'NE': 10,
            'NA': 11, 'MG': 12, 'AL': 13, 'SI': 14, 'P': 15, 'S': 16,
            'CL': 17, 'AR': 18, 'K': 19, 'CA': 20, 'BR': 35, 'I': 53,
        }
        return table.get(symbol.upper(), 6)

    def export_full_artifact(
        self,
        session_id: str,
        summary: AnalysisSummary,
        output_dir: Optional[Path] = None,
    ) -> Path:
        """Export complete analysis artifact package."""
        if output_dir is None:
            output_dir = Path.cwd() / f'artifact_{summary.structure_id}'

        output_dir.mkdir(parents=True, exist_ok=True)

        self.save_analysis_summary(session_id, summary, output_dir / 'analysis_summary.json')
        contacts_path = output_dir / 'contacts.csv'
        self.export_contacts_csv(summary.contacts, contacts_path)

        return output_dir

    def format_summary_text(self, summary: AnalysisSummary) -> str:
        """Format summary as human-readable text report."""
        lines = ['=' * 80]
        lines.append(f'ANALYSIS REPORT: {summary.structure_title}')
        lines.append('=' * 80)
        lines.append('')
        lines.append(f'Structure: {summary.structure_id}')
        lines.append(f'Ligand: {summary.ligand_name} ({summary.ligand_formula or "N/A"})')
        lines.append(f'Contacts: {summary.contact_count}')
        lines.append(f'Status: {summary.resolution_status}')
        if summary.notes:
            lines.append('')
            lines.append(f'Notes: {summary.notes}')
        lines.append('')
        lines.append('=' * 80)
        return '\n'.join(lines)
