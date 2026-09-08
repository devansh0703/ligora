"""
Batch analysis for processing multiple structures.
"""

from __future__ import annotations

import uuid
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, List, Dict, Any, Callable
from dataclasses import dataclass, field

from ..schemas import (
    Structure,
    Ligand,
    AnalysisSummary,
)
from ..workspace import WorkspaceManager
from ..parser import StructureParser
from ..ligand import LigandResolver
from ..contacts import ContactAnalyzer
from ..enrichment import EnrichmentClient
from ..config import get_config


class BatchAnalyzer:
    """
    Batch analysis processor for multiple structures.

    Features:
    - Process multiple PDB IDs or files in sequence
    - Progress reporting
    - Parallel processing (optional)
    - Batch result export
    """

    def __init__(
        self,
        workspace_manager: Optional[WorkspaceManager] = None,
        parser: Optional[StructureParser] = None,
        ligand_resolver: Optional[LigandResolver] = None,
        contact_analyzer: Optional[ContactAnalyzer] = None,
        enrichment_client: Optional[EnrichmentClient] = None,
    ):
        """Initialize the batch analyzer."""
        self.config = get_config()
        self.workspace_manager = workspace_manager or WorkspaceManager(
            self.config.workspace_dir
        )
        self.parser = parser or StructureParser()
        self.ligand_resolver = ligand_resolver or LigandResolver()
        self.contact_analyzer = contact_analyzer or ContactAnalyzer()
        self.enrichment_client = enrichment_client or EnrichmentClient()

        self._active_batch_id: Optional[str] = None
        self._batch_jobs: Dict[str, BatchJob] = {}
        self._callbacks: List[Callable] = []

    def create_batch(self, batch_id: Optional[str] = None) -> str:
        """Create a new batch analysis."""
        batch_id = batch_id or str(uuid.uuid4())
        self._batch_jobs[batch_id] = BatchJob(id=batch_id)
        self._active_batch_id = batch_id
        return batch_id

    def add_structure_to_batch(
        self,
        batch_id: str,
        structure_source: str,
        source_type: str = "local",
    ) -> Optional[BatchTask]:
        """Add a structure to the batch."""
        batch = self._batch_jobs.get(batch_id)
        if not batch:
            return None

        task = BatchTask(
            id=str(uuid.uuid4()),
            batch_id=batch_id,
            source=structure_source,
            source_type=source_type,
            status=BatchTaskStatus.PENDING,
        )
        batch.tasks.append(task)
        return task

    def run_batch(
        self,
        batch_id: str,
        on_progress: Optional[Callable] = None,
    ) -> List[BatchResult]:
        """Run the batch analysis."""
        batch = self._batch_jobs.get(batch_id)
        if not batch:
            return []

        results = []

        # Update status
        batch.status = BatchStatus.RUNNING
        batch.started_at = time.time()

        for task in batch.tasks:
            task.status = BatchTaskStatus.RUNNING

            try:
                # Load structure
                structure = self._load_structure(task)

                if not structure:
                    task.status = BatchTaskStatus.FAILED
                    task.error = "Failed to load structure"
                    results.append(BatchResult(
                        task_id=task.id,
                        structure_id=task.source,
                        summary=None,
                        success=False,
                        error=task.error,
                    ))
                    continue

                # Resolve every non-polymer instance against the live CCD
                # first (same as the app loader) so classification is real
                # data, then pick the primary ligand from classified ones.
                for ligand_instance in structure.ligands:
                    self.ligand_resolver.resolve_ligand(
                        ligand_instance,
                        structure.id if structure.id else None,
                    )

                # Get primary ligand
                ligand = self._get_primary_ligand(structure)

                if not ligand:
                    task.status = BatchTaskStatus.FAILED
                    task.error = "No ligand found"
                    results.append(BatchResult(
                        task_id=task.id,
                        structure_id=structure.id or task.source,
                        summary=None,
                        success=False,
                        error=task.error,
                    ))
                    continue

                # Run contact analysis (returns (contacts, plip_available))
                contacts, plip_available = \
                    self.contact_analyzer.analyze_contacts(
                        structure, ligand
                    )
                if not plip_available:
                    raise RuntimeError(
                        "PLIP is not installed; contact analysis cannot "
                        "run - results are never simulated")

                # Resolve ligand identity (the pre-selection pass already
                # resolved every instance; this is the resolver's cached
                # identity path for the chosen primary ligand).
                enriched_ligand = self.ligand_resolver.resolve_ligand(
                    ligand, structure.id if structure.id else None
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
                    pdb_id=structure.id if structure.id else None,
                )

                # Create analysis summary
                summary = AnalysisSummary(
                    structure_id=structure.id,
                    structure_title=structure.title,
                    source=structure.source,
                    ligand_id=ligand.id,
                    ligand_name=ligand.name,
                    ligand_formula=ligand.formula,
                    ligand_smiles=ligand.smiles,
                    resolution_status=enriched_ligand.resolution_status.value,
                    contact_count=len(contacts),
                    contacts=[
                        {
                            'id': c.id,
                            'contact_type': c.contact_type.value,
                            'distance': c.distance,
                            'protein_residue': c.protein_residue_name,
                            'residue_id': c.protein_residue_id,
                            'protein_chain': c.protein_chain_id,
                            'ligand_residue': c.ligand_residue_name,
                        }
                        for c in contacts
                    ],
                    evidence=[
                        {
                            'source': e.source,
                            'field': e.field,
                            'value': str(e.value),
                        }
                        for e in evidence
                    ],
                    job_results=[],
                    notes="",
                    created_at=datetime.now(timezone.utc).isoformat(),
                    exported_at=datetime.now(timezone.utc).isoformat(),
                )

                task.status = BatchTaskStatus.COMPLETED
                task.result = summary
                results.append(BatchResult(
                    task_id=task.id,
                    structure_id=structure.id,
                    summary=summary,
                    success=True,
                ))

            except Exception as e:
                task.status = BatchTaskStatus.FAILED
                task.error = str(e)
                results.append(BatchResult(
                    task_id=task.id,
                    structure_id=task.source,
                    summary=None,
                    success=False,
                    error=str(e),
                ))

            # Report progress
            if on_progress:
                progress = (
                    batch.tasks.index(task) + 1
                ) / len(batch.tasks)
                on_progress(progress, task)

        batch.status = BatchStatus.COMPLETED
        batch.completed_at = time.time()

        return results

    def _load_structure(self, task: 'BatchTask') -> Optional[Structure]:
        """Load a structure from task source."""
        if task.source_type == "local":
            path = Path(task.source)
            if not path.exists():
                return None

            content = path.read_text(encoding='utf-8')
            if path.suffix.lower() in ('.cif', '.mcif'):
                return self.parser.parse_mmcif(
                    content, source_id=path.stem, source='local'
                )
            else:
                return self.parser.parse_pdb(
                    content, source_id=path.stem, source='local'
                )

        elif task.source_type == "pdb_id":
            try:
                return self.parser.fetch_from_rcsb(task.source)
            except Exception:
                return None

        return None

    def _get_primary_ligand(self, structure: Structure) -> Optional[Ligand]:
        """Pick the most relevant non-polymer ligand for batch analysis.

        Classification comes from the live CCD via the resolver (pdbx_type
        codes): solvent (HETAS) and ions (HETAI) are never analysis targets.
        """
        for ligand in structure.ligands:
            if ligand.classification_hint in ('HETAS', 'HETAI'):
                continue
            if ligand.atoms:
                return ligand
        return None

    def export_batch_results(
        self,
        batch_id: str,
        output_dir: Optional[Path] = None,
    ) -> Path:
        """Export batch analysis results."""
        batch = self._batch_jobs.get(batch_id)
        if not batch:
            raise ValueError(f"Batch not found: {batch_id}")

        output_dir = output_dir or Path(self.config.base_dir) / f"batch_{batch_id}"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save summary JSON
        summary_path = output_dir / "batch_summary.json"
        summary_data = {
            'batch_id': batch.id,
            'total_tasks': len(batch.tasks),
            'completed': sum(
                1 for t in batch.tasks if t.status == BatchTaskStatus.COMPLETED),
            'failed': sum(
                1 for t in batch.tasks if t.status == BatchTaskStatus.FAILED),
            'results': [
                {
                    'task_id': r.task_id,
                    'structure_id': r.structure_id,
                    'success': r.success,
                    'error': r.error,
                }
                for r in self._get_results_for_batch(batch_id)
            ]
        }

        with open(summary_path, 'w') as f:
            json.dump(summary_data, f, indent=2)

        # Export each result
        for result in self._get_results_for_batch(batch_id):
            if result.success and result.summary:
                result_dir = output_dir / f"result_{result.task_id}"
                result_dir.mkdir(exist_ok=True)

                # Save summary
                with open(result_dir / "analysis_summary.json", 'w') as f:
                    json.dump({
                        'structure_id': result.summary.structure_id,
                        'structure_title': result.summary.structure_title,
                        'ligand_name': result.summary.ligand_name,
                        'contact_count': result.summary.contact_count,
                        'contacts': result.summary.contacts,
                        'evidence': result.summary.evidence,
                    }, f, indent=2)

                # Export contacts CSV
                csf_path = result_dir / "contacts.csv"
                with open(csf_path, 'w') as f:
                    f.write("id,contact_type,distance,protein_residue,ligand_residue\n")
                    for c in result.summary.contacts:
                        f.write(f"{c['id']},{c['contact_type']},{c['distance']},{c['protein_residue']},{c['ligand_residue']}\n")

        return output_dir

    def _get_results_for_batch(self, batch_id: str) -> List['BatchResult']:
        """Get results for a batch (one BatchResult per finished task)."""
        batch = self._batch_jobs.get(batch_id)
        if not batch:
            return []

        out: List[BatchResult] = []
        for task in batch.tasks:
            if isinstance(task.result, BatchResult):
                out.append(task.result)
            elif task.result is not None:
                # Success path stores the AnalysisSummary on the task.
                out.append(BatchResult(
                    task_id=task.id,
                    structure_id=(
                        task.result.structure_id
                        if getattr(task.result, 'structure_id', None)
                        else task.source),
                    summary=task.result,
                    success=task.status == BatchTaskStatus.COMPLETED,
                    error=task.error,
                ))
            elif task.status == BatchTaskStatus.FAILED:
                out.append(BatchResult(
                    task_id=task.id,
                    structure_id=task.source,
                    summary=None,
                    success=False,
                    error=task.error,
                ))
        return out


@dataclass
class BatchJob:
    """A batch analysis job."""
    id: str
    status: 'BatchStatus' = field(default_factory=lambda: BatchStatus.PENDING)
    tasks: List['BatchTask'] = field(default_factory=list)
    started_at: Optional[float] = None
    completed_at: Optional[float] = None


@dataclass
class BatchTask:
    """A single task in a batch."""
    id: str
    batch_id: str
    source: str
    source_type: str
    status: 'BatchTaskStatus' = field(default_factory=lambda: BatchTaskStatus.PENDING)
    result: Optional[Any] = None
    error: Optional[str] = None


@dataclass
class BatchResult:
    """Result of a batch task."""
    task_id: str
    structure_id: str
    summary: Optional[AnalysisSummary]
    success: bool
    error: Optional[str] = None


class BatchStatus:
    """Batch job status."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class BatchTaskStatus:
    """Batch task status."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"
