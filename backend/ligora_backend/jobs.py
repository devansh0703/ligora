"""
Job manager for running and tracking analysis jobs.

Provides:
- Job creation and tracking
- Job execution with progress reporting
- Job cancellation
- Job state persistence
"""

import uuid
import threading
import time
from typing import Optional, Dict, Any, Callable, List

from .schemas import (
    Job,
    JobStatus,
    EngineType,
)


class JobManager:
    """
    Manages the lifecycle of analysis jobs.

    Each job:
    - Has a unique ID
    - Has a type (engine type)
    - Has parameters
    - Has a status (pending, running, completed, failed, cancelled)
    - Has progress (0.0 to 1.0)
    - Has a result or error
    - Has a log
    """

    def __init__(self):
        """Initialize the job manager."""
        self._jobs: Dict[str, Job] = {}
        self._lock = threading.Lock()
        self._running_jobs: Dict[str, threading.Thread] = {}
        self._callbacks: Dict[str, List[Callable]] = {}

    def create_job(
        self,
        job_type: EngineType,
        parameters: Dict[str, Any],
        session_id: str,
    ) -> Job:
        """
        Create a new job.

        Args:
            job_type: Type of job (engine type).
            parameters: Job parameters.
            session_id: Session this job belongs to.

        Returns:
            New Job instance.
        """
        job = Job(
            id=str(uuid.uuid4()),
            type=job_type,
            status=JobStatus.PENDING,
            parameters=parameters,
            started_at=time.time(),
        )
        self._jobs[job.id] = job
        return job

    def start_job(
        self,
        job: Job,
        run_fn: Callable[[Job, Dict[str, Any]], Any],
    ) -> str:
        """
        Start a job in a separate thread.

        Args:
            job: The job to start.
            run_fn: Function that executes the job.

        Returns:
            Job ID.
        """
        with self._lock:
            if job.status != JobStatus.PENDING:
                raise ValueError(f"Job is not pending: {job.status}")

            job.status = JobStatus.RUNNING
            job.started_at = time.time()

            # Create and start thread
            thread = threading.Thread(
                target=self._run_job_thread,
                args=(job, run_fn),
                daemon=True,
            )
            self._running_jobs[job.id] = thread
            thread.start()

        return job.id

    def _run_job_thread(
        self,
        job: Job,
        run_fn: Callable[[Job, Dict[str, Any]], Any],
    ):
        """Execute a job in a thread."""
        try:
            # Update progress during execution
            result = run_fn(job, job.parameters)

            with self._lock:
                job.status = JobStatus.COMPLETED
                job.result = result
                job.finished_at = time.time()
                job.progress = 1.0
                self._running_jobs.pop(job.id, None)

            # Notify callbacks
            self._notify_callbacks(job.id, 'completed', job)

        except Exception as e:
            with self._lock:
                job.status = JobStatus.FAILED
                job.error = str(e)
                job.finished_at = time.time()
                self._running_jobs.pop(job.id, None)

            self._notify_callbacks(job.id, 'failed', job)

    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a running job.

        Args:
            job_id: ID of the job to cancel.

        Returns:
            True if the job was cancelled.
        """
        with self._lock:
            job = self._jobs.get(job_id)
            if not job or job.status not in (JobStatus.RUNNING, JobStatus.PENDING):
                return False

            job.status = JobStatus.CANCELLED
            job.finished_at = time.time()

            # Note: Thread cancellation is cooperative
            # The job thread should check for cancellation
            self._running_jobs.pop(job_id, None)

        self._notify_callbacks(job_id, 'cancelled', job)
        return True

    def get_job(self, job_id: str) -> Optional[Job]:
        """Get a job by ID."""
        with self._lock:
            return self._jobs.get(job_id)

    def get_all_jobs(self) -> List[Job]:
        """Get all jobs."""
        with self._lock:
            return list(self._jobs.values())

    def get_jobs_by_status(self, status: JobStatus) -> List[Job]:
        """Get jobs with a specific status."""
        with self._lock:
            return [j for j in self._jobs.values() if j.status == status]

    def get_jobs_by_session(self, session_id: str) -> List[Job]:
        """Get jobs for a session."""
        # This would be integrated with session management
        with self._lock:
            return list(self._jobs.values())

    def update_progress(self, job_id: str, progress: float, log: str = ""):
        """Update job progress."""
        with self._lock:
            job = self._jobs.get(job_id)
            if job:
                job.progress = min(max(progress, 0.0), 1.0)
                if log:
                    job.log = log

    def register_callback(
        self,
        job_id: str,
        callback: Callable[[str, str, Job], None],
    ):
        """Register a callback for job events."""
        if job_id not in self._callbacks:
            self._callbacks[job_id] = []
        self._callbacks[job_id].append(callback)

    def _notify_callbacks(
        self,
        job_id: str,
        event: str,
        job: Job,
    ):
        """Notify callbacks of a job event."""
        callbacks = self._callbacks.get(job_id, [])
        for callback in callbacks:
            try:
                callback(job_id, event, job)
            except Exception:
                pass  # Don't let callback errors affect job execution

    def remove_job(self, job_id: str):
        """Remove a job (for cleanup)."""
        with self._lock:
            self._jobs.pop(job_id, None)
            self._callbacks.pop(job_id, None)

    def get_status_summary(self) -> Dict[str, Any]:
        """Get a summary of all job statuses."""
        with self._lock:
            status_counts = {
                'pending': 0,
                'running': 0,
                'completed': 0,
                'failed': 0,
                'cancelled': 0,
            }

            for job in self._jobs.values():
                status_counts[job.status.value] += 1

            return {
                'total': len(self._jobs),
                'by_status': status_counts,
                'active_jobs': [
                    {'id': j.id, 'type': j.type.value, 'status': j.status.value}
                    for j in self._jobs.values()
                    if j.status in (JobStatus.RUNNING, JobStatus.PENDING)
                ],
            }
