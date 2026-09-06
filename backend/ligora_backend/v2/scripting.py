"""
Scripting console for Ligora.

Allows users to run Python scripts for custom analysis and automation.
"""

from typing import Optional, Dict, Any, List
from pathlib import Path

from ..config import get_config


class ScriptingConsole:
    """
    Scripting console for Ligora.

    Provides:
    - Python script execution
    - Access to analysis data
    - Custom analysis workflows
    - Automation scripts
    """

    def __init__(self):
        """Initialize the console."""
        self.config = get_config()
        self._variables: Dict[str, Any] = {}
        self._history: List[str] = []
        self._imports: Dict[str, Any] = {}

    def execute(self, code: str,
                context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Execute a Python script.

        Args:
            code: Python code to execute.
            context: Session variables (e.g. structure, selected_ligand_id)
                exposed to the script.

        Returns:
            Execution result.
        """
        try:
            # Prepare namespace
            namespace = {
                **self._imports,
                **self._variables,
                **(context or {}),
                'ligora': self._get_ligora_api(),
            }

            # Execute
            exec(compile(code, '<console>', 'exec'), namespace)

            # Update variables
            self._variables.update({
                k: v for k, v in namespace.items()
                if not k.startswith('_') and not callable(v)
                and k not in (context or {})
            })

            # Record history
            self._history.append(code)

            return {
                'success': True,
                'result': None,
                'variables': list(self._variables.keys()),
            }

        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'variables': list(self._variables.keys()),
            }

    def execute_file(self, path: str) -> Dict[str, Any]:
        """Execute a script from a file."""
        file_path = Path(path)

        if not file_path.exists():
            return {
                'success': False,
                'error': f"File not found: {path}",
            }

        code = file_path.read_text(encoding='utf-8')
        return self.execute(code)

    def get_history(self) -> List[str]:
        """Get command history."""
        return self._history.copy()

    def clear_history(self):
        """Clear command history."""
        self._history.clear()

    def get_variables(self) -> Dict[str, Any]:
        """Get current variables."""
        return self._variables.copy()

    def clear_variables(self):
        """Clear all variables."""
        self._variables.clear()

    def _get_ligora_api(self) -> Dict[str, Any]:
        """Get the Ligora API for scripting."""
        return {
            'config': self.config.to_dict(),
            'version': '0.1.0',
        }
