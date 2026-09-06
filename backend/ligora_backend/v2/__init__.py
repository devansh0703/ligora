"""
V2 features for Ligora.

Includes:
- Batch analysis
- 2D ligand editor sync
- Water network analysis
- Scripting console
- Result diff/compare
"""

from .batch import BatchAnalyzer
from .editors import Ligand2DEditor
from .water import WaterNetworkAnalyzer
from .scripting import ScriptingConsole
from .comparison import ResultComparator

__all__ = [
    "BatchAnalyzer",
    "Ligand2DEditor",
    "WaterNetworkAnalyzer",
    "ScriptingConsole",
    "ResultComparator",
]
