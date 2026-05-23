"""
Compatibility entry point for the modern FRAP-Toolbox app.

The current GUI direction is a local Streamlit browser app because it is easier
to install and run across lab workstations than a native desktop bundle while
the Python port is still evolving.
"""

from __future__ import annotations

from ..app_launcher import main

__all__ = ["main"]
