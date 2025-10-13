"""Runtime configuration holder.

This module provides a simple, importable place to store parsed
command-line arguments and a chosen QC program. It is intentionally
minimal so existing modules can either import `config` and read
`config.args`/`config.QC_program` or (temporarily) DATS can
monkey-patch other modules at startup to preserve backward
compatibility.
"""
from typing import Any, Optional


class Config:
    def __init__(self) -> None:
        # Will be set to the argparse.Namespace returned by parser.parse_args()
        self.args: Optional[Any] = None
        # String like 'ORCA' or 'G16'
        self.QC_program: Optional[str] = None


# Module-level instance used by other modules: `from src.config import config`
config = Config()
