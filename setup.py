"""Builds the `trash_py._ext` C extension. Project metadata lives in
`pyproject.toml`; this file exists only because pyproject.toml does not
yet support `ext_modules` declaratively across the build backends we
care about.
"""
from __future__ import annotations

import sys

from setuptools import Extension, setup


def _compile_args() -> list[str]:
    if sys.platform == "win32":
        return ["/O2"]
    return ["-O3"]


setup(
    ext_modules=[
        Extension(
            "trash_py._ext",
            sources=["src/trash_py/_ext.c"],
            extra_compile_args=_compile_args(),
        )
    ]
)
