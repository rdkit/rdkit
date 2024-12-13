"""The pytest entry point.

This is pytest.__main__ patched to work when RDKit Python modules
are installed out-of-tree"""

import sys
import re
from pathlib import Path
import _pytest.pathlib
import importlib
import importlib.util

if not hasattr(_pytest.pathlib.import_path, "_before_rdkit_patch"):
  orig_import_path = _pytest.pathlib.import_path

  def import_path_rdkit_patch(p, *, mode, root):
    if not hasattr(import_path_rdkit_patch, "_rdkit_path"):
      setattr(import_path_rdkit_patch, "_rdkit_path", Path(importlib.import_module("rdkit").__file__).parent)
      setattr(import_path_rdkit_patch, "PY_EXT_REGEX", re.compile(r"\.py$"))
    rdkit_install_path = import_path_rdkit_patch._rdkit_path.joinpath(p.relative_to(root))
    if p.name != "conftest.py" and rdkit_install_path.exists():
      p = rdkit_install_path
    else:
      spec = importlib.util.spec_from_file_location(p.name, p.as_posix())
      mod = importlib.util.module_from_spec(spec)
      mod_path = import_path_rdkit_patch.PY_EXT_REGEX.sub("", p.relative_to(root.parent).as_posix().replace("/", "."))
      sys.modules[mod_path] = mod
    return _pytest.pathlib.import_path._before_rdkit_patch(p, mode=mode, root=root)

  setattr(import_path_rdkit_patch, "_before_rdkit_patch", orig_import_path)
  _pytest.pathlib.import_path = import_path_rdkit_patch

import pytest

if __name__ == "__main__":
    raise SystemExit(pytest.console_main())
