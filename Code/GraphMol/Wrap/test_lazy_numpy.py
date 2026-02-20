#
#  Copyright (C) 2025 RDKit contributors
#         All Rights Reserved
#
"""Tests that 'from rdkit import Chem' does not eagerly load numpy,
and that core Chem functionality works before numpy is loaded.

Because numpy's import state is process-global and cannot be unloaded,
import-ordering tests must run in a **fresh subprocess**.
"""
import subprocess
import sys
import textwrap
import unittest


def _run_snippet(code: str) -> subprocess.CompletedProcess:
  """Run *code* in a clean Python subprocess and return the result."""
  return subprocess.run(
    [sys.executable, "-c", textwrap.dedent(code)],
    capture_output=True,
    text=True,
    timeout=60,
  )


class TestLazyNumpy(unittest.TestCase):

  def test_chem_import_does_not_load_numpy(self):
    """Importing rdkit.Chem must not pull numpy into sys.modules."""
    result = _run_snippet("""\
      import sys
      from rdkit import Chem
      # numpy must not have been imported as a side-effect
      if "numpy" in sys.modules:
          sys.exit("FAIL: numpy was loaded by 'from rdkit import Chem'")
    """)
    self.assertEqual(result.returncode, 0, result.stderr or result.stdout)

  def test_basic_smiles_roundtrip_without_numpy(self):
    """MolFromSmiles / MolToSmiles must work before numpy is loaded."""
    result = _run_snippet("""\
      import sys
      from rdkit import Chem

      mol = Chem.MolFromSmiles("c1ccccc1")
      assert mol is not None, "MolFromSmiles returned None"

      smi = Chem.MolToSmiles(mol)
      assert smi == "c1ccccc1", f"unexpected SMILES: {smi}"

      assert "numpy" not in sys.modules, "numpy crept in during SMILES ops"
    """)
    self.assertEqual(result.returncode, 0, result.stderr or result.stdout)

  def test_mol_operations_without_numpy(self):
    """Core mol operations (AddHs, sanitize, substruct) work pre-numpy."""
    result = _run_snippet("""\
      import sys
      from rdkit import Chem

      mol = Chem.MolFromSmiles("CCO")
      assert mol is not None

      # AddHs / RemoveHs
      molH = Chem.AddHs(mol)
      assert molH.GetNumAtoms() == 9  # 3 heavy + 6 H
      mol2 = Chem.RemoveHs(molH)
      assert mol2.GetNumAtoms() == 3

      # Substructure match
      query = Chem.MolFromSmarts("[OH]")
      assert mol.HasSubstructMatch(query)

      # MolToMolBlock (V2000)
      mb = Chem.MolToMolBlock(mol)
      assert "V2000" in mb

      # InChI round-trip
      inchi = Chem.MolToInchi(mol)
      assert inchi.startswith("InChI=")

      assert "numpy" not in sys.modules, "numpy crept in during mol ops"
    """)
    self.assertEqual(result.returncode, 0, result.stderr or result.stdout)

  def test_numpy_loads_on_demand(self):
    """Calling a numpy-returning API (GetDistanceMatrix) loads numpy lazily."""
    result = _run_snippet("""\
      import sys
      from rdkit import Chem

      assert "numpy" not in sys.modules, "numpy loaded too early"

      mol = Chem.MolFromSmiles("CCO")
      dm = Chem.GetDistanceMatrix(mol)

      assert "numpy" in sys.modules, "numpy should be loaded after GetDistanceMatrix"
      assert dm.shape == (3, 3), f"unexpected shape: {dm.shape}"
      assert dm[0, 2] == 2.0, f"unexpected distance: {dm[0, 2]}"
    """)
    self.assertEqual(result.returncode, 0, result.stderr or result.stdout)

  def test_conformer_positions_loads_numpy(self):
    """GetPositions() / SetPositions() load numpy on first call."""
    result = _run_snippet("""\
      import sys
      from rdkit import Chem
      from rdkit.Chem import rdchem

      assert "numpy" not in sys.modules, "numpy loaded too early"

      mol = Chem.MolFromSmiles("C")
      conf = rdchem.Conformer(mol.GetNumAtoms())
      conf.SetAtomPosition(0, (1.0, 2.0, 3.0))
      mol.AddConformer(conf, assignId=True)

      pos = mol.GetConformer().GetPositions()
      assert "numpy" in sys.modules, "numpy should load after GetPositions"
      assert pos.shape == (1, 3), f"unexpected shape: {pos.shape}"
      assert abs(pos[0, 0] - 1.0) < 1e-6
    """)
    self.assertEqual(result.returncode, 0, result.stderr or result.stdout)

  def test_datastructs_convert_loads_numpy(self):
    """DataStructs.ConvertToNumpyArray loads numpy on demand."""
    result = _run_snippet("""\
      import sys
      from rdkit import Chem, DataStructs

      assert "numpy" not in sys.modules, "numpy loaded too early"

      fp = Chem.RDKFingerprint(Chem.MolFromSmiles("c1ccccc1"))

      import numpy as np
      arr = np.zeros(len(fp), dtype=int)
      DataStructs.ConvertToNumpyArray(fp, arr)
      assert arr.sum() > 0, "fingerprint should have bits set"
    """)
    self.assertEqual(result.returncode, 0, result.stderr or result.stdout)

  def test_adjacency_matrix_loads_numpy(self):
    """GetAdjacencyMatrix loads numpy on demand and returns correct result."""
    result = _run_snippet("""\
      import sys
      from rdkit import Chem

      assert "numpy" not in sys.modules, "numpy loaded too early"

      mol = Chem.MolFromSmiles("CC")
      adj = Chem.GetAdjacencyMatrix(mol)

      assert "numpy" in sys.modules, "numpy should load after GetAdjacencyMatrix"
      assert adj[0, 1] == 1
      assert adj[1, 0] == 1
      assert adj[0, 0] == 0
    """)
    self.assertEqual(result.returncode, 0, result.stderr or result.stdout)


if __name__ == '__main__':
  unittest.main()
