import csv
import io
import logging
import os

import freewilson as fw

from rdkit import Chem, rdBase

PATH = os.path.join(os.path.dirname(fw.__file__), 'data')
assert os.path.exists(PATH), PATH


def test_chembl():
  logging.getLogger().setLevel(logging.INFO)
  smilesfile = os.path.join(PATH, "CHEMBL2321810.smi")
  scaffoldfile = os.path.join(PATH, "CHEMBL2321810_scaffold.mol")
  csvfile = os.path.join(PATH, "CHEMBL2321810_act.csv")
  assert os.path.exists(smilesfile)
  mols = []
  for line in open(smilesfile):
    smiles, name = line.strip().split()
    m = Chem.MolFromSmiles(smiles)
    m.SetProp("_Name", name)
    mols.append(m)

  scaffold = Chem.MolFromMolBlock(open(scaffoldfile).read())
  data = {k: float(v) for k, v in list(csv.reader(open(csvfile)))[1:]}

  scores = [data[m.GetProp("_Name")] for m in mols]
  assert mols and len(mols) == len(scores)

  with rdBase.BlockLogs():
    free = fw.FWDecompose(scaffold, mols, scores)
  # let's make sure the r squared is decent
  assert free.r2 > 0.8

  # assert we get something
  preds = list(fw.FWBuild(free))
  assert len(preds)

  # check to see that the prediction filters work
  preds2 = list(fw.FWBuild(free, pred_filter=lambda x: x > 8))
  assert len(preds2)
  assert len([p for p in preds if p.prediction > 8]) == len(list(preds2))

  # check to see that the R groups are output in order, i.e. R10 after R3
  s = io.StringIO()
  fw.predictions_to_csv(s, free, preds2)
  assert s.getvalue()

  s2 = io.StringIO(s.getvalue())
  for i, row in enumerate(csv.reader(s2)):
    if i == 0:
      assert row == ['smiles', 'prediction', 'Core_smiles', 'R1_smiles', 'R3_smiles', 'R10_smiles']
  assert i > 0


def test_multicore():
  # test that we can add rgroups for later cores and not throw an exception
  scaffolds = [Chem.MolFromSmiles("c1ccccc1[*].NC=O"), Chem.MolFromSmiles("C1CCCCC1")]
  mols = [
    Chem.MolFromSmiles(x) for x in [
      'c1ccccc1CC2CNC2C(=O)N', 'Cc1ccccc1CC2CNC2C(=O)N', 'Cc1ccccc1CC2CNCC(=O)NC2',
      'C3c1ccccc1CC2CNC2C(=O)N3', 'C1CCCCC1F', 'ClC1CCCCC1F'
    ]
  ]
  decomp = fw.FWDecompose(scaffolds, mols, [1, 2, 3, 4, 5, 6])
  s = io.StringIO()
  fw.predictions_to_csv(s, decomp, fw.FWBuild(decomp))
