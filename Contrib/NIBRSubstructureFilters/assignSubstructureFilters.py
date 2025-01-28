import argparse
import operator
import sys
from collections import Counter, defaultdict, namedtuple

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import FilterCatalog, RDConfig, rdMolDescriptors

FilterMatch = namedtuple(
  'FilterMatch',
  ('SubstructureMatches', 'Min_N_O_filter', 'Frac_N_O', 'Covalent', 'SpecialMol', 'SeverityScore'))


# Build the filter catalog using the RDKit filterCatalog module
def buildFilterCatalog():

  inhousefilter = pd.read_csv(
    f'{RDConfig.RDContribDir}/NIBRSubstructureFilters/SubstructureFilter_HitTriaging_wPubChemExamples.csv'
  )
  inhouseFiltersCat = FilterCatalog.FilterCatalog()
  for i in range(inhousefilter.shape[0]):
    mincount = 1
    if inhousefilter['MIN_COUNT'][i] != 0:
      mincount = int(inhousefilter['MIN_COUNT'][i])
    pname = inhousefilter['PATTERN_NAME'][i]
    sname = inhousefilter['SET_NAME'][i]
    pname_final = '{0}_min({1})__{2}__{3}__{4}'.format(pname, mincount,
                                                       inhousefilter['SEVERITY_SCORE'][i],
                                                       inhousefilter['COVALENT'][i],
                                                       inhousefilter['SPECIAL_MOL'][i])
    fil = FilterCatalog.SmartsMatcher(pname_final, inhousefilter['SMARTS'][i], mincount)
    inhouseFiltersCat.AddEntry(FilterCatalog.FilterCatalogEntry(pname_final, fil))
    inhouseFiltersCat.GetEntry(i).SetProp('Scope', sname)
  return inhouseFiltersCat


# Assign substructure filters and fraction of Nitrogen and Oxygen atoms
def assignFilters(data, nameSmilesColumn='smiles'):

  results = []

  inhouseFiltersCat = buildFilterCatalog()

  NO_filter = '[#7,#8]'
  sma = Chem.MolFromSmarts(NO_filter, mergeHs=True)

  for smi in data[nameSmilesColumn]:
    qc, NO_filter, fracNO, co, sc, sm = [np.nan] * 6

# The following files require numpy and were explicitely checked for their compatibility.

    try:
      mol = Chem.MolFromSmiles(smi)

      # fraction of N and O atoms
      numHeavyAtoms = mol.GetNumHeavyAtoms()
      numNO = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#7,#8]')))
      fracNO = float(numNO) / numHeavyAtoms

      # all substructure filters
      entries = inhouseFiltersCat.GetMatches(mol)
      if len(list(entries)):
        # initialize empty lists
        fs, sev, cov, spm = ([] for _ in range(4))
        # get the matches
        for entry in entries:
          pname = entry.GetDescription()
          n, s, c, m = pname.split('__')
          fs.append(entry.GetProp("Scope") + '_' + n)
          sev.append(int(s))
          cov.append(int(c))
          spm.append(int(m))
        # concatenate all matching filters
        qc = ' | '.join(fs)
        # assign overall severity
        if sev.count(2):
          sc = 10
        else:
          sc = sum(sev)
        # get number of covalent flags and special molecule flags
        co = sum(cov)
        sm = sum(spm)
      # if non of the filters matches
      else:
        qc = 'no match'
        sc = 0
        co = 0
        sm = 0

      # special NO filter
      if not mol.HasSubstructMatch(sma):
        NO_filter = 'no_oxygen_or_nitrogen'
      else:
        NO_filter = 'no match'
    except Exception:
      print("Failed on compound {0}\n".format(smi))
      pass
    results.append(FilterMatch(qc, NO_filter, fracNO, co, sm, sc))
  return results


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('--data', type=str, required=True,
                      help='Please specify the path to your data file. Required format: csv')
  parser.add_argument('--smilesColumn', type=str, required=True,
                      help='Please specify the name of your SMILES column.')
  parser.add_argument('--result', type=str, required=True,
                      help='Please specify the name of your result file.')
  parser.add_argument('--verbose', type=bool, default=1, help='Generate output? Default: False')
  args = parser.parse_args()

  if args.verbose:
    print('---> Reading data')
  datafile = args.data
  try:
    data = pd.read_csv(datafile)
  except Exception:
    if args.verbose:
      print('Data could not be read. Please check your file.')
    sys.exit()

  smiCol = args.smilesColumn

  if args.verbose:
    print('---> Apply filters to data')
  try:
    results = assignFilters(data, nameSmilesColumn=smiCol)
  except Exception:
    if args.verbose:
      print('Smiles column does not exist. Please check.')
    sys.exit()

  df_tmp = pd.DataFrame.from_records(results, columns=FilterMatch._fields)

  data = data.merge(df_tmp, how='left', left_index=True, right_index=True)
  data.to_csv(args.result, index=False)

  if args.verbose:
    print('---> Done')
