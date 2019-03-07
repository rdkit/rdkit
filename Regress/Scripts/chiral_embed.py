
from rdkit import Chem
from rdkit.Chem.PyMol import MolViewer
from rdkit.Chem import AllChem
import gzip

try:
  v = MolViewer()
  v.DeleteAll()
except Exception:
  v = None

outf = open('chembl_20_chiral.problems.smi', 'w+')
noconff = open('chembl_20_chiral.noconfs.smi', 'w+')
for i, line in enumerate(gzip.open('../Data/chembl_20_chiral.smi.gz')):
  line = line.strip().decode().split(' ')
  mol = Chem.MolFromSmiles(line[0])
  if not mol:
    continue
  cents = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
  if len([y for x, y in cents if y == '?']):
    continue
  nm = line[1]
  csmi = Chem.MolToSmiles(mol, True)
  for j in range(100):
    mh = Chem.AddHs(mol)
    ok = AllChem.EmbedMolecule(mh, randomSeed=j + 1)
    if ok >= 0:
      Chem.AssignAtomChiralTagsFromStructure(mh)
      newm = Chem.RemoveHs(mh)
      smi = Chem.MolToSmiles(newm, True)
      if smi != csmi:
        print('%d %d %s:\n%s\n%s' % (i, j, nm, csmi, smi))
        print('%s %s %d' % (line[0], line[1], j + 1), file=outf)

        if v is not None:
          v.ShowMol(mh, name='%s-%d' % (nm, j), showOnly=False)
          break  # move immediately onto the next molecule
    else:
      print('noconf %d %d %s: %s' % (i, j, nm, line[0]))
      print('%s %s %d' % (line[0], line[1], j + 1), file=noconff)

  print('Done with mol %d' % i)
