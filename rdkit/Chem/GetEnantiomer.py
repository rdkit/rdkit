from rdkit import Chem
from rdkit.Chem import AllChem

def GetEnantiomer(m):
  """ returns the enantiomer of a molecule. If an enantiomer does not exist, return None

   Arguments:
      - m: the molecule to work with

    >>> from rdkit import Chem
    >>> from rdkit.Chem.GetEnantiomer import GetEnantiomer
    >>> m = Chem.MolFromSmiles('CN1CCC[C@H]1C2=CN=CC=C2') # L-(-)-Nicotine
    >>> Chem.MolToSmiles(GetEnantiomer(m)) # (+)-Nicotine
    CN1CCC[C@@H]1C2=CN=CC=C2

    If an enantiomer does not exist (i.e. the molecule is not chiral), None will be returned

    >>> m = Chem.MolFromSmiles('[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O') # meso-Tartaric Acid
    >>> GetEnantiomer(m)
    None

  """

  
  m1 = Chem.AddHs(m)
  AllChem.EmbedMolecule(m1) # Make sure molecule is in 3D
  m_b = Chem.MolToMolBlock(m1)

  # Reflect the molecule across the origin
  m_b = m_b.replace("-", "+").replace("    ", "   -").replace("+", " ")
  m2 = Chem.MolFromMolBlock(m_b)

  # If an enantiomer does not exist, return None
  if Chem.MolToSmiles(m1) == Chem.MolToSmiles(Chem.AddHs(m2)):
      return None

  return m2
