from rdkit import Chem
from rdkit.Chem.GetEnantiomer import GetEnantiomer

def GetChirality(m):
  """ returns the chirality of a molecule. Return True if the molecule is chiral and False if it is not. 

   Arguments:
      - m: the molecule to work with

    >>> from rdkit import Chem
    >>> from rdkit.Chem.GetChirality import GetChirality
    >>> m = Chem.MolFromSmiles('CN1CCC[C@H]1C2=CN=CC=C2') # L-(-)-Nicotine
    >>> GetChirality(m)
    True

    If a molecule is achiral (not chiral), the function will return False

    >>> m = Chem.MolFromSmiles('[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O') # meso-Tartaric Acid
    >>> GetChirality(m)
    False

  """

  return GetEnantiomer(m) is not None
