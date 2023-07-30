from rdkit import Chem
from rdkit.Chem import AllChem
from logging import warning

def _GetEnantiomer(m, forceHs: bool = False, silent: bool = False):
   m1 = Chem.AddHs(m) if forceHs else rdchem.Mol(m)
   AllChem.EmbedMolecule(m1) # Make sure molecule is in 3D
   m_b = Chem.MolToMolBlock(m1)

   # Reflect the molecule across the origin
   m_b = m_b.replace("-", "+").replace("    ", "   -").replace("+", " ")
   m2 = Chem.MolFromMolBlock(m_b)

   # If an enantiomer does not exist, return None
   m1_str = Chem.MolToInchi(m1)
   m2_str = Chem.MolToInchi(Chem.AddHs(m2)) if forceHs else Chem.MolToSmiles(m2)

   if m1_str == m2_str:
      if not silent:
         warning(" No enantiomer found for molecule")
      return None
   elif m1_str != m2_str:
      return m2

def IsMolChiral(m, forceHs: bool = True) -> bool:
   """ returns whether the molecule is chiral. Return True if the molecule is chiral and False if it is not. 

   Arguments:
      - m: the molecule to work with
      - forceHs: force the addition of hydrogens when constructing the enantiomer

    >>> from rdkit import Chem
    >>> from rdkit.Chem.MolChiral import IsMolChiral
    >>> m = Chem.MolFromSmiles('CN1CCC[C@H]1C2=CN=CC=C2') # L-(-)-Nicotine
    >>> IsMolChiral(m)
    True

    If a molecule is achiral (not chiral), the function will return False

    >>> m = Chem.MolFromSmiles('[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O') # meso-Tartaric Acid
    >>> IsMolChiral(m)
    False

   """
   return _GetEnantiomer(m, forceHs=forceHs, silent=True) is not None

def GetEnantiomer(m, forceHs: bool = True, silent=False):
   """ returns the enantiomer of a molecule. If an enantiomer does not exist, return None

   Arguments:
      - m: the molecule to work with
      - forceHs: force the addition of hydrogens when constructing the enantiomer
      - silent: print warnings

    >>> from rdkit import Chem
    >>> from rdkit.Chem.MolChiral import GetEnantiomer
    >>> m = Chem.MolFromSmiles('CN1CCC[C@H]1C2=CN=CC=C2') # L-(-)-Nicotine
    >>> Chem.MolToSmiles(GetEnantiomer(m)) # (+)-Nicotine
    CN1CCC[C@@H]1C2=CN=CC=C2

    If an enantiomer does not exist (i.e. the molecule is not chiral), None will be returned

    >>> m = Chem.MolFromSmiles('[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O') # meso-Tartaric Acid
    >>> GetEnantiomer(m)
    None

   """

   return _GetEnantiomer(m, forceHs=forceHs, silent=silent)
