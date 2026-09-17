from collections.abc import Iterable, Sequence
from typing import overload

import rdkit.Chem.rdGeneralizedSubstruct
import rdkit.Chem.rdTautomerQuery
import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs


class MolHolderBase:
    """
    Base class for holding molecules used in the Substructure Library.
    Instantiations of this class are passed into the SubstructureLibrary.
    The API is quite simple:
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule
      GetMol(idx) -> return the molecule at index idx
    """

    def __len__(self) -> int: ...

    def AddMol(self, m: rdkit.Chem.rdchem.Mol) -> int:
        """Adds molecule to the molecule holder"""

    def GetMol(self, arg1: int) -> rdkit.Chem.rdchem.Mol | None:
        """
        Returns a particular molecule in the molecule holder

          ARGUMENTS:
            - idx: which molecule to return

            - sanitize: if sanitize is False, return the internal molecule state [default True]

          NOTE: molecule indices start at 0
        """

class MolHolder(MolHolderBase):
    """
    Holds raw in-memory molecules
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule
      GetMol(idx,sanitize=True) -> return the molecule at index idx
    """

    def __init__(self) -> None: ...

class CachedMolHolder(MolHolderBase):
    """
    Holds molecules in their binary representation.
    This allows more molecules to be held in memory at a time
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

      AddBinary(data) -> adds a picked molecule molecule to the molecule holder, returns index of molecule
                         The data is stored as-is, no checking is done for validity.
      GetMol(idx) -> return the molecule at index idx
    """

    def __init__(self) -> None: ...

    def AddBinary(self, pickle: bytes) -> int:
        """
        Add a binary pickle to the molecule holder, no checking is done on the input data
        """

class CachedSmilesMolHolder(MolHolderBase):
    """
    Holds molecules as smiles string
    This allows more molecules to be held in memory at a time
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

      AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule
                           The smiles is stored as-is, no checking is done for validity.
      GetMol(idx) -> return the molecule at index idx
    """

    def __init__(self) -> None: ...

    def AddSmiles(self, smiles: str) -> int:
        """
        Add a trusted smiles string to the molecule holder, no checking is done on the input data
        """

class CachedTrustedSmilesMolHolder(MolHolderBase):
    """
    Holds molecules as trusted smiles string
    This allows more molecules to be held in memory at a time and avoids RDKit sanitization
    overhead.
    See: http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

      AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule
                           The smiles is stored as-is, no checking is done for validity.
      GetMol(idx,s) -> return the molecule at index idx,
                  note, only light sanitization is done here, for instance
                  the molecules RingInfo is not initialized
    """

    def __init__(self) -> None: ...

    def AddSmiles(self, smiles: str) -> int:
        """
        Add a trusted smiles string to the molecule holder, no checking is done on the input data
        """

class FPHolderBase:
    def __len__(self) -> int: ...

    def AddMol(self, m: rdkit.Chem.rdchem.Mol) -> int:
        """
        Adds a molecule to the fingerprint database, returns the index of the new pattern
        """

    def AddFingerprint(self, v: rdkit.DataStructs.cDataStructs.ExplicitBitVect) -> int:
        """
        Adds a raw bit vector to the fingerprint database, returns the index of the supplied pattern
        """

    def GetFingerprint(self, idx: int) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
        """Return the bit vector at the specified index"""

    def PassesFilter(self, idx: int, query: rdkit.DataStructs.cDataStructs.ExplicitBitVect) -> bool:
        """
        Returns True if the specified index passes the filter supplied by the query bit vector
        """

    def MakeFingerprint(self, mol: rdkit.Chem.rdchem.Mol) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
        """Compute the query bits for the holder"""

class PatternHolder(FPHolderBase):
    """
    Holds fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, numBits: int) -> None: ...

class TautomerPatternHolder(FPHolderBase):
    """
    Holds tautomeric fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.
    These fingerprints are designed to be used with TautomerQueries.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, numBits: int) -> None: ...

class KeyHolderBase:
    def __len__(self) -> int: ...

    def AddMol(self, m: rdkit.Chem.rdchem.Mol) -> int:
        """
        Adds a molecule to the fingerprint database, returns the index of the new pattern
        """

    def AddKey(self, arg1: str) -> int:
        """Add a key to the key holder, must be manually synced"""

    def GetKey(self, arg1: int) -> str:
        """Return the key at the specified index"""

    def GetKeys(self, indices: Sequence[int]) -> list[str]:
        """
        Returns the keys for the given indices as return by GetMatches

          ARGUMENTS:
            - indices: The indices of the keys
        """

class KeyFromPropHolder(KeyHolderBase):
    """
    Holds keys to return external references to the molecules in the molholder.
    By default use the _Name property but can be overridden to be any property
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, propname: str) -> None: ...

    def GetPropName(self) -> str:
        """Return the key for the given molecule index"""

class SubstructLibrary:
    """
    SubstructLibrary: This provides a simple API for substructure searching large datasets
    The SubstructLibrary takes full advantage of available threads during the search operation.
    Basic operation is simple

    >>> from __future__ import print_function
    >>> import os
    >>> from rdkit import Chem, RDConfig
    >>> from rdkit.Chem import rdSubstructLibrary
    >>> library = rdSubstructLibrary.SubstructLibrary()
    >>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir,
    ...                               'NCI', 'first_200.props.sdf')):
    ...   idx = library.AddMol(mol)
    >>> core = Chem.MolFromSmarts('CCCCOC')
    >>> indices = library.GetMatches(core)
    >>> len(indices)
    11

    Substructure matching options can be sent into GetMatches:

    >>> indices = library.GetMatches(core, useChirality=False)
    >>> len(indices)
    11

    Controlling the number of threads or the maximum number of matches returned:
    is also available (the default is to run on all cores)

    >>> indices = library.GetMatches(core, numThreads=2, maxResults=10)
    >>> len(indices)
    10

    Working on larger datasets:

    Molecules are fairly large objects and will limit the number that can be kept in memory.
    To assist this we supply three other molecule holders:
      CachedMolHolder - stores molecules as their pickled representation

      CachedSmilesMolHolder - stores molecules internally as smiles strings

      CachedTrustedSmilesMolHolder = excepts (and stores) molecules as trusted smiles strings

    Using Pattern fingerprints as a pre-filter:
    Pattern fingerprints provide an easy way to indicate whether the substructure search should be
    be done at all.  This is particularly useful with the Binary and Smiles based molecule holders
    as they have an expensive molecule creation step in addition to the substructure searching step

    >>> library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.CachedSmilesMolHolder(),
    ...                                               rdSubstructLibrary.PatternHolder())
    >>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir,
    ...                               'NCI', 'first_200.props.sdf')):
    ...   idx = library.AddMol(mol)
    >>> indices = library.GetMatches(core)
    >>> len(indices)
    11

    This (obviously) takes longer to initialize.  However, both the molecule and pattern
    holders can be populated with raw data, a simple example is below:

    >>> import csv
    >>> molholder = rdSubstructLibrary.CachedSmilesMolHolder()
    >>> pattern_holder = rdSubstructLibrary.PatternHolder()
    >>> with open(os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.tpsa.csv')) as inf:
    ...   for i, row in enumerate(csv.reader(inf)):
    ...     if i:
    ...       idx = molholder.AddSmiles(row[0])
    ...       idx2 = pattern_holder.AddFingerprint(
    ...           pattern_holder.MakeFingerprint(Chem.MolFromSmiles(row[0])))
    ...       assert idx==idx2
    >>> library = rdSubstructLibrary.SubstructLibrary(molholder,pattern_holder)
    >>> indices = library.GetMatches(core)
    >>> len(indices)
    11

    Finally, the KeyFromPropHolder can be used to use external keys such as
    compound names.  By default the holder uses the '_Name' property but can
    be changed to any property.

    >>> library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.MolHolder(), rdSubstructLibrary.KeyFromPropHolder())
    >>> m = Chem.MolFromSmiles('CCC')
    >>> m.SetProp('_Name', 'Z11234')
    >>> idx = library.AddMol(m)
    >>> indices = library.GetMatches(m)
    >>> list(library.GetKeyHolder().GetKeys(indices))
    ['Z11234']
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, arg: bytes, /) -> None: ...

    @overload
    def __init__(self, arg: str, /) -> None: ...

    @overload
    def __init__(self, molecules: object) -> None: ...

    @overload
    def __init__(self, molecules: object, fingerprints_or_keys: object | None) -> None: ...

    @overload
    def __init__(self, molecules: object, fingerprints: object | None, keys: object | None) -> None: ...

    def GetMolHolder(self) -> MolHolderBase: ...

    def GetFpHolder(self) -> FPHolderBase: ...

    def GetKeyHolder(self) -> KeyHolderBase: ...

    def AddMol(self, mol: rdkit.Chem.rdchem.Mol) -> int:
        """Adds a molecule to the substruct library"""

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.Mol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.Mol, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.Mol, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.MolBundle, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.MolBundle, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.MolBundle, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdchem.MolBundle, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> list[int]: ...

    @overload
    def GetMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]:
        """
        Get the matches for the query.

         Arguments:
          - query:      substructure query
          - numThreads: number of threads to use, -1 means all threads
          - maxResults: maximum number of results to return
        """

    @overload
    def GetMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> list[int]:
        """
        Get the matches for the query.

         Arguments:
          - query:      substructure query
          - startIdx:   index to search from
          - endIdx:     index (non-inclusize) to search to
          - numThreads: number of threads to use, -1 means all threads
          - maxResults: maximum number of results to return
        """

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.Mol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.Mol, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.Mol, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.MolBundle, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.MolBundle, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.MolBundle, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdchem.MolBundle, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...

    @overload
    def CountMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int:
        """
        Get the matches for the query.

         Arguments:
          - query:      substructure query
          - numThreads: number of threads to use, -1 means all threads
        """

    @overload
    def CountMatches(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> int:
        """
        Get the matches for the query.

         Arguments:
          - query:      substructure query
          - startIdx:   index to search from
          - endIdx:     index (non-inclusize) to search to
          - numThreads: number of threads to use, -1 means all threads
        """

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.Mol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.Mol, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.Mol, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdTautomerQuery.TautomerQuery, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.MolBundle, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.MolBundle, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.MolBundle, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdchem.MolBundle, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...

    @overload
    def HasMatch(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool:
        """
        Get the matches for the query.

         Arguments:
          - query:      substructure query
          - numThreads: number of threads to use, -1 means all threads
        """

    @overload
    def HasMatch(self, query: rdkit.Chem.rdGeneralizedSubstruct.ExtendedQueryMol, startIdx: int, endIdx: int, parameters: rdkit.Chem.rdchem.SubstructMatchParameters, numThreads: int = -1) -> bool:
        """
        Get the matches for the query.

         Arguments:
          - query:      substructure query
          - startIdx:   index to search from
          - endIdx:     index (non-inclusize) to search to
          - numThreads: number of threads to use, -1 means all threads
        """

    def GetMol(self, idx: int) -> rdkit.Chem.rdchem.Mol | None:
        """
        Returns a particular molecule in the molecule holder

          ARGUMENTS:
            - idx: which molecule to return

          NOTE: molecule indices start at 0
        """

    def SetSearchOrder(self, seq: Iterable[int] | None) -> None:
        """
        Sets the search order for the library

          ARGUMENTS:
            - order: sequence of molecule indices

          NOTE: molecule indices start at 0
        """

    def GetSearchOrder(self) -> tuple:
        """
        Returns the search order for the library

          NOTE: molecule indices start at 0
        """

    def __len__(self) -> int: ...

    def ToStream(self, stream: object) -> None:
        """
        Serialize a substructure library to a python text stream.
        The stream can be a file in text mode or an io.StringIO type object

          ARGUMENTS:
            - stream: a text or text stream like object

          >>> from rdkit.Chem import rdSubstructLibrary
          >>> import io
          >>> lib = rdSubstructLibrary.SubstructLibrary()
          >>> stream = io.StringIO()
          >>> lib.ToStream(stream)

           or
          >>> with open('rdkit.sslib', 'w') as stream:
          ...  lib.ToStream(stream)
        """

    def InitFromStream(self, stream: object) -> None:
        """
        Deserialize a substructure library from a python bytes stream.
        Python doesn't allow seeking operations inside a unicode or string stream anymore
        so this requires opening a file in binary mode or using an io.ByteIO type object

          ARGUMENTS:
            - stream: a binary stream like object

          SubstructLibrary.Serialize already writes a binary stream

          >>> from rdkit.Chem import rdSubstructLibrary
          >>> import io
          >>> lib = rdSubstructLibrary.SubstructLibrary()
          >>> stream = io.BytesIO( lib.Serialize() )
          >>> lib.InitFromStream(stream)

           remember to write to text and read from a binary stream
          >>> with open('rdkit.sslib', 'w') as f: lib.ToStream(f)
          >>> with open('rdkit.sslib', 'rb') as f: lib.InitFromStream(f)
        """

    def Serialize(self) -> bytes: ...

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

def SubstructLibraryCanSerialize() -> bool:
    """
    Returns True if the SubstructLibrary is serializable (requires boost serialization
    """

@overload
def AddPatterns(sslib: SubstructLibrary, numThreads: int = 1) -> None: ...

@overload
def AddPatterns(sslib: SubstructLibrary, patterns: FPHolderBase, numThreads: int = 1) -> None:
    """
    Add pattern fingerprints to the given library, use numThreads=-1 to use all available cores
    """
