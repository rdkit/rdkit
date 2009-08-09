rdkit_python_extension(rdchem 
  rdchem.cpp Table.cpp Atom.cpp Bond.cpp Mol.cpp
  Conformer.cpp RingInfo.cpp EditableMol.cpp 
  DEST Chem 
  LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol 
     RDGeometry DataStructs RDGeneral RDBoost )

rdkit_python_extension(rdmolops
  rdmolops.cpp MolOps.cpp
  DEST Chem 
  LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol 
     RDGeometry DataStructs RDGeneral RDBoost )

rdkit_python_extension(rdmolfiles
  rdmolfiles.cpp SDMolSupplier.cpp TDTMolSupplier.cpp
  SmilesMolSupplier.cpp SmilesWriter.cpp SDWriter.cpp TDTWriter.cpp  
  DEST Chem 
  LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol 
     RDGeometry DataStructs RDGeneral RDBoost )

