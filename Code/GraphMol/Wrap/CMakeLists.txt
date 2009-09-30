rdkit_python_extension(rdchem 
  rdchem.cpp Table.cpp Atom.cpp Bond.cpp Mol.cpp
  Conformer.cpp RingInfo.cpp EditableMol.cpp 
  DEST Chem 
  LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol 
     RDGeometryLib RDGeneral RDBoost )

rdkit_python_extension(rdmolops
  rdmolops.cpp MolOps.cpp
  DEST Chem 
  LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol 
     RDGeometryLib RDGeneral RDBoost )

rdkit_python_extension(rdmolfiles
  rdmolfiles.cpp SDMolSupplier.cpp 
  SmilesMolSupplier.cpp SmilesWriter.cpp SDWriter.cpp
  DEST Chem 
  LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol 
     RDGeometryLib RDGeneral RDBoost )

add_test( pyGraphMolWrap python ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py )
add_test( pyTestConformerWrap python ${CMAKE_CURRENT_SOURCE_DIR}/testConformer.py )
