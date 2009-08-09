add_library(GraphMol SHARED Atom.cpp QueryAtom.cpp QueryBond.cpp Bond.cpp MolOps.cpp
  FindRings.cpp ROMol.cpp RWMol.cpp PeriodicTable.cpp atomic_data.cpp
  QueryOps.cpp MolPickler.cpp Canon.cpp AtomIterators.cpp BondIterators.cpp
  Aromaticity.cpp Kekulize.cpp MolDiscriminators.cpp ConjugHybrid.cpp
  AddHs.cpp RankAtoms.cpp Matrices.cpp Chirality.cpp RingInfo.cpp
  Conformer.cpp )
INSTALL(TARGETS GraphMol LIBRARY DESTINATION ${RDKit_LibDir})

add_subdirectory(Depictor)
add_subdirectory(FileParsers)
add_subdirectory(SmilesParse)
add_subdirectory(Substruct)
add_subdirectory(ChemReactions)
add_subdirectory(ChemTransforms)
add_subdirectory(Wrap)

rdkit_test(graphmolTest1 test1.cpp  LINK_LIBRARIES 
  GraphMol RDGeneral RDGeometry)
rdkit_test(graphmolcpTest cptest.cpp  LINK_LIBRARIES 
  GraphMol RDGeneral RDGeometry)
rdkit_test(graphmolqueryTest querytest.cpp  LINK_LIBRARIES
   GraphMol RDGeneral RDGeometry)
rdkit_test(graphmolMolOpsTest molopstest.cpp  LINK_LIBRARIES 
  GraphMol FileParsers SmilesParse RDGeneral RDGeometry)
rdkit_test(graphmoltestCanon testCanon.cpp  LINK_LIBRARIES 
  GraphMol SubstructMatch SmilesParse RDGeneral RDGeometry)
rdkit_test(graphmoltestChirality testChirality.cpp  LINK_LIBRARIES 
  GraphMol FileParsers SmilesParse RDGeneral RDGeometry)
rdkit_test(graphmoltestPickler testPickler.cpp  LINK_LIBRARIES 
  GraphMol SubstructMatch FileParsers SmilesParse RDGeneral RDGeometry)
rdkit_test(graphmolIterTest itertest.cpp  LINK_LIBRARIES 
  SmilesParse GraphMol RDGeneral RDGeometry)


