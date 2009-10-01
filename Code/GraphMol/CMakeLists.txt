set(lib_srcs Atom.cpp QueryAtom.cpp QueryBond.cpp Bond.cpp MolOps.cpp
  FindRings.cpp ROMol.cpp RWMol.cpp PeriodicTable.cpp atomic_data.cpp
  QueryOps.cpp MolPickler.cpp Canon.cpp AtomIterators.cpp BondIterators.cpp
  Aromaticity.cpp Kekulize.cpp MolDiscriminators.cpp ConjugHybrid.cpp
  AddHs.cpp RankAtoms.cpp Matrices.cpp Chirality.cpp RingInfo.cpp
  Conformer.cpp )

if(MSVC)
add_library(GraphMol ${lib_srcs} )
else(MSVC)
add_library(GraphMol SHARED ${lib_srcs} )
target_link_libraries(GraphMol RDGeometryLib RDGeneral )
endif(MSVC)
INSTALL(TARGETS GraphMol 
  DESTINATION ${RDKit_LibDir}
)

add_subdirectory(Depictor)
add_subdirectory(FileParsers)
add_subdirectory(SmilesParse)
add_subdirectory(Substruct)
add_subdirectory(ChemReactions)
add_subdirectory(ChemTransforms)
add_subdirectory(Wrap)

rdkit_test(graphmolTest1 test1.cpp  LINK_LIBRARIES 
  GraphMol RDGeneral RDGeometryLib)
rdkit_test(graphmolcpTest cptest.cpp  LINK_LIBRARIES 
  GraphMol RDGeneral RDGeometryLib)
rdkit_test(graphmolqueryTest querytest.cpp  LINK_LIBRARIES
   GraphMol RDGeneral RDGeometryLib)
rdkit_test(graphmolMolOpsTest molopstest.cpp  LINK_LIBRARIES 
  FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib)
rdkit_test(graphmoltestCanon testCanon.cpp  LINK_LIBRARIES 
  SubstructMatch SmilesParse GraphMol RDGeneral RDGeometryLib)
rdkit_test(graphmoltestChirality testChirality.cpp  LINK_LIBRARIES 
   FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib)
rdkit_test(graphmoltestPickler testPickler.cpp  LINK_LIBRARIES 
  SubstructMatch FileParsers SmilesParse GraphMol RDGeneral RDGeometryLib)
rdkit_test(graphmolIterTest itertest.cpp  LINK_LIBRARIES 
  SmilesParse GraphMol RDGeneral RDGeometryLib)


