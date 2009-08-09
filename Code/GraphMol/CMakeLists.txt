add_library(GraphMol Atom.cpp QueryAtom.cpp QueryBond.cpp Bond.cpp MolOps.cpp
  FindRings.cpp ROMol.cpp RWMol.cpp PeriodicTable.cpp atomic_data.cpp
  QueryOps.cpp MolPickler.cpp Canon.cpp AtomIterators.cpp BondIterators.cpp
  Aromaticity.cpp Kekulize.cpp MolDiscriminators.cpp ConjugHybrid.cpp
  AddHs.cpp RankAtoms.cpp Matrices.cpp Chirality.cpp RingInfo.cpp
  Conformer.cpp )

add_executable(test1 test1.cpp)
target_link_libraries(test1  GraphMol RDGeneral RDGeometry )
add_test(test1 ${EXECUTABLE_OUTPUT_PATH}/test1 )


add_subdirectory(Depictor)
add_subdirectory(FileParsers)
add_subdirectory(SmilesParse)
add_subdirectory(Substruct)
add_subdirectory(ChemReactions)

add_subdirectory(Wrap)


