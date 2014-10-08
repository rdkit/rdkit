rdkit_library(GraphMol 
              Atom.cpp QueryAtom.cpp QueryBond.cpp Bond.cpp 
              MolOps.cpp FindRings.cpp ROMol.cpp RWMol.cpp PeriodicTable.cpp 
              atomic_data.cpp QueryOps.cpp MolPickler.cpp Canon.cpp 
              AtomIterators.cpp BondIterators.cpp Aromaticity.cpp Kekulize.cpp 
              MolDiscriminators.cpp ConjugHybrid.cpp AddHs.cpp RankAtoms.cpp 
              Matrices.cpp Chirality.cpp RingInfo.cpp Conformer.cpp
              Renumber.cpp
              SHARED 
              LINK_LIBRARIES RDGeometryLib RDGeneral 
                 ${RDKit_THREAD_LIBS})

rdkit_headers(Atom.h
              atomic_data.h
              AtomIterators.h
              Bond.h
              BondIterators.h
              Canon.h
              Chirality.h
              Conformer.h
              GraphMol.h
              MolOps.h
              MolPickler.h
              PeriodicTable.h
              QueryAtom.h
              QueryBond.h
              QueryOps.h
              RankAtoms.h
              RDKitBase.h
              RDKitQueries.h
              RingInfo.h
              Rings.h
              ROMol.h
              RWMol.h
              SanitException.h
              MonomerInfo.h
              DEST GraphMol)

add_subdirectory(Depictor)
add_subdirectory(SmilesParse)
add_subdirectory(FileParsers)
add_subdirectory(Substruct)
add_subdirectory(ChemReactions)
add_subdirectory(ChemTransforms)

add_subdirectory(Subgraphs)
add_subdirectory(FragCatalog)
add_subdirectory(Descriptors)

add_subdirectory(Fingerprints)
add_subdirectory(PartialCharges)

add_subdirectory(MolTransforms)
add_subdirectory(ForceFieldHelpers)
add_subdirectory(DistGeomHelpers)

add_subdirectory(MolAlign)
add_subdirectory(MolChemicalFeatures)
add_subdirectory(ShapeHelpers)
add_subdirectory(MolCatalog)

add_subdirectory(MolDrawing)

add_subdirectory(FMCS)



if(RDK_BUILD_SLN_SUPPORT AND NOT CMAKE_COMPILER_IS_CLANG)
add_subdirectory(SLNParse)
endif()

add_subdirectory(Wrap)

rdkit_test(graphmolTest1 test1.cpp LINK_LIBRARIES GraphMol RDGeometryLib RDGeneral)

rdkit_test(graphmolcpTest cptest.cpp LINK_LIBRARIES  SmilesParse GraphMol RDGeometryLib RDGeneral)

rdkit_test(graphmolqueryTest querytest.cpp LINK_LIBRARIES SubstructMatch SmilesParse GraphMol RDGeometryLib RDGeneral)

rdkit_test(graphmolMolOpsTest molopstest.cpp 
           LINK_LIBRARIES SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral)

rdkit_test(graphmoltestCanon testCanon.cpp  
           LINK_LIBRARIES SubstructMatch SmilesParse GraphMol RDGeometryLib RDGeneral)

rdkit_test(graphmoltestChirality testChirality.cpp  
           LINK_LIBRARIES SmilesParse FileParsers GraphMol RDGeometryLib RDGeneral)

rdkit_test(graphmoltestPickler testPickler.cpp  
           LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol RDGeometryLib RDGeneral)

rdkit_test(graphmolIterTest itertest.cpp LINK_LIBRARIES SmilesParse GraphMol RDGeometryLib RDGeneral)


