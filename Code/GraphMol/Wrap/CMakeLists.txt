remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
rdkit_python_extension(rdchem
                       rdchem.cpp Table.cpp Atom.cpp Bond.cpp Mol.cpp
                       Conformer.cpp RingInfo.cpp EditableMol.cpp MonomerInfo.cpp
                       ResonanceMolSupplier.cpp MolBundle.cpp StereoGroup.cpp SubstanceGroup.cpp
                       Chirality.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       SmilesParse ChemTransforms SubstructMatch GraphMol
                       MarvinParser)

rdkit_python_extension(rdmolops
                       rdmolops.cpp MolOps.cpp
                       ChiralityOps.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ChemReactions Depictor
                       FileParsers SubstructMatch Fingerprints ChemTransforms
                       Subgraphs SmilesParse MolTransforms GraphMol)

rdkit_python_extension(rdqueries
                       rdqueries.cpp Queries.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       GraphMol)

if(RDK_BUILD_MAEPARSER_SUPPORT)
    set (maesupplier MaeMolSupplier.cpp MaeWriter.cpp)
    include_directories(${maeparser_INCLUDE_DIRS})
endif()

if(RDK_BUILD_COMPRESSED_SUPPLIERS)
    set(compressedsupplier CompressedSDMolSupplier.cpp)
endif()
set(rdmolfiles_sources rdmolfiles.cpp
                       ForwardSDMolSupplier.cpp ${maesupplier}
                       SDMolSupplier.cpp TDTMolSupplier.cpp ${compressedsupplier}
                       SmilesMolSupplier.cpp SmilesWriter.cpp SDWriter.cpp
                       TDTWriter.cpp
                       PDBWriter.cpp
                       MultithreadedSmilesMolSupplier.cpp
                       MultithreadedSDMolSupplier.cpp)

rdkit_python_extension(rdmolfiles
                       ${rdmolfiles_sources}
                       DEST Chem
                       LINK_LIBRARIES  SubstructMatch
                       SmilesParse FileParsers GraphMol
                       MarvinParser
                       RDGeometryLib RDGeneral RDBoost)

rdkit_python_extension(rdtrajectory
                       Trajectory.cpp rdTrajectory.cpp
                       DEST Chem
                       LINK_LIBRARIES Trajectory RDGeometryLib RDGeneral RDBoost GraphMol)

if(RDK_BUILD_COMPRESSED_SUPPLIERS)
  set_target_properties(rdmolfiles PROPERTIES DEFINE_SYMBOL SUPPORT_COMPRESSED_SUPPLIERS)
endif(RDK_BUILD_COMPRESSED_SUPPLIERS)

add_pytest(pyGraphMolWrap
         ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)

add_pytest(pyTestConformerWrap
         ${CMAKE_CURRENT_SOURCE_DIR}/testConformer.py)

add_pytest(pyTestTrajectory
         ${CMAKE_CURRENT_SOURCE_DIR}/testTrajectory.py)

add_pytest(pyTestSGroups
         ${CMAKE_CURRENT_SOURCE_DIR}/testSGroups.py)

add_pytest(pyTestSCSR
         ${CMAKE_CURRENT_SOURCE_DIR}/testSCSR.py)

add_pytest(pySubstructureMatch
         ${CMAKE_CURRENT_SOURCE_DIR}/testSubstructureMatch.py)



if (RDK_TEST_MULTITHREADED)
  add_pytest(pyTestThreads
           ${CMAKE_CURRENT_SOURCE_DIR}/testThreads.py)
  add_pytest(pyTestMultithreadedMolSupplier
         ${CMAKE_CURRENT_SOURCE_DIR}/testMultithreadedMolSupplier.py)
endif (RDK_TEST_MULTITHREADED)

add_pytest(pyTestPropertyLists
         ${CMAKE_CURRENT_SOURCE_DIR}/testPropertyLists.py)

add_pytest(pyCDXMLTest
         ${CMAKE_CURRENT_SOURCE_DIR}/test_cdxml.py)

add_pytest(pySubsetTest
         ${CMAKE_CURRENT_SOURCE_DIR}/test_subset.py)

add_pytest(pyTestLazyNumpy
         ${CMAKE_CURRENT_SOURCE_DIR}/test_lazy_numpy.py)
