remove_definitions(-DRDKIT_GRAPHMOL_BUILD)
rdkit_python_extension(rdchem
                       rdchem.cpp Table.cpp Atom.cpp Bond.cpp Mol.cpp
                       Conformer.cpp RingInfo.cpp EditableMol.cpp MonomerInfo.cpp
                       ResonanceMolSupplier.cpp MolBundle.cpp StereoGroup.cpp Sgroup.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       SmilesParse ChemTransforms SubstructMatch GraphMol RDGeometryLib RDGeneral RDBoost)

rdkit_python_extension(rdmolops
                       rdmolops.cpp MolOps.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       ChemReactions MolDraw2D Depictor
                       FileParsers SubstructMatch Fingerprints ChemTransforms
                       Catalogs Subgraphs SmilesParse MolTransforms GraphMol EigenSolvers
RDGeometryLib DataStructs RDGeneral RDBoost  )

rdkit_python_extension(rdqueries
                       rdqueries.cpp Queries.cpp
                       DEST Chem
                       LINK_LIBRARIES
                       GraphMol RDGeometryLib DataStructs RDGeneral RDBoost  )

if(RDK_BUILD_COORDGEN_SUPPORT)
    set (maesupplier MaeMolSupplier.cpp)
    include_directories(${maeparser_INCLUDE_DIRS})
endif()

if(RDK_BUILD_COMPRESSED_SUPPLIERS)
find_package(Boost 1.56.0 COMPONENTS iostreams REQUIRED)
set(rdmolfiles_sources rdmolfiles.cpp
                       ForwardSDMolSupplier.cpp ${maesupplier}
                       SDMolSupplier.cpp TDTMolSupplier.cpp
                       SmilesMolSupplier.cpp SmilesWriter.cpp SDWriter.cpp
                       TDTWriter.cpp CompressedSDMolSupplier.cpp
                       PDBWriter.cpp )
else(RDK_BUILD_COMPRESSED_SUPPLIERS)
set(rdmolfiles_sources rdmolfiles.cpp
                       ForwardSDMolSupplier.cpp ${maesupplier}
                       SDMolSupplier.cpp TDTMolSupplier.cpp
                       SmilesMolSupplier.cpp SmilesWriter.cpp SDWriter.cpp
                       TDTWriter.cpp
                       PDBWriter.cpp )
endif(RDK_BUILD_COMPRESSED_SUPPLIERS)

rdkit_python_extension(rdmolfiles
                       ${rdmolfiles_sources}
                       DEST Chem
                       LINK_LIBRARIES  SubstructMatch
		       SmilesParse FileParsers GraphMol
                       RDGeometryLib RDGeneral RDBoost)

rdkit_python_extension(rdtrajectory
                       Trajectory.cpp rdTrajectory.cpp
                       DEST Chem
                       LINK_LIBRARIES Trajectory RDGeometryLib RDGeneral GraphMol)

if(RDK_BUILD_COMPRESSED_SUPPLIERS)
  set_target_properties(rdmolfiles PROPERTIES DEFINE_SYMBOL SUPPORT_COMPRESSED_SUPPLIERS )
endif(RDK_BUILD_COMPRESSED_SUPPLIERS)

add_pytest(pyGraphMolWrap
         ${CMAKE_CURRENT_SOURCE_DIR}/rough_test.py)

add_pytest(pyTestConformerWrap
         ${CMAKE_CURRENT_SOURCE_DIR}/testConformer.py)

add_pytest(pyTestTrajectory
         ${CMAKE_CURRENT_SOURCE_DIR}/testTrajectory.py)

add_pytest(pyTestSGroups
         ${CMAKE_CURRENT_SOURCE_DIR}/testSGroups.py)

if (RDK_TEST_MULTITHREADED)
  add_pytest(pyTestThreads
           ${CMAKE_CURRENT_SOURCE_DIR}/testThreads.py)
endif (RDK_TEST_MULTITHREADED)

add_pytest(pyTestPropertyLists
         ${CMAKE_CURRENT_SOURCE_DIR}/testPropertyLists.py)
