rdkit_library(GraphMol
        Atom.cpp QueryAtom.cpp QueryBond.cpp Bond.cpp
        MolOps.cpp FindRings.cpp ROMol.cpp RWMol.cpp PeriodicTable.cpp
        atomic_data.cpp QueryOps.cpp MolPickler.cpp Canon.cpp
        AtomIterators.cpp BondIterators.cpp Aromaticity.cpp Kekulize.cpp
        ConjugHybrid.cpp AddHs.cpp
        Matrices.cpp Chirality.cpp RingInfo.cpp Conformer.cpp
        Renumber.cpp AdjustQuery.cpp Resonance.cpp StereoGroup.cpp
        new_canon.cpp SubstanceGroup.cpp FindStereo.cpp MonomerInfo.cpp
              NontetrahedralStereo.cpp Atropisomers.cpp
              WedgeBonds.cpp
        SHARED
              LINK_LIBRARIES RDGeometryLib RDGeneral )
target_compile_definitions(GraphMol PRIVATE RDKIT_GRAPHMOL_BUILD)
if (RDK_USE_URF)
    target_link_libraries(GraphMol PUBLIC ${RDK_URF_LIBS})
    if (NOT WIN32 AND RDK_INSTALL_STATIC_LIBS)
        target_link_libraries(GraphMol_static PUBLIC ${RDK_URF_LIBS}_static)
    endif ()
endif ()

rdkit_headers(Atom.h
        atomic_data.h
        AtomIterators.h
        Atropisomers.h
        Bond.h
        BondIterators.h
        Canon.h
        Chirality.h
        Conformer.h
        details.h
        GraphMol.h
        MolOps.h
        MolPickler.h
        PeriodicTable.h
        QueryAtom.h
        QueryBond.h
        QueryOps.h
        RDKitBase.h
        RDKitQueries.h
        Resonance.h
        RingInfo.h
        Rings.h
        ROMol.h
        RWMol.h
        SanitException.h
        SubstanceGroup.h
        StereoGroup.h
        MonomerInfo.h
        new_canon.h
        MolBundle.h
        DEST GraphMol)

add_subdirectory(SmilesParse)
add_subdirectory(Depictor)
add_subdirectory(MarvinParse)
add_subdirectory(FileParsers)
add_subdirectory(Substruct)
add_subdirectory(GenericGroups)
add_subdirectory(ChemReactions)
add_subdirectory(ChemTransforms)
add_subdirectory(TautomerQuery)

add_subdirectory(Subgraphs)
add_subdirectory(FilterCatalog)
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

add_subdirectory(MolDraw2D)

add_subdirectory(FMCS)
add_subdirectory(MolHash)
add_subdirectory(MMPA)
add_subdirectory(RascalMCES)

add_subdirectory(CIPLabeler)
add_subdirectory(Deprotect)

if (RDK_BUILD_STRUCTCHECKER_SUPPORT)
    add_subdirectory(StructChecker)
endif ()

add_subdirectory(ReducedGraphs)
add_subdirectory(Trajectory)

add_subdirectory(SubstructLibrary)
add_subdirectory(RGroupDecomposition)

if (RDK_BUILD_XYZ2MOL_SUPPORT)
    add_subdirectory(DetermineBonds)
endif ()

if (RDK_BUILD_MOLINTERCHANGE_SUPPORT)
    add_subdirectory(MolInterchange)
endif (RDK_BUILD_MOLINTERCHANGE_SUPPORT)

if (RDK_BUILD_SLN_SUPPORT AND NOT CMAKE_COMPILER_IS_CLANG)
    if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 4.8)
        message("Skipping SLN build with gcc version <= 4.8")
    else ()
        add_subdirectory(SLNParse)
    endif ()
endif ()

if (RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif ()

add_subdirectory(MolStandardize)
add_subdirectory(ScaffoldNetwork)
add_subdirectory(MolEnumerator)
add_subdirectory(Abbreviations)
add_subdirectory(GeneralizedSubstruct)


rdkit_test(graphmolTest1 test1.cpp LINK_LIBRARIES FileParsers SmilesParse GraphMol
        SubstructMatch ForceFieldHelpers ForceField)

rdkit_test(graphmolcpTest cptest.cpp LINK_LIBRARIES SmilesParse GraphMol)

rdkit_test(graphmolqueryTest querytest.cpp LINK_LIBRARIES SubstructMatch SmilesParse GraphMol)

rdkit_test(graphmolMolOpsTest molopstest.cpp
        LINK_LIBRARIES SubstructMatch FileParsers
        SmilesParse GraphMol)

rdkit_test(graphmoltestChirality testChirality.cpp
        LINK_LIBRARIES SmilesParse FileParsers
        GraphMol)

rdkit_test(graphmoltestPickler testPickler.cpp
        LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol)

rdkit_test(graphmoltestPicklerGlobalSetting testPicklerGlobalSettings.cpp
        LINK_LIBRARIES SmilesParse SubstructMatch FileParsers GraphMol)

rdkit_test(graphmolIterTest itertest.cpp LINK_LIBRARIES SmilesParse GraphMol)

rdkit_test(hanoiTest hanoitest.cpp LINK_LIBRARIES
        SubstructMatch SmilesParse FileParsers GraphMol
        )

rdkit_test(graphmolMemTest1 memtest1.cpp LINK_LIBRARIES SmilesParse GraphMol)

rdkit_test(resMolSupplierTest resMolSupplierTest.cpp
        LINK_LIBRARIES SmilesParse GraphMol SubstructMatch FileParsers)

rdkit_test(molBundleTest testMolBundle.cpp
        LINK_LIBRARIES SmilesParse GraphMol SubstructMatch FileParsers)

rdkit_test(testSGroup testSGroup.cpp LINK_LIBRARIES FileParsers GraphMol)

rdkit_test(test-valgrind test-valgrind.cpp LINK_LIBRARIES SmilesParse FileParsers GraphMol)

rdkit_catch_test(graphmolTestsCatch catch_graphmol.cpp
        LINK_LIBRARIES SubstructMatch FileParsers SmilesParse GraphMol)

rdkit_catch_test(graphmolOrganometallicsCatch catch_organometallics.cpp
        LINK_LIBRARIES FileParsers SmilesParse GraphMol)

rdkit_catch_test(graphmolSGroupCatch catch_sgroups.cpp
        LINK_LIBRARIES SmilesParse FileParsers GraphMol)

rdkit_catch_test(graphmolAdjustQueryCatch catch_adjustquery.cpp
        LINK_LIBRARIES SubstructMatch FileParsers SmilesParse GraphMol)

rdkit_catch_test(chiralityTestsCatch catch_chirality.cpp
        LINK_LIBRARIES CIPLabeler FileParsers SmilesParse GraphMol)

rdkit_catch_test(canonTestsCatch catch_canon.cpp
        LINK_LIBRARIES FileParsers SmilesParse GraphMol)

rdkit_catch_test(moliteratorTestsCatch catch_moliterators.cpp
        LINK_LIBRARIES SubstructMatch SmilesParse GraphMol)

rdkit_catch_test(queryTestsCatch catch_queries.cpp
        LINK_LIBRARIES SmilesParse GraphMol)

rdkit_catch_test(molbundleTestsCatch catch_molbundle.cpp
        LINK_LIBRARIES SmilesParse GraphMol)

rdkit_catch_test(pickleTestsCatch catch_pickles.cpp
        LINK_LIBRARIES FileParsers SmilesParse GraphMol)

rdkit_catch_test(tableTestsCatch catch_periodictable.cpp
        LINK_LIBRARIES GraphMol)

rdkit_catch_test(molopsTestsCatch catch_molops.cpp
        LINK_LIBRARIES SmilesParse GraphMol)

