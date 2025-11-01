include_directories(${RDKit_ExternalDir})
include_directories(${RDKit_ExternalDir}/rapidjson/rapidjson-1.1.0/include)

if(RDK_BUILD_MINIMAL_LIB)
    set(MINIMAL_LIB_LIBRARIES "MolInterchange;Abbreviations;"
        "CIPLabeler;MolDraw2D;Depictor;"
        "Descriptors;SubstructMatch;FileParsers;"
        "SmilesParse;GraphMol;RDGeometryLib;"
        "RDGeneral;RGroupDecomposition;Fingerprints")
    if(RDK_BUILD_INCHI_SUPPORT)
        add_definitions(-DRDK_BUILD_INCHI_SUPPORT)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};RDInchiLib")
    endif()
    if(RDK_BUILD_MINIMAL_LIB_RXN)
        add_definitions(-DRDK_BUILD_MINIMAL_LIB_RXN)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};ChemReactions")
    endif()
    if(RDK_BUILD_MINIMAL_LIB_SUBSTRUCTLIBRARY)
        add_definitions(-DRDK_BUILD_MINIMAL_LIB_SUBSTRUCTLIBRARY)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};SubstructLibrary")
    endif()
    if(RDK_BUILD_MINIMAL_LIB_MCS)
        add_definitions(-DRDK_BUILD_MINIMAL_LIB_MCS)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};FMCS")
    endif()
    if(RDK_BUILD_MINIMAL_LIB_MMPA)
        add_definitions(-DRDK_BUILD_MINIMAL_LIB_MMPA)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};MMPA")
    endif()
    if(RDK_BUILD_MINIMAL_LIB_MOLZIP)
        add_definitions(-DRDK_BUILD_MINIMAL_LIB_MOLZIP)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};ChemTransforms")
    endif()
    if(RDK_BUILD_MINIMAL_LIB_RGROUPDECOMP)
        add_definitions(-DRDK_BUILD_MINIMAL_LIB_RGROUPDECOMP)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};RGroupDecomposition")
    endif()
    if(RDK_BUILD_MINIMAL_LIB_SCAFFOLDNETWORK)
        add_definitions(-DRDK_BUILD_MINIMAL_LIB_SCAFFOLDNETWORK)
        set(MINIMAL_LIB_LIBRARIES "${MINIMAL_LIB_LIBRARIES};ScaffoldNetwork")
    endif()
    if(RDK_BUILD_FREETYPE_SUPPORT)
        if( ${CMAKE_SYSTEM_NAME} MATCHES "Emscripten")
            set(USE_FLAGS "-s USE_FREETYPE=1")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USE_FLAGS}")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${USE_FLAGS}")
        endif()
    endif()
    add_executable(RDKit_minimal jswrapper.cpp minilib.cpp JSONParsers.cpp)
    target_link_libraries(RDKit_minimal ${MINIMAL_LIB_LIBRARIES})

    set_target_properties(RDKit_minimal PROPERTIES LINK_FLAGS "--bind")
endif(RDK_BUILD_MINIMAL_LIB)

if(RDK_BUILD_CFFI_LIB)
    set(CMAKE_C_STANDARD 99)

    set(LIBS_TO_USE 
        MolStandardize DistGeomHelpers ForceFieldHelpers DistGeometry
        ForceField Alignment
        MolInterchange Abbreviations CIPLabeler 
        MolDraw2D Depictor Descriptors
        SubstructMatch FileParsers SmilesParse GraphMol
        RDGeometryLib RDGeneral RGroupDecomposition Fingerprints)
    if(RDK_BUILD_INCHI_SUPPORT)
        add_definitions(-DRDK_BUILD_INCHI_SUPPORT)
        list(APPEND LIBS_TO_USE RDInchiLib)
    endif()
    if(RDK_URF_LIBS)
        list(APPEND LIBS_TO_USE RingDecomposerLib)
    endif()
    
    if(RDK_CFFI_STATIC AND ((NOT WIN32) OR (WIN32 AND RDK_INSTALL_DLLS_MSVC)))
        set(staticLibSuffix "_static")
        set(tmpLibs "")
        foreach(lib ${LIBS_TO_USE})
            set(tmpLibs "${tmpLibs}${lib}${staticLibSuffix};")
        endforeach()
        set(LIBS_TO_USE ${tmpLibs})
    endif()
    
    add_library(rdkitcffi SHARED cffiwrapper.cpp JSONParsers.cpp)
    target_link_libraries(rdkitcffi PUBLIC rdkit_base)
    target_link_libraries(rdkitcffi PUBLIC ${LIBS_TO_USE})
    INSTALL(TARGETS rdkitcffi EXPORT rdkit-targets
            DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
            COMPONENT runtime )
    set_target_properties(rdkitcffi PROPERTIES
                            OUTPUT_NAME "rdkitcffi"
                            VERSION "${RDKit_ABI}.${RDKit_Year}.${RDKit_Month}.${RDKit_Revision}"
                            VERSION ${RDKit_VERSION}
                            SOVERSION ${RDKit_ABI} )
    set_target_properties(rdkitcffi PROPERTIES
                            ARCHIVE_OUTPUT_DIRECTORY ${RDK_ARCHIVE_OUTPUT_DIRECTORY}
                            RUNTIME_OUTPUT_DIRECTORY ${RDK_RUNTIME_OUTPUT_DIRECTORY}
                            LIBRARY_OUTPUT_DIRECTORY ${RDK_LIBRARY_OUTPUT_DIRECTORY})
if(MSVC OR WIN32)
  set_target_properties(rdkitcffi PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()


    add_executable(cffi_test cffi_test.c)
    target_link_libraries(cffi_test rdkitcffi)
    set_target_properties(cffi_test PROPERTIES LINKER_LANGUAGE CXX)

#if(NOT MSVC)
    # doesn't work as a test on windows because the DLL needs to be either in the PATH OR
	# in the same dir as the executable
    add_test(cffi_test ${EXECUTABLE_OUTPUT_PATH}/cffi_test)
#endif
endif(RDK_BUILD_CFFI_LIB)
