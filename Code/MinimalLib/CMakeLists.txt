include_directories(${RDKit_ExternalDir})
include_directories(${RDKit_ExternalDir}/rapidjson-1.1.0/include)

if(RDK_BUILD_MINIMAL_LIB)
    if(RDK_BUILD_FREETYPE_SUPPORT)
        if( ${CMAKE_SYSTEM_NAME} MATCHES "Emscripten")
            set(USE_FLAGS "-s USE_FREETYPE=1")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USE_FLAGS}")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${USE_FLAGS}")
        endif()
    endif()
    add_executable(RDKit_minimal jswrapper.cpp minilib.cpp)
    target_link_libraries(RDKit_minimal MolInterchange_static Abbreviations_static
        CIPLabeler_static MolDraw2D_static Depictor_static RDInchiLib_static
        SubstructMatch_static SubstructLibrary_static FileParsers_static
        SmilesParse_static GraphMol_static RDGeometryLib_static RDGeneral_static)

    set_target_properties(RDKit_minimal PROPERTIES LINK_FLAGS "--bind")
endif(RDK_BUILD_MINIMAL_LIB)

if(RDK_BUILD_CFFI_LIB)

    set(CMAKE_C_STANDARD 99)
    if((NOT WIN32) OR (WIN32 AND RDK_INSTALL_DLLS_MSVC))
        set(LIBS_TO_USE 
            MolStandardize_static DistGeomHelpers_static ForceFieldHelpers_static DistGeometry_static
            ForceField_static Alignment_static
            MolInterchange_static Abbreviations_static CIPLabeler_static 
            MolDraw2D_static Depictor_static 
            RDInchiLib_static SubstructMatch_static FileParsers_static 
            SmilesParse_static GraphMol_static RingDecomposerLib_static RDGeometryLib_static RDGeneral_static )
        if(RDK_URF_LIBS)
            list(APPEND LIBS_TO_USE RingDecomposerLib_static)
        endif()
    else()
        set(LIBS_TO_USE 
        MolStandardize DistGeomHelpers ForceFieldHelpers DistGeometry
        ForceField Alignment
        MolInterchange Abbreviations CIPLabeler 
        MolDraw2D Depictor 
        RDInchiLib SubstructMatch FileParsers 
        SmilesParse GraphMol RingDecomposerLib RDGeometryLib RDGeneral )
        if(RDK_URF_LIBS)
            list(APPEND LIBS_TO_USE RingDecomposerLib)
        endif()
    endif()
    
    add_library(rdkitcffi SHARED cffiwrapper.cpp)
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
#if(NOT MSVC)
    # doesn't work as a test on windows because the DLL needs to be either in the PATH OR
	# in the same dir as the executable
    add_test(cffi_test ${EXECUTABLE_OUTPUT_PATH}/cffi_test)
#endif
endif(RDK_BUILD_CFFI_LIB)
