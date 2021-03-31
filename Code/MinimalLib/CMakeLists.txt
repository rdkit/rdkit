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
    target_link_libraries(RDKit_minimal MolInterchange_static Abbreviations_static CIPLabeler_static MolDraw2D_static Depictor_static RDInchiLib_static SubstructMatch_static FileParsers_static 
        SmilesParse_static GraphMol_static RDGeometryLib_static RDGeneral_static)

    set_target_properties(RDKit_minimal PROPERTIES LINK_FLAGS "--bind")
endif(RDK_BUILD_MINIMAL_LIB)

if(RDK_BUILD_CFFI_LIB)
    set(CMAKE_C_STANDARD 99)
    #rdkit _ library(rdkitcffi cffiwrapper.cpp 
    #  SHARED
    #  LINK_LIBRARIES MolInterchange_static Abbreviations_static CIPLabeler_static MolDraw2D_static Depictor_static 
    #    RDInchiLib_static SubstructMatch_static FileParsers_static 
    #    SmilesParse_static GraphMol_static RDGeometryLib_static RDGeneral_static )

    set(LIBS_TO_USE MolInterchange_static Abbreviations_static CIPLabeler_static MolDraw2D_static Depictor_static 
        RDInchiLib_static SubstructMatch_static FileParsers_static 
        SmilesParse_static GraphMol_static RDGeometryLib_static RDGeneral_static )
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


    add_executable(cffi_test cffi_test.c)
    target_link_libraries(cffi_test rdkitcffi)
    add_test(cffi_test ${EXECUTABLE_OUTPUT_PATH}/cffi_test)

endif(RDK_BUILD_CFFI_LIB)
