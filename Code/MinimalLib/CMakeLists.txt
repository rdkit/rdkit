include_directories(${RDKit_ExternalDir})
include_directories(${RDKit_ExternalDir}/rapidjson-1.1.0/include)

if(RDK_BUILD_FREETYPE_SUPPORT)
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Emscripten")
        set(USE_FLAGS "-s USE_FREETYPE=1")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USE_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${USE_FLAGS}")
    endif()
endif()
add_executable(RDKit_minimal jswrapper.cpp minilib.cpp)
target_link_libraries(RDKit_minimal Abbreviations_static CIPLabeler_static MolDraw2D_static Depictor_static RDInchiLib_static SubstructMatch_static FileParsers_static 
    SmilesParse_static GraphMol_static RDGeometryLib_static RDGeneral_static)

set_target_properties(RDKit_minimal PROPERTIES LINK_FLAGS "--bind")

