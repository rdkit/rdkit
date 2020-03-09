include_directories(${RDKit_ExternalDir})
include_directories(${RDKit_ExternalDir}/rapidjson-1.1.0/include)

add_executable(RDKit_minimal jswrapper.cpp minilib.cpp)
target_link_libraries(RDKit_minimal MolDraw2D_static Depictor_static RDInchiLib_static SubstructMatch_static FileParsers_static 
    SmilesParse_static GraphMol_static RDGeometryLib_static RDGeneral_static)

set_target_properties(RDKit_minimal PROPERTIES LINK_FLAGS "--bind")

