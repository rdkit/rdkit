#rdkit_library(RDKit_Minimal
#              minilib.cpp
#              LINK_LIBRARIES SubstructMatch FileParsers SmilesParse GraphMol RDGeometryLib RDGeneral
#                 )

add_executable(RDKit_minimal jswrapper.cpp minilib.cpp)
target_link_libraries(RDKit_minimal SubstructMatch_static FileParsers_static 
    SmilesParse_static GraphMol_static RDGeometryLib_static RDGeneral_static)

set(CMAKE_)
set_target_properties(RDKit_minimal PROPERTIES LINK_FLAGS "--bind")

