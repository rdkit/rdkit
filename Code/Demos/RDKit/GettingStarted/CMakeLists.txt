cmake_minimum_required (VERSION 2.6)
project (GettingStartedDemo)

find_package (RDKit 1.2011.09 REQUIRED)

include_directories (${RDKit_INCLUDE_DIRS})

add_executable (sample sample.cpp)
target_link_libraries (sample ChemReactions FileParsers
   SmilesParse Depictor SubstructMatch GraphMol
   DataStructs RDGeometryLib RDGeneral
                      )
