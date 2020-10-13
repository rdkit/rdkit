cmake_minimum_required (VERSION 3.5)
project (GettingStartedDemo)

find_package (RDKit 1.2011.09 REQUIRED)

add_executable (sample sample.cpp)
target_link_libraries (sample ChemReactions FileParsers SmilesParse
                       Depictor SubstructMatch GraphMol RDGeneral)
