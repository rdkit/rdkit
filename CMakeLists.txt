cmake_minimum_required(VERSION 2.6)
project (RDKit)

enable_testing()

add_subdirectory(External)
add_subdirectory(Code)
add_subdirectory(Projects)
add_subdirectory(rdkit)
