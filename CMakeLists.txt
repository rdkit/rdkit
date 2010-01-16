cmake_minimum_required(VERSION 2.6)
project (RDKit)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)


enable_testing()

add_subdirectory(External)
add_subdirectory(Code)
add_subdirectory(Projects)
add_subdirectory(rdkit)
