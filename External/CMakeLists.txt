add_subdirectory(INCHI-API)
add_subdirectory(AvalonTools)
add_subdirectory(FreeSASA)
add_subdirectory(CoordGen)
add_subdirectory(YAeHMOP)
if(RDK_USE_URF)
#  add_subdirectory(RingDecomposerLib)
include(ExternalProject)

ExternalProject_Add(RingDecomposerLib_project
  GIT_REPOSITORY https://github.com/greglandrum/RingDecomposerLib.git
  GIT_TAG master
  PATCH_COMMAND ""
  BUILD_ALWAYS 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/RingDecomposerLib"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
  TEST_COMMAND "")
endif(RDK_USE_URF)
