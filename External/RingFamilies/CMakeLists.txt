if(RDK_USE_URF)

include(ExternalProject)

ExternalProject_Add(RingDecomposerLib_project
    GIT_REPOSITORY https://github.com/greglandrum/RingDecomposerLib.git
    GIT_TAG master
    PATCH_COMMAND ""
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/RingDecomposerLib"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
    TEST_COMMAND "")

set(RDK_URF_LIBS RingDecomposerLibStatic
    CACHE STRING "the external libraries" FORCE)

rdkit_library(URFLib URF.cpp SHARED LINK_LIBRARIES ${RDK_URF_LIBS} )
add_dependencies(URFLib RingDecomposerLib_project)

endif(RDK_USE_URF)
