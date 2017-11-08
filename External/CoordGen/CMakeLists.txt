add_custom_target(coordgen_support ALL)

if(RDK_BUILD_COORDGEN_SUPPORT)

    file(GLOB MAESOURCES "maeparser/*.cpp")
    rdkit_library(maeparser ${MAESOURCES} SHARED )
    install(TARGETS maeparser DESTINATION ${RDKit_LibDir})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/maeparser)

    file(GLOB CGSOURCES "coordgenlibs/*.cpp")
    rdkit_library(coordgenlib ${CGSOURCES} SHARED LINK_LIBRARIES maeparser)
    install(TARGETS coordgenlib DESTINATION ${RDKit_LibDir})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/coordgenlibs)

  #rdkit_library(RDCoordGenLib coordgen.cpp SHARED LINK_LIBRARIES ${INCHI_LIBRARIES}
  #              GraphMol RDGeneral Depictor SubstructMatch SmilesParse
  #              ${RDKit_THREAD_LIBS})
  rdkit_headers(CoordGen.h DEST GraphMol)
  #if(RDK_BUILD_PYTHON_WRAPPERS)
  #  add_subdirectory(Wrap)
  #endif(RDK_BUILD_PYTHON_WRAPPERS)

  #add_definitions(-DBUILD_INCHI_SUPPORT)

  rdkit_test(testCoordGen test.cpp
    LINK_LIBRARIES
    ${INCHI_LIBRARIES}
    coordgenlib maeparser
    FileParsers SmilesParse SubstructMatch GraphMol
    RDGeneral DataStructs RDGeneral RDGeometryLib
    ${RDKit_THREAD_LIBS})
endif(RDK_BUILD_COORDGEN_SUPPORT)

#if(RDK_BUILD_PYTHON_WRAPPERS)
#  # we need an inchi.py... copy in the appropriate one
#  if(RDK_BUILD_INCHI_SUPPORT)
#   add_custom_command(TARGET inchi_support
#                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/inchi.py
#                     ${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py)
#  else(RDK_BUILD_INCHI_SUPPORT)
#    add_custom_command(TARGET inchi_support
#                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/python/noinchi.py
#                     ${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py)
#  endif(RDK_BUILD_INCHI_SUPPORT)
#  get_directory_property(extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
#  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
#          "${extra_clean_files};${CMAKE_SOURCE_DIR}/rdkit/Chem/inchi.py")
#endif(RDK_BUILD_PYTHON_WRAPPERS)
