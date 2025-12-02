if(NOT RDK_BUILD_PUBCHEMSHAPE_SUPPORT)
  return()
endif()

if(NOT DEFINED PUBCHEMSHAPE_DIR)
  set(PUBCHEMSHAPE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/pubchem-align3d")
  set(fileToCheck "${PUBCHEMSHAPE_DIR}/shape_functions1.cpp")
  set(needDownload "TRUE")
  if(EXISTS "${fileToCheck}")
    set(needDownload "FALSE")
  endif()
else()
  set(needDownload "FALSE")
endif()

if(needDownload)
  set(PUBCHEM_COMMIT_SHA daefab3)
  if(NOT DEFINED PUBCHEMSHAPE_URL)
    set(PUBCHEMSHAPE_URL "https://github.com/ncbi/pubchem-align3d/archive/${PUBCHEM_COMMIT_SHA}.tar.gz")
  endif()
  if(NOT DEFINED PUBCHEMSHAPE_MD5SUM)
    set(PUBCHEMSHAPE_MD5SUM "a78a6273851709083ac3803d62cb5ae6")
  endif()
  if(NOT DEFINED PUBCHEMSHAPE_BASE)
    string(REGEX REPLACE "^.*/" "" PUBCHEMSHAPE_BASE "${PUBCHEMSHAPE_URL}")
  endif()
  downloadAndCheckMD5(${PUBCHEMSHAPE_URL} "${CMAKE_CURRENT_SOURCE_DIR}/${PUBCHEMSHAPE_BASE}" ${PUBCHEMSHAPE_MD5SUM})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf
    ${CMAKE_CURRENT_SOURCE_DIR}/${PUBCHEM_COMMIT_SHA}.tar.gz
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  file(GLOB tar_dirname ${CMAKE_CURRENT_SOURCE_DIR}/pubchem-align3d-${PUBCHEM_COMMIT_SHA}*)
  execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${tar_dirname}
       ${PUBCHEMSHAPE_DIR}
       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

endif()

# simple patch for a typo in the pubchem_align library
file(READ ${PUBCHEMSHAPE_DIR}/shape_neighbor.cpp FILE_CONTENTS)
string(REPLACE "memcpy( qom, old_quattrans, 7 * sizeof( float ) );" "memcpy( qom, old_quattrans, 7 * sizeof( double ) );" FILE_CONTENTS "${FILE_CONTENTS}")
file(WRITE ${PUBCHEMSHAPE_DIR}/shape_neighbor.cpp "${FILE_CONTENTS}")

rdkit_library(pubchem_align3d ./pubchem-align3d/shape_functions1.cpp
            ./pubchem-align3d/shape_functions2.cpp ./pubchem-align3d/shape_neighbor.cpp SHARED)
if((MSVC AND RDK_INSTALL_DLLS_MSVC) OR ((NOT MSVC) AND WIN32))
  set_target_properties(pubchem_align3d PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
endif()

rdkit_library(PubChemShape PubChemShape.cpp SHARED
 LINK_LIBRARIES pubchem_align3d SmilesParse SubstructMatch)
target_compile_definitions(PubChemShape PRIVATE RDKIT_PUBCHEMSHAPE_BUILD)

rdkit_headers(PubChemShape.hpp DEST GraphMol)

rdkit_catch_test(shape_test test.cpp LINK_LIBRARIES PubChemShape FileParsers MolAlign MolTransforms)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
