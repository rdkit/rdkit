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
  set(PUBCHEM_COMMIT_SHA 169a241)
  if(NOT DEFINED PUBCHEMSHAPE_URL)
    set(PUBCHEMSHAPE_URL "https://github.com/ncbi/pubchem-align3d/archive/${PUBCHEM_COMMIT_SHA}.tar.gz")
  endif()
  if(NOT DEFINED PUBCHEMSHAPE_MD5SUM)
    set(PUBCHEMSHAPE_MD5SUM "97ae9388968b2c4983b94632aa668d95")
  endif()
  if(NOT DEFINED PUBCHEMSHAPE_BASE)
    string(REGEX REPLACE "^.*/" "" PUBCHEMSHAPE_BASE "${PUBCHEMSHAPE_URL}")
  endif()
  downloadAndCheckMD5(${PUBCHEMSHAPE_URL} "${CMAKE_CURRENT_SOURCE_DIR}/${PUBCHEMSHAPE_BASE}" ${PUBCHEMSHAPE_MD5SUM})
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf  
    ${CMAKE_CURRENT_SOURCE_DIR}/${PUBCHEM_COMMIT_SHA}.tar.gz 
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} -E rename "${CMAKE_CURRENT_SOURCE_DIR}/pubchem-align3d-${PUBCHEM_COMMIT_SHA}*" 
       ${PUBCHEMSHAPE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

rdkit_library(pubchem_align3d ./pubchem-align3d/shape_functions1.cpp 
            ./pubchem-align3d/shape_functions2.cpp ./pubchem-align3d/shape_neighbor.cpp)

rdkit_library(PubChemShape PubChemShape.cpp SHARED
 LINK_LIBRARIES pubchem_align3d SmilesParse SubstructMatch)
target_compile_definitions(PubChemShape PRIVATE RDKIT_PUBCHEMSHAPE_BUILD)

rdkit_headers(PubChemShape.h DEST GraphMol)

rdkit_catch_test(shape_test test.cpp LINK_LIBRARIES PubChemShape FileParsers)

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
