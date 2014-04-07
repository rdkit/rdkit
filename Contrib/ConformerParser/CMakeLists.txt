cmake_minimum_required (VERSION 2.8) 
set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package (RDKit REQUIRED)

set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../Code/cmake/Modules)
set(RDKit_LibDir ${PROJECT_SOURCE_DIR}/../../lib)
include_directories (${RDKIT_INCLUDE_DIR})
include(RDKitUtils)

macro(local_library)
  PARSE_ARGUMENTS(RDKLIB
    "LINK_LIBRARIES;DEST"
    "SHARED"
    ${ARGN})
  CAR(RDKLIB_NAME ${RDKLIB_DEFAULT_ARGS})
  CDR(RDKLIB_SOURCES ${RDKLIB_DEFAULT_ARGS})
  if(MSVC)
    add_library(${RDKLIB_NAME} ${RDKLIB_SOURCES})
    target_link_libraries(${RDKLIB_NAME} ${Boost_SYSTEM_LIBRARY} )
    INSTALL(TARGETS ${RDKLIB_NAME} EXPORT ${RDKit_EXPORTED_TARGETS}
            DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
            COMPONENT dev )
  else(MSVC)
    # we're going to always build in shared mode since we
    # need exceptions to be (correctly) catchable across
    # boundaries. As of now (June 2010), this doesn't work
    # with g++ unless libraries are shared.
      add_library(${RDKLIB_NAME} SHARED ${RDKLIB_SOURCES})
      INSTALL(TARGETS ${RDKLIB_NAME} EXPORT ${RDKit_EXPORTED_TARGETS}
              DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
              COMPONENT runtime )
      if(RDK_INSTALL_STATIC_LIBS)
        add_library(${RDKLIB_NAME}_static ${RDKLIB_SOURCES})
        INSTALL(TARGETS ${RDKLIB_NAME}_static EXPORT ${RDKit_EXPORTED_TARGETS}
                DESTINATION ${RDKit_LibDir}/${RDKLIB_DEST}
                COMPONENT dev )
      endif(RDK_INSTALL_STATIC_LIBS)
    IF(RDKLIB_LINK_LIBRARIES)
      target_link_libraries(${RDKLIB_NAME} ${RDKLIB_LINK_LIBRARIES})
    ENDIF(RDKLIB_LINK_LIBRARIES)
  endif(MSVC)
  if(WIN32)
    set_target_properties(${RDKLIB_NAME} PROPERTIES 
                          OUTPUT_NAME "${RDKLIB_NAME}" 
                          VERSION "${RDKit_ABI}.${RDKit_Year}.${RDKit_Month}")
  else(WIN32)
    set_target_properties(${RDKLIB_NAME} PROPERTIES 
                          OUTPUT_NAME ${RDKLIB_NAME} 
                          )
  endif(WIN32)			  

endmacro(local_library)

find_package(PythonLibs)
if(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
  include_directories(${PYTHON_INCLUDE_PATH})
else(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
  include_directories(${PYTHON_INCLUDE_DIR})
endif(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_NUMPY_INCLUDE_PATH})
find_package(Boost 1.42.0 COMPONENTS python REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

local_library(ConformerParser ConformerParser.cpp
              LINK_LIBRARIES ${RDKIT_LIBARIES})

rdkit_test(testConformerParser test.cpp 
           LINK_LIBRARIES ConformerParser
           ${RDKIT_LIBRARIES} 
)

add_subdirectory(Wrap)