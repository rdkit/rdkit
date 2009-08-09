cmake_minimum_required(VERSION 2.6)
project (RDKit)

include(BoostUtils.cmake)
macro(rdkit_python_extension)
  PARSE_ARGUMENTS(RDKPY
    "LINK_LIBRARIES;DEPENDS;DEST"
    ""
    ${ARGN}
    )
  CAR(RDKPY_NAME ${RDKPY_DEFAULT_ARGS})
  CDR(RDKPY_SOURCES ${RDKPY_DEFAULT_ARGS})
  PYTHON_ADD_MODULE(${RDKPY_NAME} ${RDKPY_SOURCES}  )
  set_target_properties(${RDKPY_NAME} PROPERTIES PREFIX "")
  
  target_link_libraries(${RDKPY_NAME} ${RDKPY_LINK_LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

  INSTALL(TARGETS ${RDKPY_NAME} 
      LIBRARY DESTINATION ${RDKit_PythonDir}/${RDKPY_DEST} )
endmacro(rdkit_python_extension)




# pull in python:
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
link_directories(${PYTHON_LIBRARIES})

# pull in boost:
find_package(Boost COMPONENTS python thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

enable_testing()

# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing")
endif()

set(RDKit_PythonDir $ENV{RDBASE}/rdkit)

include_directories(.)
add_subdirectory(RDGeneral)
add_subdirectory(RDBoost)
add_subdirectory(DataStructs)
add_subdirectory(Geometry)
add_subdirectory(GraphMol)

