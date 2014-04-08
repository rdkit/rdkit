cmake_minimum_required (VERSION 2.8)

set(RDKit_Year "2014")
set(RDKit_Month "03")
set(RDKit_ABI "1")

set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package (RDKit REQUIRED)

set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../Code/cmake/Modules)
set(RDKit_LibDir ${PROJECT_SOURCE_DIR}/../../lib)
include_directories (${RDKIT_INCLUDE_DIR})

include(RDKitUtils)
install(EXPORT ${RDKit_EXPORTED_TARGETS} DESTINATION ${RDKit_LibDir})

set(RDK_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(RDK_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(RDK_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(RDK_PYTHON_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/rdkit")

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

rdkit_library(ConformerParser ConformerParser.cpp
              LINK_LIBRARIES ${RDKIT_LIBARIES})

set(RDK_BUILD_CPP_TESTS ON)              

enable_testing()

rdkit_test(testConformerParser test.cpp 
           LINK_LIBRARIES ConformerParser
           ${RDKIT_LIBRARIES} 
)

add_subdirectory(Wrap)