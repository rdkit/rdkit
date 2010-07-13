cmake_minimum_required(VERSION 2.6)
project (RDKit)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

enable_testing()

set(CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Code/cmake/Modules/")

option(RDK_BUILD_SWIG_WRAPPERS "build the experimental SWIG wrappers" OFF )
option(RDK_BUILD_PYTHON_WRAPPERS "build the standard python wrappers" ON )
option(RDK_BUILD_COMPRESSED_SUPPLIERS "build in support for compressed MolSuppliers" OFF )

#------- 
# Config variables:
set(RDKit_CodeDir "${CMAKE_SOURCE_DIR}/Code")
set(RDKit_ExternalDir "${CMAKE_SOURCE_DIR}/External")
set(RDKit_BinDir "bin")
set(RDKit_LibDir "lib${LIB_SUFFIX}")
set(RDKit_ShareDir "share/RDKit")

# defines macros: rdkit_python_extension, rdkit_test
include(RDKitUtils)

# extra boost versions
if(MSVC)
SET(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0" "1.41" "1.41.0" "1.40" "1.40.0")
endif(MSVC)

if(RDK_BUILD_PYTHON_WRAPPERS)
  #------- 
  # pull in python:
  find_package(PythonLibs)
  if(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
    include_directories(${PYTHON_INCLUDE_PATH})
  else(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
    include_directories(${PYTHON_INCLUDE_DIR})
  endif(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
  link_directories(${PYTHON_LIBRARIES})
  find_package(NumPy REQUIRED)
  include_directories(${PYTHON_NUMPY_INCLUDE_PATH})
  find_package(Boost 1.39.0 COMPONENTS python REQUIRED)

# Determine correct installation directory for Python bindings
# From http://plplot.svn.sourceforge.net/viewvc/plplot/trunk/cmake/modules/python.cmake?revision=11014
  execute_process(
    COMMAND
    ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}')"
    OUTPUT_VARIABLE PYTHON_INSTDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  message("Python Install directory ${PYTHON_INSTDIR}")

  install(DIRECTORY rdkit DESTINATION ${PYTHON_INSTDIR}
    PATTERN ".svn" EXCLUDE )


else(RDK_BUILD_PYTHON_WRAPPERS)
  find_package(Boost 1.39.0 REQUIRED)
endif(RDK_BUILD_PYTHON_WRAPPERS)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${RDKit_CodeDir})

# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing -fPIC")
endif()

install(DIRECTORY Data DESTINATION ${RDKit_ShareDir}/Data
    PATTERN ".svn" EXCLUDE)
install(DIRECTORY Docs DESTINATION ${RDKit_ShareDir}/Docs
    PATTERN ".svn" EXCLUDE)
install(DIRECTORY Projects DESTINATION ${RDKit_ShareDir}/Projects
    PATTERN ".svn" EXCLUDE)
install(DIRECTORY Contrib DESTINATION ${RDKit_ShareDir}/Contrib
    PATTERN ".svn" EXCLUDE)
install(FILES README license.txt DESTINATION ${RDKit_ShareDir})


add_subdirectory(External)
add_subdirectory(Code)
add_subdirectory(Projects)
add_subdirectory(rdkit)
