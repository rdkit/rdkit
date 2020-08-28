cmake_minimum_required(VERSION 3.5)

project (RDKit)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Enable CTest
include(CTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Code/cmake/Modules/")

add_library(rdkit_base INTERFACE)
add_dependencies(rdkit_base catch)

option(RDK_BUILD_SWIG_WRAPPERS "build the SWIG wrappers" OFF )
option(RDK_BUILD_PYTHON_WRAPPERS "build the standard python wrappers" ON )
option(RDK_BUILD_COMPRESSED_SUPPLIERS "build in support for compressed MolSuppliers" OFF )
option(RDK_BUILD_ANI "build support for ANI style ForceField" OFF)
option(RDK_BUILD_INCHI_SUPPORT "build the rdkit inchi wrapper" OFF )
option(RDK_BUILD_AVALON_SUPPORT "install support for the avalon toolkit. Use the variable AVALONTOOLS_DIR to set the location of the source." OFF )
option(RDK_BUILD_PGSQL "build the PostgreSQL cartridge" OFF )
option(RDK_BUILD_RPATH_SUPPORT "build shared libraries using rpath" OFF)
option(RDK_PGSQL_STATIC "statically link rdkit libraries into the PostgreSQL cartridge" ON )
option(RDK_BUILD_CONTRIB "build the Contrib directory" OFF )
option(RDK_INSTALL_INTREE "install the rdkit in the source tree (former behavior)" ON )
option(RDK_INSTALL_DLLS_MSVC "install the rdkit DLLs when using MSVC" OFF)
option(RDK_INSTALL_STATIC_LIBS "install the rdkit static libraries" ON )
option(RDK_INSTALL_PYTHON_TESTS "install the rdkit Python tests with the wrappers" OFF )
option(RDK_BUILD_THREADSAFE_SSS "enable thread-safe substructure searching" ON )
option(RDK_BUILD_SLN_SUPPORT "include support for the SLN format" ON )
option(RDK_TEST_MULTITHREADED "run some tests of multithreading" ON )
option(RDK_BUILD_SWIG_JAVA_WRAPPER "build the SWIG JAVA wrappers (does nothing if RDK_BUILD_SWIG_WRAPPERS is not set)" ON )
option(RDK_BUILD_SWIG_CSHARP_WRAPPER "build the experimental SWIG C# wrappers (does nothing if RDK_BUILD_SWIG_WRAPPERS is not set)" OFF )
option(RDK_SWIG_STATIC "statically link rdkit libraries into the SWIG wrappers" ON )
option(RDK_TEST_MMFF_COMPLIANCE "run MMFF compliance tests (requires tar/gzip)" ON )
option(RDK_BUILD_CPP_TESTS "build the c++ tests (disabing can speed up builds" ON)
option(RDK_USE_FLEXBISON "use flex/bison, if available, to build the SMILES/SMARTS/SLN parsers" OFF)
option(RDK_TEST_COVERAGE "Use G(L)COV to compute test coverage" OFF)
option(RDK_USE_BOOST_SERIALIZATION "Use the boost serialization library if available" ON)
option(RDK_USE_BOOST_STACKTRACE "use boost::stacktrace to do more verbose invariant output (linux only)" ON)
option(RDK_BUILD_TEST_GZIP "Build the gzip'd stream test" OFF)
option(RDK_OPTIMIZE_POPCNT "Use SSE4.2 popcount instruction while compiling." ON)
option(RDK_USE_STRICT_ROTOR_DEFINITION "Use the most strict rotatable bond definition" ON)
option(RDK_BUILD_DESCRIPTORS3D "Build the 3D descriptors calculators, requires Eigen3 to be installed" ON)
option(RDK_BUILD_FREESASA_SUPPORT "build the rdkit freesasa wrapper" OFF )
option(RDK_BUILD_COORDGEN_SUPPORT "build the rdkit coordgen wrapper" ON )
option(RDK_BUILD_MAEPARSER_SUPPORT "build the rdkit MAE parser wrapper" ON )
option(RDK_BUILD_MOLINTERCHANGE_SUPPORT "build in support for CommonChem molecule interchange" ON )
option(RDK_BUILD_YAEHMOP_SUPPORT "build support for the YAeHMOP wrapper" OFF)
option(RDK_BUILD_STRUCTCHECKER_SUPPORT "build in support for the StructChecker alpha (not recommended, use the MolVS integration instead)" OFF )
option(RDK_USE_URF "Build support for Florian Flachsenberg's URF library" ON)
option(RDK_INSTALL_DEV_COMPONENT "install libraries and headers" ON)
option(RDK_USE_BOOST_REGEX "use boost::regex instead of std::regex (needed for systems with g++-4.8)" OFF)
option(RDK_USE_BOOST_IOSTREAMS "use boost::iostreams" ON)
option(RDK_BUILD_MINIMAL_LIB "build the minimal RDKit wrapper (for the JS bindings)" OFF)
option(RDK_BUILD_FUZZ_TARGETS "build the fuzz targets" OFF)

set(RDK_BOOST_VERSION "1.58.0")

if(NOT MSVC)
  if(RDK_OPTIMIZE_POPCNT)
    if(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64")
      message("RDK_OPTIMIZE_POPCNT is not available on aarch64")
    elseif(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "ppc64le")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcntb")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
    endif()
  endif()
endif()

if(NOT UNIX OR APPLE)
  if(RDK_USE_BOOST_STACKTRACE)
    MESSAGE("Disabling boost::stacktrace on non-linux platform")
    set(RDK_USE_BOOST_STACKTRACE OFF)
  endif(RDK_USE_BOOST_STACKTRACE)
endif()

if(RDK_BUILD_SWIG_WRAPPERS!=ON)
  set(RDK_BUILD_SWIG_JAVA_WRAPPER OFF)
  set(RDK_BUILD_SWIG_CSHARP_WRAPPER OFF)
endif()

if(RDK_SWIG_STATIC AND RDK_BUILD_SWIG_WRAPPERS)
  if(NOT RDK_INSTALL_STATIC_LIBS)
    message("Enabling RDK_INSTALL_STATIC_LIBS because RDK_SWIG_STATIC is set.")
    set(RDK_INSTALL_STATIC_LIBS ON CACHE BOOL "install the rdkit static libraries" FORCE)
  endif(NOT RDK_INSTALL_STATIC_LIBS)
endif()

if(NOT (MSVC AND (NOT RDK_INSTALL_DLLS_MSVC)))
  add_definitions(-DRDKIT_DYN_LINK)
endif()

if(RDK_PGSQL_STATIC AND RDK_BUILD_PGSQL)
  if(NOT RDK_INSTALL_STATIC_LIBS)
    message("Enabling RDK_INSTALL_STATIC_LIBS because RDK_PGSQL_STATIC is set.")
    set(RDK_INSTALL_STATIC_LIBS ON CACHE BOOL "install the rdkit static libraries" FORCE)
  endif(NOT RDK_INSTALL_STATIC_LIBS)
endif()


include(TestBigEndian)
TEST_BIG_ENDIAN(RDK_BIG_ENDIAN)

# At build time put runtime binaries in the bin subdirectory
set(RDK_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(RDK_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(RDK_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
# Also place the python extension where the source tree would expect them
set(RDK_PYTHON_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/rdkit")

#-------
# Config variables:
set(RDKit_Year "2020")
set(RDKit_Month "09")
set(RDKit_Revision "1")
set(RDKit_RevisionModifier "dev1")
set(RDKit_ABI "1")

# we need an integer version of the month later, so remove the zero padding
# if there is any
string(REGEX REPLACE "^0" "" RDKit_intMonth ${RDKit_Month} )

set(RDKit_CodeDir "${CMAKE_CURRENT_SOURCE_DIR}/Code")
set(RDKit_ExternalDir "${CMAKE_CURRENT_SOURCE_DIR}/External")

#include catch
add_subdirectory(${RDKit_ExternalDir}/catch)
include_directories(${CATCH_INCLUDE_DIR})
message("CATCH: ${CATCH_INCLUDE_DIR}")

if(RDK_INSTALL_INTREE)
  set(RDKit_BinDir "${CMAKE_SOURCE_DIR}/bin")
  set(RDKit_LibDir "${CMAKE_SOURCE_DIR}/lib")
  set(RDKit_HdrDir "Code") # used in rdkit-config.cmake, path prefix not needed.
  set(RDKit_ShareDir "${CMAKE_SOURCE_DIR}")
else(RDK_INSTALL_INTREE)
  set(RDKit_BinDir "bin")
  set(RDKit_LibDir "lib${LIB_SUFFIX}")
  set(RDKit_HdrDir "include/rdkit")
  set(RDKit_ShareDir "share/RDKit")
endif(RDK_INSTALL_INTREE)

if(RDK_BUILD_RPATH_SUPPORT)
  # use, i.e. don't skip the full RPATH for the build tree
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  message("CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system
  # directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
  "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
     SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  ENDIF("${isSystemDir}" STREQUAL "-1")
endif()



if(RDK_BUILD_SWIG_WRAPPERS)
  set(RDKit_JavaLibDir "${RDKit_ExternalDir}/java_lib")
endif(RDK_BUILD_SWIG_WRAPPERS)


if(RDK_BUILD_COORDGEN_SUPPORT)
    add_definitions(-DRDK_BUILD_COORDGEN_SUPPORT)
    include_directories(${RDKit_ExternalDir})
endif()

if(RDK_BUILD_MAEPARSER_SUPPORT)
    add_definitions(-DRDK_BUILD_MAEPARSER_SUPPORT)
    include_directories(${RDKit_ExternalDir})
endif()

if(RDK_USE_URF)
  ADD_DEFINITIONS("-DRDK_USE_URF")
  include_directories(${RDKit_ExternalDir}/RingFamilies/RingDecomposerLib/src/RingDecomposerLib)
endif()


if(CMAKE_SIZEOF_VOID_P MATCHES 4)
  target_compile_definitions(rdkit_base INTERFACE "-DRDK_32BIT_BUILD")
else()
  target_compile_definitions(rdkit_base INTERFACE "-DRDK_64BIT_BUILD")
endif()

if(MINGW)
  target_compile_definitions(rdkit_base INTERFACE "-DBOOST_SYSTEM_NO_DEPRECATED")
endif(MINGW)

# defines macros: rdkit_python_extension, rdkit_test
include(RDKitUtils)
install(EXPORT ${RDKit_EXPORTED_TARGETS} DESTINATION ${RDKit_LibDir}/cmake/rdkit)
install(TARGETS rdkit_base EXPORT ${RDKit_EXPORTED_TARGETS}
        COMPONENT dev )

# create and install package configuration and version files
configure_file (
    ${RDKit_SOURCE_DIR}/rdkit-config.cmake.in
    ${RDKit_BINARY_DIR}/rdkit-config.cmake @ONLY)

configure_file (
    ${RDKit_SOURCE_DIR}/rdkit-config-version.cmake.in
    ${RDKit_BINARY_DIR}/rdkit-config-version.cmake @ONLY)

install(FILES
    ${RDKit_BINARY_DIR}/rdkit-config.cmake
    ${RDKit_BINARY_DIR}/rdkit-config-version.cmake
    DESTINATION ${RDKit_LibDir}/cmake/rdkit)

# disable some warnings that we don't care about
if(MSVC)
  # These should occasionally be enabled to check what's going on
  # 4267 conversion with possible loss of data.
  # 4305 truncation from 'double' to 'const float'
  # 4244 conversion from 'uint64_t' to 'unsigned int', possible loss of data
  target_compile_options(rdkit_base INTERFACE "/wd4267" "/wd4305" "/wd4244")
  target_compile_definitions(rdkit_base INTERFACE "_CRT_SECURE_NO_WARNINGS")
endif(MSVC)


# extra boost versions
if(MSVC)
  # FIX: do we still need this?
  target_compile_definitions(rdkit_base INTERFACE "-DBOOST_ALL_NO_LIB")
endif(MSVC)

if(RDK_BUILD_INCHI_SUPPORT)
  find_package(Inchi)
endif(RDK_BUILD_INCHI_SUPPORT)

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_library(rdkit_py_base INTERFACE)

  #-------
  # pull in python:
  find_package(PythonInterp)
  if (PYTHONINTERP_FOUND AND NOT Python_ADDITIONAL_VERSIONS)
    set(Python_ADDITIONAL_VERSIONS "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
  endif (PYTHONINTERP_FOUND AND NOT Python_ADDITIONAL_VERSIONS)

  if(PYTHON_VERSION_MAJOR EQUAL 2)
    message(FATAL_ERROR "ERROR: The RDKit no longer supports Python 2.\nIf you have Python 3 installed, you can tell cmake where to find it using the -DPYTHON_EXECUTABLE argument.")
  endif()

  find_package(PythonLibs)
  target_include_directories(rdkit_base INTERFACE ${PYTHON_INCLUDE_DIR})

  # determine linkage of python
  execute_process(
      COMMAND
      ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_config_var('Py_ENABLE_SHARED'))"
      OUTPUT_VARIABLE Py_ENABLE_SHARED
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(WIN32 OR "${Py_ENABLE_SHARED}" STREQUAL "1")
    target_link_libraries(rdkit_py_base INTERFACE ${PYTHON_LIBRARIES} )
  endif()
  
  find_package(NumPy REQUIRED)
  target_include_directories(rdkit_base INTERFACE ${PYTHON_NUMPY_INCLUDE_PATH})

  if(PYTHON_VERSION_MAJOR EQUAL 3)
    # Find boost-python3 using name specified as command line option then fall back to commonly used names
    set(RDK_BOOST_PYTHON3_NAME "python3" CACHE STRING "Name of the boost python3 library. If installed as libboost_python-xxx.so, use python-xxx.")
    list(APPEND Boost_Python_Names "${RDK_BOOST_PYTHON3_NAME}" "python-py3${PYTHON_VERSION_MINOR}" "python3")
  endif(PYTHON_VERSION_MAJOR EQUAL 3)

  # Boost 1.67+ uses a version suffix like "python36" or "python27"
  list(APPEND Boost_Python_Names "python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")

  # Try each potential boost-python name until one works
  foreach(Boost_Python_LibN ${Boost_Python_Names})
    find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS "${Boost_Python_LibN}" QUIET)
    if(Boost_FOUND)
      set(Boost_Python_Lib ${Boost_Python_LibN})
      break()
    endif()
  endforeach()
  # Finally just try "python" and hope it is a compatible version
  if(NOT Boost_FOUND)
    find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS python REQUIRED)
    set(Boost_Python_Lib "python")
  endif()

  target_link_libraries(rdkit_py_base INTERFACE Boost::${Boost_Python_Lib})

  if(RDK_INSTALL_INTREE)
    set(RDKit_PythonDir "${CMAKE_SOURCE_DIR}/rdkit")
  else(RDK_INSTALL_INTREE)
    if (NOT PYTHON_INSTDIR)
      # Determine correct installation directory for Python bindings
      # From http://plplot.svn.sourceforge.net/viewvc/plplot/trunk/cmake/modules/python.cmake?revision=11014
      execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; from distutils import sysconfig; print(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
        OUTPUT_VARIABLE PYTHON_INSTDIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    endif (NOT PYTHON_INSTDIR)
    message("Python Install directory ${PYTHON_INSTDIR}")
    set(RDKit_PythonDir "${PYTHON_INSTDIR}/rdkit")
    if(RDK_INSTALL_PYTHON_TESTS)
      install(DIRECTORY rdkit DESTINATION ${PYTHON_INSTDIR}
        COMPONENT python
        PATTERN ".svn" EXCLUDE
        PATTERN "CMake*" EXCLUDE
        PATTERN "Basement" EXCLUDE
        )
    else(RDK_INSTALL_PYTHON_TESTS)
      install(DIRECTORY rdkit DESTINATION ${PYTHON_INSTDIR}
        COMPONENT python
        PATTERN ".svn" EXCLUDE
        PATTERN "test_data" EXCLUDE
        PATTERN "testData" EXCLUDE
        PATTERN "test_list*" EXCLUDE
        PATTERN "CMake*" EXCLUDE
        PATTERN "Basement" EXCLUDE
        PATTERN "UnitTest*" EXCLUDE
        )
    endif(RDK_INSTALL_PYTHON_TESTS)

  endif(RDK_INSTALL_INTREE)

  if(NOT WIN32)
    # See https://bugs.python.org/msg277944
    # The "command to create shared modules". Used as variable in the "Makefile (and similar) templates to build python modules"
    # for both in-python and third party modules. Initialized to hold the value which works for third party modules to link
    # against the _installed_ python.
    # We strip off the first word though (which will be the compiler name).
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_config_var('LDSHARED').split(' ', 1)[1])"
        OUTPUT_VARIABLE PYTHON_LDSHARED
        OUTPUT_STRIP_TRAILING_WHITESPACE
     )
     message("PYTHON Py_ENABLE_SHARED: ${Py_ENABLE_SHARED}")
     message("PYTHON USING LINK LINE: ${PYTHON_LDSHARED}")
  endif()

  install(TARGETS rdkit_py_base EXPORT ${RDKit_EXPORTED_TARGETS}
          COMPONENT dev )

else(RDK_BUILD_PYTHON_WRAPPERS)
  find_package(Boost ${RDK_BOOST_VERSION} REQUIRED)
endif(RDK_BUILD_PYTHON_WRAPPERS)

# Detect clang, which masquerades as gcc.  CMake 2.6 doesn't know how to
# detect it.
string(REGEX MATCH "clang" CMAKE_COMPILER_IS_CLANG "${CMAKE_C_COMPILER}")

if(RDK_BUILD_DESCRIPTORS3D)
  find_package(Eigen3)
  if(EIGEN3_FOUND)
    target_compile_definitions(rdkit_base INTERFACE "-DRDK_HAS_EIGEN3" "-DRDK_BUILD_DESCRIPTORS3D")
    target_link_libraries(rdkit_base INTERFACE Eigen3::Eigen)
    set(RDK_HAS_EIGEN ON)
  else(EIGEN3_FOUND)
    message("Eigen3 not found, disabling the Descriptors3D build.")
    set(RDK_BUILD_DESCRIPTORS3D OFF)
  endif(EIGEN3_FOUND)
endif(RDK_BUILD_DESCRIPTORS3D)

if(RDK_BUILD_ANI)
  find_package(Eigen3)
  if(EIGEN3_FOUND)
    set(RDK_HAS_EIGEN ON)
    if(NOT WIN32)
      target_compile_definitions(rdkit_base INTERFACE "-DRDK_HAS_EIGEN3" "-DRDK_BUILD_ANI")
      target_link_libraries(rdkit_base INTERFACE Eigen3::Eigen)
      set(ANIParamsDir "${CMAKE_CURRENT_SOURCE_DIR}/Data/ANIParams")
      if(NOT EXISTS "${ANIParamsDir}/ANI-1x/selfEnergies")
        set(RELEASE_NO "1.0")
        downloadAndCheckMD5("https://github.com/rdkit/rdkit_ANI_weights/archive/v${RELEASE_NO}.tar.gz"
                "${CMAKE_CURRENT_SOURCE_DIR}/ANIWeights-v${RELEASE_NO}.tar.gz" "")
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar zxf
          ${CMAKE_CURRENT_SOURCE_DIR}/ANIWeights-v${RELEASE_NO}.tar.gz
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        file(RENAME "${CMAKE_CURRENT_SOURCE_DIR}/rdkit_ANI_weights-${RELEASE_NO}/Params" "${CMAKE_CURRENT_SOURCE_DIR}/Data/ANIParams")
        file(REMOVE_RECURSE "${CMAKE_CURRENT_SOURCE_DIR}/rdkit_ANI_weights-${RELEASE_NO}")
        file(REMOVE "${CMAKE_CURRENT_SOURCE_DIR}/ANIWeights-v${RELEASE_NO}.tar.gz")
      endif()
    endif()
  else(EIGEN3_FOUND)
    message("Eigen3 not found, disabling the ANI build.")
    set(RDK_BUILD_ANI OFF)
  endif(EIGEN3_FOUND)
  if(WIN32)
    message("ANI not supported on Windows")
    set(RDK_BUILD_ANI OFF)
  endif(WIN32)
endif(RDK_BUILD_ANI)

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)
if(RDK_BUILD_THREADSAFE_SSS)
  target_compile_definitions(rdkit_base INTERFACE -DRDK_THREADSAFE_SSS)
  if(RDK_TEST_MULTITHREADED)
    target_compile_definitions(rdkit_base INTERFACE "-DRDK_TEST_MULTITHREADED")
  endif()
  target_link_libraries(rdkit_base INTERFACE Threads::Threads)
else()
  if(RDK_TEST_MULTITHREADED)
    message("RDK_TEST_MULTITHREADED does not make sense without RD_BUILD_THREADSAFE_SSS, disabling it.")
    set(RDK_TEST_MULTITHREADED OFF)
  endif()
endif()

if(RDK_USE_BOOST_SERIALIZATION)
    set(T_LIBS ${Boost_LIBRARIES})
    find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS serialization)
    if (Boost_SERIALIZATION_LIBRARY)
      set(Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})
    else()
      set(Boost_LIBRARIES ${T_LIBS})
    endif()
  target_compile_definitions(rdkit_base INTERFACE -DRDK_USE_BOOST_SERIALIZATION -DBOOST_SERIALIZATION_DYN_LINK)
endif()


# set the boost include directories and linkage:
target_include_directories(rdkit_base INTERFACE
                           $<BUILD_INTERFACE:${RDKit_CodeDir}>
                           $<INSTALL_INTERFACE:${RDKit_HdrDir}>
                           ${Boost_INCLUDE_DIR}
                           )
target_link_libraries(rdkit_base INTERFACE ${BOOST_LIBRARIES})


if(RDK_BUILD_PGSQL)
  find_package(PostgreSQL REQUIRED)
endif(RDK_BUILD_PGSQL)

# setup our compiler flags:
if (RDK_TEST_COVERAGE)
 if (NOT RDK_USE_FLEXBISON)
   message(FATAL_ERROR, "Test coverage doesn't current work unless FLEX and BISON are run.")
 endif(NOT RDK_USE_FLEXBISON)

 message("====== Installing test coverage support ======")
 message("  To run:")
 message("   make install")
 message("   make RDKit_coverage")
 message("")
 message(" open <build_dir>/coverage/index.html")
 message("")
 message(" If any of the RDKit tests fail, coverage will probably not be generated.")
 message("====== Installing test coverage support ======")
 message("")

 INCLUDE(CodeCoverage)
 SET(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
 SET(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
 message("== setup_target_for_coverage(${PROJECT_NAME}_coverage test coverage)")
 setup_target_for_coverage(${PROJECT_NAME}_coverage ctest coverage)

else(RDK_TEST_COVERAGE)
 if(CMAKE_COMPILER_IS_GNUCXX)
   # We don't on C_FLAGS warnings to keep Avalon compiling quiet
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated -Wno-unused-function -fno-strict-aliasing")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-unused-function -fno-strict-aliasing -Wall -Wextra")
  if (NOT CYGWIN)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
 endif()
 if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   # We don't on C_FLAGS warnings to keep Avalon compiling quiet
   SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-array-bounds -fPIC -Wno-parentheses -Wno-logical-op-parentheses -Wno-format")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-deprecated -Wno-unused-function -fno-strict-aliasing -Wno-format -Wno-logical-op-parentheses -fPIC")
   if(APPLE)
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
   endif()
 endif()
endif(RDK_TEST_COVERAGE)

if(NOT RDK_INSTALL_INTREE)
  install(DIRECTORY Data DESTINATION
      ${RDKit_ShareDir}
      COMPONENT data
      PATTERN ".svn" EXCLUDE
      )
  install(DIRECTORY Docs DESTINATION
      ${RDKit_ShareDir}
      COMPONENT docs
      PATTERN ".svn" EXCLUDE
      PATTERN "_build" EXCLUDE
      )
  install(DIRECTORY Contrib DESTINATION
      ${RDKit_ShareDir}
      COMPONENT extras
      PATTERN ".svn" EXCLUDE
      )
  if(RDK_BUILD_PYTHON_WRAPPERS)
    install(DIRECTORY Projects DESTINATION
        ${RDKit_ShareDir}
        COMPONENT extras
        PATTERN ".svn" EXCLUDE
        PATTERN "test_list*" EXCLUDE
        PATTERN "CMake*" EXCLUDE
        )
    if(RDK_INSTALL_PYTHON_TESTS)
      install(DIRECTORY Scripts DESTINATION
        ${RDKit_ShareDir}
        COMPONENT data
        PATTERN ".svn" EXCLUDE
        )
    else(RDK_INSTALL_PYTHON_TESTS)
      install(DIRECTORY Scripts DESTINATION
        ${RDKit_ShareDir}
        COMPONENT data
        PATTERN ".svn" EXCLUDE
        PATTERN "run_python_tests.py" EXCLUDE
        )
    endif(RDK_INSTALL_PYTHON_TESTS)


  endif(RDK_BUILD_PYTHON_WRAPPERS)
  install(FILES README license.txt
          DESTINATION ${RDKit_ShareDir}
          COMPONENT base
          )
endif(NOT RDK_INSTALL_INTREE)

if (RDK_SQUASH_MVC_SECURE_WARNINGS)
  MESSAGE("== Squashing MSVC Security warnings (do not use during development)")
  target_compile_definitions(rdkit_base INTERFACE "-D_SCL_SECURE_NO_WARNINGS")
  target_compile_definitions(rdkit_base INTERFACE "/wd4996")
endif(RDK_SQUASH_MVC_SECURE_WARNINGS)

if(RDK_USE_STRICT_ROTOR_DEFINITION)
    MESSAGE("== Using strict rotor definition")
    target_compile_definitions(rdkit_base INTERFACE "-DRDK_USE_STRICT_ROTOR_DEFINITION")
endif()

if (MSVC)
  # disable warnings:
  #  - 4267: conversion from 'size_t' to 'unsigned int', possible loss of data
  add_definitions( "/wd4267" )
endif(MSVC)



add_subdirectory(External)
add_subdirectory(Code)

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(Projects)
  add_subdirectory(rdkit)
endif(RDK_BUILD_PYTHON_WRAPPERS)

if(RDK_BUILD_CONTRIB)
  add_subdirectory(Contrib)
endif(RDK_BUILD_CONTRIB)

# Memory testing setup
FIND_PROGRAM(MEMORYCHECK_COMMAND valgrind)
CONFIGURE_FILE(CTestCustom.ctest.in ${RDKit_BINARY_DIR}/CTestCustom.ctest)

# Packaging
SET(CPACK_GENERATOR "TGZ;DEB;RPM")
SET(CPACK_RPM_COMPONENT_INSTALL ON)
SET(CPACK_DEB_COMPONENT_INSTALL ON)
SET(CPACK_MONOLITHIC_INSTALL OFF)

SET(CPACK_COMPONENTS_ALL runtime base data docs dev python extras)
set(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
set(CPACK_COMPONENT_BASE_GROUP "Runtime")
set(CPACK_COMPONENT_DATA_GROUP "Runtime")
set(CPACK_COMPONENT_DOCS_GROUP "Runtime")
set(CPACK_COMPONENT_DEV_GROUP "Development")
set(CPACK_COMPONENT_PYTHON_GROUP "Python")
set(CPACK_COMPONENT_EXTRAS_GROUP "Extras")

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "RDKit")
SET(CPACK_PACKAGE_VENDOR "rdkit.org")
SET(CPACK_PACKAGE_CONTACT "greg.landrum@gmail.com")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR ${RDKit_Year})
SET(CPACK_PACKAGE_VERSION_MINOR ${RDKit_Month})
SET(CPACK_PACKAGE_VERSION_PATCH ${RDKit_Revision}${RDKit_RevisionModifier})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
#IF(WIN32 AND NOT UNIX)
#  # There is a bug in NSI that does not handle full unix paths properly. Make
#  # sure there is at least one set of four (4) backlasshes.
#  SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
#  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} My Famous Project")
#  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.my-project-home-page.org")
#  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.my-personal-home-page.com")
#  SET(CPACK_NSIS_CONTACT "me@my-personal-home-page.com")
#  SET(CPACK_NSIS_MODIFY_PATH ON)
#ELSE(WIN32 AND NOT UNIX)
#  SET(CPACK_STRIP_FILES "bin/MyExecutable")
#  SET(CPACK_SOURCE_STRIP_FILES "")
#ENDIF(WIN32 AND NOT UNIX)
#SET(CPACK_PACKAGE_EXECUTABLES "MyExecutable" "My Executable")


SET(CPACK_SET_DESTDIR ON)

INCLUDE(CPack)
