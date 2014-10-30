cmake_minimum_required(VERSION 2.6)
if(CMAKE_MAJOR_VERSION EQUAL 2 AND ((CMAKE_MINOR_VERSION EQUAL 8 AND CMAKE_PATCH_VERSION GREATER 1) OR CMAKE_MINOR_VERSION GREATER 8))
  cmake_policy(SET CMP0015 OLD)
endif()

project (RDKit)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

enable_testing()

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/Code/cmake/Modules/")

option(RDK_BUILD_SWIG_WRAPPERS "build the SWIG wrappers" OFF )
option(RDK_BUILD_PYTHON_WRAPPERS "build the standard python wrappers" ON )
option(RDK_BUILD_COMPRESSED_SUPPLIERS "build in support for compressed MolSuppliers" OFF )
option(RDK_BUILD_INCHI_SUPPORT "build the rdkit inchi wrapper" OFF )
option(RDK_BUILD_AVALON_SUPPORT "install support for the avalon toolkit. Use the variable AVALONTOOLS_DIR to set the location of the source." OFF )
option(RDK_BUILD_CONTRIB "build the Contrib directory" OFF )
option(RDK_INSTALL_INTREE "install the rdkit in the source tree (former behavior)" ON )
option(RDK_INSTALL_STATIC_LIBS "install the rdkit static libraries" ON )
option(RDK_BUILD_THREADSAFE_SSS "enable thread-safe substructure searching (requires boost.thread)" OFF )
option(RDK_BUILD_SLN_SUPPORT "include support for the SLN format" ON )
option(RDK_TEST_MULTITHREADED "run some tests of multithreading" OFF )
option(RDK_BUILD_SWIG_JAVA_WRAPPER "build the SWIG JAVA wrappers (does nothing if RDK_BUILD_SWIG_WRAPPERS is not set)" ON )
option(RDK_BUILD_SWIG_CSHARP_WRAPPER "build the experimental SWIG C# wrappers (does nothing if RDK_BUILD_SWIG_WRAPPERS is not set)" OFF )
option(RDK_TEST_MMFF_COMPLIANCE "run MMFF compliance tests (requires tar/gzip)" ON )
option(RDK_BUILD_CPP_TESTS "build the c++ tests (disabing can speed up builds" ON)
option(RDK_USE_FLEXBISON "use flex/bison, if available, to build the SMILES/SMARTS/SLN parsers" OFF)

if(RDK_BUILD_SWIG_WRAPPERS!=ON)
  set(RDK_BUILD_SWIG_JAVA_WRAPPER OFF)
  set(RDK_BUILD_SWIG_CSHARP_WRAPPER OFF)
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
set(RDKit_Year "2015")
set(RDKit_Month "03")
set(RDKit_Revision "1pre")
set(RDKit_ABI "1")

set(RDKit_CodeDir "${CMAKE_SOURCE_DIR}/Code")
set(RDKit_ExternalDir "${CMAKE_SOURCE_DIR}/External")

if(RDK_INSTALL_INTREE)
  set(RDKit_BinDir "${CMAKE_SOURCE_DIR}/bin")
  set(RDKit_LibDir "${CMAKE_SOURCE_DIR}/lib")
  set(RDKit_HdrDir "Code") # used in rdkit-config.cmake, path prefix not needed.
else(RDK_INSTALL_INTREE)
  set(RDKit_BinDir "bin")
  set(RDKit_LibDir "lib${LIB_SUFFIX}")
  set(RDKit_HdrDir "include/rdkit")
  set(RDKit_ShareDir "share/RDKit")
endif(RDK_INSTALL_INTREE)

if(RDK_BUILD_SWIG_WRAPPERS)
  set(RDKit_JavaLibDir "${RDKit_ExternalDir}/java_lib")
endif(RDK_BUILD_SWIG_WRAPPERS)

if(CMAKE_SIZEOF_VOID_P MATCHES 4)
  ADD_DEFINITIONS("-DRDK_32BIT_BUILD")
else()
  ADD_DEFINITIONS("-DRDK_64BIT_BUILD")
endif()

# defines macros: rdkit_python_extension, rdkit_test
include(RDKitUtils)
install(EXPORT ${RDKit_EXPORTED_TARGETS} DESTINATION ${RDKit_LibDir})

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
    DESTINATION ${RDKit_LibDir})

# extra boost versions
if(MSVC)
  SET(Boost_ADDITIONAL_VERSIONS "1.48" "1.48.0" "1.45" "1.45.0" "1.44" "1.44.0" "1.43" "1.43.0" "1.42" "1.42.0" "1.41" "1.41.0" "1.40" "1.40.0")
endif(MSVC)

if(RDK_BUILD_PYTHON_WRAPPERS)
  #-------
  # pull in python:
  find_package(PythonInterp)
  find_package(PythonLibs)
  if(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
    include_directories(${PYTHON_INCLUDE_PATH})
  else(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
    include_directories(${PYTHON_INCLUDE_DIR})
  endif(CMAKE_MAJOR_VERSION EQUAL 2 AND CMAKE_MINOR_VERSION EQUAL 6)
  find_package(NumPy REQUIRED)
  include_directories(${PYTHON_NUMPY_INCLUDE_PATH})

  if(PYTHON_VERSION_MAJOR EQUAL 3)
    # Find boost-python3 using name specified as command line option then fall back to commonly used names
    set(RDK_BOOST_PYTHON3_NAME "python3" CACHE STRING "Name of the boost python3 library. If installed as libboost_python-xxx.so, use python-xxx.")
    foreach(Boost_Python_Lib "${RDK_BOOST_PYTHON3_NAME}" "python-py3${PYTHON_VERSION_MINOR}" "python3")
      find_package(Boost 1.45.0 COMPONENTS "${Boost_Python_Lib}" QUIET)
      if(Boost_FOUND)
        break()
      endif()
    endforeach()
    # If boost-python3 not found, just find boost-python and hope that it is compatible with python3
    if(NOT Boost_FOUND)
      find_package(Boost 1.45.0 COMPONENTS python REQUIRED)
    endif()
  else(PYTHON_VERSION_MAJOR EQUAL 3)
     find_package(Boost 1.39.0 COMPONENTS python REQUIRED)
  endif(PYTHON_VERSION_MAJOR EQUAL 3)


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
    set(RDKit_PythonDir "${PYTHON_INSTDIR}/rdkit")
  endif(RDK_INSTALL_INTREE)

else(RDK_BUILD_PYTHON_WRAPPERS)
  find_package(Boost 1.39.0 REQUIRED)
endif(RDK_BUILD_PYTHON_WRAPPERS)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${RDKit_CodeDir})


# Detect clang, which masquerades as gcc.  CMake 2.6 doesn't know how to
# detect it.
string(REGEX MATCH "clang" CMAKE_COMPILER_IS_CLANG "${CMAKE_C_COMPILER}")

if(RDK_TEST_MULTITHREADED)
  ADD_DEFINITIONS("-DRDK_TEST_MULTITHREADED")
endif()

find_package (Threads)
if(RDK_BUILD_THREADSAFE_SSS)
  set(T_LIBS ${Boost_LIBRARIES})
  find_package(Boost 1.39.0 COMPONENTS thread system REQUIRED)
  set(Boost_LIBRARIES ${T_LIBS} ${Boost_LIBRARIES})
  ADD_DEFINITIONS("-DRDK_THREADSAFE_SSS")
  set(RDKit_THREAD_LIBS ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
else()
  set(RDKit_THREAD_LIBS ${CMAKE_THREAD_LIBS_INIT})
endif()


# setup our compiler flags:
if(CMAKE_COMPILER_IS_GNUCXX)
  ADD_DEFINITIONS("-Wno-deprecated -Wno-unused-function -fno-strict-aliasing -fPIC")
endif()
if(CMAKE_COMPILER_IS_CLANG)
  ADD_DEFINITIONS("-Wno-array-bounds -fPIC")
endif()

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
  endif(RDK_BUILD_PYTHON_WRAPPERS)
  install(FILES README license.txt
          DESTINATION ${RDKit_ShareDir}
          COMPONENT base
          )
endif(NOT RDK_INSTALL_INTREE)

add_subdirectory(External)
add_subdirectory(Code)

if(RDK_BUILD_PYTHON_WRAPPERS)
  add_subdirectory(Projects)
  add_subdirectory(rdkit)
endif(RDK_BUILD_PYTHON_WRAPPERS)

if(RDK_BUILD_CONTRIB)
  add_subdirectory(Contrib)
endif(RDK_BUILD_CONTRIB)

# Packaging
SET(CPACK_GENERATOR "TGZ;DEB;RPM")
SET(CPACK_RPM_COMPONENT_INSTALL ON)
SET(CPACK_DEB_COMPONENT_INSTALL ON)
SET(CPACK_MONOLITHIC_INSTALL OFF)

SET(CPACK_COMPONENTS_ALL headers static_libs dynamic_libs python data docs extras base)
set(CPACK_COMPONENT_DYNAMIC_LIBS_GROUP "Runtime")
set(CPACK_COMPONENT_BASE_GROUP "Runtime")
set(CPACK_COMPONENT_DATA_GROUP "Runtime")
set(CPACK_COMPONENT_DOCS_GROUP "Runtime")
set(CPACK_COMPONENT_HEADERS_GROUP "Development")
set(CPACK_COMPONENT_STATIC_LIBS_GROUP "Development")
set(CPACK_COMPONENT_PYTHON_GROUP "Python")
set(CPACK_COMPONENT_EXTRAS_GROUP "Extras")


SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "RDKit")
SET(CPACK_PACKAGE_VENDOR "rdkit.org")
SET(CPACK_PACKAGE_CONTACT "greg.landrum@gmail.com")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "2011")
SET(CPACK_PACKAGE_VERSION_MINOR "6")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
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
