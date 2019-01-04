if(RDK_BUILD_YAEHMOP_SUPPORT)

if(MSVC)
set(YAEHMOP_DIR "c:/Users/glandrum/Code/yaehmop/tightbind/")
set(YAEHMOP_LIB "${YAEHMOP_DIR}/build/Debug/yaehmop_eht.lib")
else(MSVC)
set(YAEHMOP_DIR "/mnt/c/Users/glandrum/Code/yaehmop/tightbind/")
set(YAEHMOP_LIB "${YAEHMOP_DIR}/build_lin/libyaehmop_eht.a")
endif(MSVC)
#find_library(YAEHMOP_LIB yaehmop_eht PATHS ${YAEHMOP_DIR} ${YAEHMOP_DIR}/build/Release)
#if(NOT YAEHMOP_LIB)
#message(FATAL_ERROR "YAeHMOP library not found")
#endif(NOT YAEHMOP_LIB)
set(YAEHMOP_INCLUDE_PATH ${YAEHMOP_DIR})

if(NOT DEFINED YAEHMOP_DIR)
  set(YAEHMOP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/yaehmop")
else()
  string(REGEX REPLACE "\\\\" "/" YAEHMOP_DIR ${YAEHMOP_DIR})
endif()

add_definitions(-DRDK_BUILD_YAEHMOP_SUPPORT)

include(ExternalProject)

ExternalProject_Add(yaehmop_project
  GIT_REPOSITORY https://github.com/greglandrum/yaehmop.git
  GIT_TAG master
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/yaehmop"
  CMAKE_ARGS -DUSE_BLAS_LAPACK=OFF  
  TEST_COMMAND "")
ExternalProject_Get_Property(yaehmop_project install_dir)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${YAEHMOP_INCLUDE_PATH})

add_definitions(-DRDKIT_AVALONLIB_BUILD)
rdkit_library(EHTLib EHTTools.cpp SHARED LINK_LIBRARIES ${YAEHMOP_LIB} GraphMol DataStructs RDGeometryLib RDGeneral )
rdkit_headers(EHTTools.h DEST GraphMol)
rdkit_catch_test(testEHTLib1 test1.cpp
           LINK_LIBRARIES EHTLib ${YAEHMOP_LIB}
	   FileParsers SmilesParse GraphMol DataStructs RDGeometryLib RDGeneral)


endif(RDK_BUILD_YAEHMOP_SUPPORT)
