if(RDK_BUILD_YAEHMOP_SUPPORT)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-DRDK_BUILD_YAEHMOP_SUPPORT)

include(ExternalProject)

if(CMAKE_COMPILER_IS_GNUCXX AND NOT CYGWIN)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
endif()

include_directories( ${RDKit_ExternalDir}/YAeHMOP )

ExternalProject_Add(yaehmop_project
  GIT_REPOSITORY https://github.com/greglandrum/yaehmop.git
  GIT_TAG master
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/yaehmop"
  SOURCE_SUBDIR "tightbind"
  CMAKE_ARGS -DUSE_BLAS_LAPACK=OFF -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
  TEST_COMMAND "")

include_directories(${PROJECT_BINARY_DIR}/include)
link_directories(${PROJECT_BINARY_DIR}/lib)

set(EHT_PARAM_FILE ${CMAKE_CURRENT_SOURCE_DIR}/yaehmop/tightbind/eht_parms.dat )
install(FILES ${EHT_PARAM_FILE}
        DESTINATION ${RDKit_ShareDir}/Data
        COMPONENT data)

message("include_dirs: ${PROJECT_BINARY_DIR}/include")
message("link_dirs: ${PROJECT_BINARY_DIR}/lib")

add_definitions(-DRDKIT_AVALONLIB_BUILD)
rdkit_library(EHTLib EHTTools.cpp SHARED LINK_LIBRARIES yaehmop_eht GraphMol DataStructs RDGeometryLib RDGeneral )
add_dependencies(EHTLib yaehmop_project)
rdkit_headers(EHTTools.h DEST GraphMol)
rdkit_catch_test(testEHTLib1 test1.cpp
           LINK_LIBRARIES EHTLib yaehmop_eht 
	   FileParsers SmilesParse GraphMol DataStructs RDGeometryLib RDGeneral)

add_subdirectory(Wrap)

endif(RDK_BUILD_YAEHMOP_SUPPORT)
