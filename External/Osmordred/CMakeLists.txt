if(NOT RDK_BUILD_OSMORDRED_SUPPORT)
  return()
endif()

find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
    PATH_SUFFIXES /eigen
    REQUIRED
)
include_directories(Osmordred ${EIGEN3_INCLUDE_DIR})

# The FindLAPACK distributed with cmake does not find the include directories for some
#  bizarre reason.  So use the wozune distributed with eigen
LIST(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)

if(NOT LAPACKE_FOUND)
  message(fatal_error "-- LAPACKE not found")
endif()

message(status "-- LAPACKE INCLUDE ${LAPACKE_INCLUDE_DIRS}")
message(status "-- LAPACKE LIB ${LAPACKE_LIBRARIES}")

include_directories(Osmordred LAPACK::LAPACK ${LAPACK_INCLUDE_DIRS})
rdkit_library(Osmordred
    Osmordred.cpp
    LINK_LIBRARIES DataStructs Descriptors Fingerprints GraphMol
                   PartialCharges RDGeneral SmilesParse Subgraphs SubstructMatch
		   ${LAPACKE_LIBRARIES})

target_compile_definitions(Osmordred PRIVATE RDKIT_OSMORDRED_BUILD)

if(MSVC)
else()
 target_compile_options(Osmordred PRIVATE -w -Wno-error)
      if(RDK_INSTALL_STATIC_LIBS)
            target_compile_options(Osmordred_static PRIVATE -w -Wno-error)
      endif(RDK_INSTALL_STATIC_LIBS)
endif()

rdkit_headers(Osmordred.h DEST Osmordred)

rdkit_catch_test(testOsmordred test_osmordred.cpp LINK_LIBRARIES Osmordred )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
