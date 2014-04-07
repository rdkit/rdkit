# FindRDKit.cmake
# Placed in the public domain by NextMove Software in 2013
# Try to find RDKit headers and libraries
# Defines:
#
#  RDKIT_FOUND - system has RDKit
#  RDKIT_INCLUDE_DIR - the RDKit include directory
#  RDKIT_LIBRARIES - Link these to use RDKit

if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
  # in cache already or user-specified
  set(RDKIT_FOUND TRUE)

else()

  if(NOT RDKIT_INCLUDE_DIR)
    if(WIN32)
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
        ${RDKIT_DIR}\\Code
        $ENV{RDKIT_INCLUDE_DIR}
        $ENV{RDKIT_INCLUDE_PATH}
        $ENV{RDKIT_BASE}\\Code
        $ENV{RDBASE}\\Code
        C:\\RDKit\\include
        C:\\RDKit\\Code
      )
    else()
      find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
        PATHS
          ${RDKIT_DIR}/Code
          $ENV{RDKIT_INCLUDE_DIR}
          $ENV{RDKIT_INCLUDE_PATH}
          $ENV{RDKIT_BASE}/Code
          $ENV{RDBASE}/Code
          /usr/local/rdkit/include/Code
          /usr/local/rdkit/include
          /usr/local/rdkit/Code
          ~/rdkit/Code
      )
    endif()
    if(RDKIT_INCLUDE_DIR)
       message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_DIR}")
    endif()
  endif()

  if(NOT RDKIT_LIBRARIES)
    find_library(FILEPARSERS_LIB NAMES FileParsers
      PATHS
        ${RDKIT_DIR}/lib
        $ENV{RDKIT_LIB_DIR}
        $ENV{RDKIT_LIB_PATH}
        $ENV{RDKIT_LIBRARIES}
        $ENV{RDKIT_BASE}/lib
        $ENV{RDBASE}/lib
        /usr/local/rdkit/lib
        ~/rdkit/lib
        $ENV{LD_LIBRARY_PATH}
    )
    if(FILEPARSERS_LIB)
       GET_FILENAME_COMPONENT(RDKIT_LIBRARY_DIR ${FILEPARSERS_LIB} PATH)
       message(STATUS "Found RDKit libraries at ${RDKIT_LIBRARY_DIR}")

      # Note that the order of the following libraries is significant!!
      find_library(SMILESPARSE_LIB NAMES SmilesParse
                                   HINTS ${RDKIT_LIBRARY_DIR})
      find_library(DEPICTOR_LIB NAMES Depictor
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(GRAPHMOL_LIB NAMES GraphMol
                                HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGEOMETRYLIB_LIB NAMES RDGeometryLib
                                     HINTS ${RDKIT_LIBRARY_DIR})
      find_library(RDGENERAL_LIB NAMES RDGeneral
                                 HINTS ${RDKIT_LIBRARY_DIR})
      set (RDKIT_LIBRARIES ${FILEPARSERS_LIB} ${SMILESPARSE_LIB}
              ${DEPICTOR_LIB} ${GRAPHMOL_LIB} ${RDGEOMETRYLIB_LIB}
              ${RDGENERAL_LIB})
    endif()
    if(RDKIT_LIBRARIES)
            message(STATUS "Found RDKit library files at ${RDKIT_LIBRARIES}")
    endif()
  endif()

  if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
    set(RDKIT_FOUND TRUE)
  endif()

  mark_as_advanced(RDKIT_INCLUDE_DIR RDKIT_LIBRARIES)
endif()