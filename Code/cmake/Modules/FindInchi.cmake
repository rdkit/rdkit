# - Try to find Inchi lib
# Once done this will define
#
#  INCHI_FOUND - system has inchi lib
#  INCHI_INCLUDE_DIR - the inchi include directory
#  INCHI_LIBRARIES - the inchi library

# Taken from the open babel distribution
#  https://openbabel.svn.sourceforge.net/svnroot/openbabel/openbabel/trunk/cmake/modules/FindInchi.cmake
# Copyright (c) 2010 Marcus D. Hanwell, <marcus@cryos.org>
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if(INCHI_INCLUDE_DIR AND INCHI_LIBRARY)
  # in cache already
  set(INCHI_FOUND TRUE)
else()
  find_path(INCHI_INCLUDE_DIR NAMES inchi_api.h PATHS /usr/include/inchi )
  find_library(INCHI_LIBRARY NAMES inchi Inchi)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Inchi
	  "Could NOT find InChI in system locations" INCHI_LIBRARY
	  INCHI_INCLUDE_DIR)
  set(INCHI_INCLUDE_DIRS ${INCHI_INCLUDE_DIR} )
  set(INCHI_LIBRARIES ${INCHI_LIBRARY} )
  mark_as_advanced(INCHI_INCLUDE_DIR INCHI_LIBRARY)
endif()
set(CUSTOM_INCHI_PATH "${CMAKE_CURRENT_SOURCE_DIR}/External/INCHI-API")
# check whether we have custom InChI source
if(EXISTS ${CUSTOM_INCHI_PATH}/src/INCHI_BASE/src/ichican2.c)
  message("CUSTOM_INCHI_PATH = ${CUSTOM_INCHI_PATH}")
  message(STATUS  "Found InChI software locally")
  if(INCHI_FOUND)
    message(WARNING "** Local InChI software takes precedence when both system InChI and local InChI are found")
  endif(INCHI_FOUND)
else(EXISTS ${CUSTOM_INCHI_PATH}/src/INCHI_BASE/src/ichican2.c)
  if (INCHI_FOUND)
    include_directories(${INCHI_INCLUDE_DIR})
  else (INCHI_FOUND)
    # system InChI is missing, download it
    if(NOT DEFINED INCHI_URL)
      #set(INCHI_URL "https://www.inchi-trust.org/download/106/INCHI-1-SRC.zip")
      set(INCHI_URL "https://rdkit.org/downloads/INCHI-1-SRC.zip")
    endif()
    if(NOT DEFINED INCHI_MD5SUM)
      set(INCHI_MD5SUM "f2efa0c58cef32915686c04d7055b4e9")
    endif()
    if(NOT DEFINED INCHI_BASE)
      string(REGEX REPLACE "^.*/" "" INCHI_BASE "${INCHI_URL}")
    endif()
    downloadAndCheckMD5(${INCHI_URL} "${CUSTOM_INCHI_PATH}/${INCHI_BASE}" ${INCHI_MD5SUM})
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf
      ${CUSTOM_INCHI_PATH}/INCHI-1-SRC.zip
      WORKING_DIRECTORY ${CUSTOM_INCHI_PATH})
    execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_directory
        "${CUSTOM_INCHI_PATH}/INCHI-1-SRC" "${CUSTOM_INCHI_PATH}/src")
  endif(INCHI_FOUND)
endif(EXISTS ${CUSTOM_INCHI_PATH}/src/INCHI_BASE/src/ichican2.c)
if(EXISTS ${CUSTOM_INCHI_PATH}/src/INCHI_BASE/src/ichican2.c)
  set(INCHI_LIBRARIES Inchi)
endif(EXISTS ${CUSTOM_INCHI_PATH}/src/INCHI_BASE/src/ichican2.c)
if((NOT EXISTS ${CUSTOM_INCHI_PATH}/src/INCHI_BASE/src/ichican2.c)
  AND (NOT INCHI_FOUND))
  message(WARNING  "** NO INCHI SOFTWARE FOUND\n"
          "InChI support will be turned off. If you want to add InChI support, please follow the instructions inside $RDBASE/External/INCHI-API/README to download a copy of InChI software and then rerun cmake.")
  set(RDK_BUILD_INCHI_SUPPORT OFF)
endif()
