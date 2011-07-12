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
  find_package_handle_standard_args(INCHI
	  "Could NOT find InChI in system locations" INCHI_LIBRARY
	  INCHI_INCLUDE_DIR)
  set(INCHI_INCLUDE_DIRS ${INCHI_INCLUDE_DIR} )
  set(INCHI_LIBRARIES ${INCHI_LIBRARY} )
  mark_as_advanced(INCHI_INCLUDE_DIR INCHI_LIBRARY)
endif()
