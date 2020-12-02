# Grabbed from here:
# https://github.com/avaucher/cmake-find-module-imported-target/blob/master/cmake/FindCairo.cmake

# - Try to find Cairo
# Once done, this will define
#
#  CAIRO_FOUND - system has Cairo
#  CAIRO_INCLUDE_DIRS - the Cairo include directories
#  CAIRO_LIBRARIES - link these to use Cairo
#
# and the following imported targets::
#
#  Cairo::Cairo - The Cairo library
#
#
# Original work Copyright (C) 2012 Raphael Kubo da Costa <rakuco@webkit.org>
# Modified work Copyright (C) 2017 Alain Vaucher
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1.  Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 2.  Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND ITS CONTRIBUTORS ``AS
# IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

FIND_PACKAGE(PkgConfig)
PKG_CHECK_MODULES(PC_CAIRO cairo QUIET)

FIND_PATH(CAIRO_INCLUDE_DIRS
    NAMES cairo.h
    HINTS ${CAIRO_INCLUDE_DIR}
          ${PC_CAIRO_INCLUDEDIR}
          ${PC_CAIRO_INCLUDE_DIRS}
    PATH_SUFFIXES cairo
)

FIND_LIBRARY(CAIRO_LIBRARIES
    NAMES cairo
    HINTS ${CAIRO_LIBRARY_DIR}
          ${PC_CAIRO_LIBDIR}
          ${PC_CAIRO_LIBRARY_DIRS}
)

IF (CAIRO_INCLUDE_DIRS)
    IF (EXISTS "${CAIRO_INCLUDE_DIRS}/cairo-version.h")
        FILE(READ "${CAIRO_INCLUDE_DIRS}/cairo-version.h" CAIRO_VERSION_CONTENT)

        STRING(REGEX MATCH "#define +CAIRO_VERSION_MAJOR +([0-9]+)" _dummy "${CAIRO_VERSION_CONTENT}")
        SET(CAIRO_VERSION_MAJOR "${CMAKE_MATCH_1}")

        STRING(REGEX MATCH "#define +CAIRO_VERSION_MINOR +([0-9]+)" _dummy "${CAIRO_VERSION_CONTENT}")
        SET(CAIRO_VERSION_MINOR "${CMAKE_MATCH_1}")

        STRING(REGEX MATCH "#define +CAIRO_VERSION_MICRO +([0-9]+)" _dummy "${CAIRO_VERSION_CONTENT}")
        SET(CAIRO_VERSION_MICRO "${CMAKE_MATCH_1}")

        SET(CAIRO_VERSION "${CAIRO_VERSION_MAJOR}.${CAIRO_VERSION_MINOR}.${CAIRO_VERSION_MICRO}")
    ENDIF ()
ENDIF ()

SET(VERSION_OK TRUE)
IF (Cairo_FIND_VERSION)
    IF (Cairo_FIND_VERSION_EXACT)
        IF (NOT "${Cairo_FIND_VERSION}" VERSION_EQUAL "${CAIRO_VERSION}")
            SET(VERSION_OK FALSE)
        ENDIF ()
    ELSE ()
        IF ("${Cairo_FIND_VERSION}" VERSION_GREATER "${CAIRO_VERSION}")
            SET(VERSION_OK FALSE)
        ENDIF ()
    ENDIF ()
ENDIF ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cairo
                                  FOUND_VAR Cairo_FOUND
                                  REQUIRED_VARS CAIRO_INCLUDE_DIRS CAIRO_LIBRARIES VERSION_OK
                                  )

if(Cairo_FOUND AND NOT TARGET Cairo::Cairo)
  add_library(Cairo::Cairo INTERFACE IMPORTED)
  set_target_properties(Cairo::Cairo PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES  "${CAIRO_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES       "${CAIRO_LIBRARIES}")
endif()
