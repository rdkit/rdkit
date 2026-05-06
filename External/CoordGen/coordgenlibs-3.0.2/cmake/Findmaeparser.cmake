# Try to find Schrodinger's MAEParser libraries.
#
# Different version handling is not yet supported
#
# Once found, this will find and define the following variables:
#
# maeparser_INCLUDE_DIRS  - maeparser's includes directory
# maeparser_LIBRARIES     - maeparser's shared libraries
#
#

include(FindPackageHandleStandardArgs)

find_path(maeparser_INCLUDE_DIRS
    NAMES "maeparser/Reader.hpp"
    HINTS ${maeparser_INCLUDE_DIRS} ${maeparser_DIR}
    PATH_SUFFIXES "include"
    DOC "include path for maeparser"
)
message(STATUS "maeparser include dir set as '${maeparser_INCLUDE_DIRS}'")

find_library(maeparser_LIBRARIES
    NAMES maeparser
    HINTS ${maeparser_LIBRARIES} ${maeparser_DIR}
    PATH_SUFFIXES "lib"
    DOC "libraries for maeparser"
)
message(STATUS "maeparser libraries set as '${maeparser_LIBRARIES}'")

find_package_handle_standard_args(maeparser FOUND_VAR maeparser_FOUND
                                  REQUIRED_VARS maeparser_INCLUDE_DIRS
                                  maeparser_LIBRARIES)
