# Try to find Schrodinger's CoorgGen libraries.
#
# Different version handling is not yet supported
#
# Once found, this will find and define the following variables:
#
# coordgen_INCLUDE_DIRS   - CoordGen's includes directory
# coordgen_LIBRARIES      - CoordGen's shared libraries
#
#

include(FindPackageHandleStandardArgs)

find_path(coordgen_INCLUDE_DIRS
    NAMES "coordgen/sketcherMinimizer.h"
    HINTS ${COORDGEN_DIR} ${coordgen_DIR}
    PATH_SUFFIXES "include"
    DOC "include path for coordgen"
)
message("-- coordgen include dir set as ${coordgen_INCLUDE_DIRS}")

find_library(coordgen_LIBRARIES
    NAMES coordgen coordgenlibs
    HINTS ${COORDGEN_DIR} ${coordgen_DIR}
    PATH_SUFFIXES "lib"
    DOC "libraries for coordgen"
)
message("-- coordgen libraries set as '${coordgen_LIBRARIES}'")

find_package_handle_standard_args(coordgen FOUND_VAR coordgen_FOUND
                                  REQUIRED_VARS coordgen_INCLUDE_DIRS
                                  coordgen_LIBRARIES)



