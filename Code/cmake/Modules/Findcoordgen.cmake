# Try to find Schrodinger's CoorgGen libraries.
#
# Different version handling is not yet supported
#
# Once found, this will find and define the following variables:
#
# coordgen_INCLUDE_DIRS   - CoordGen's includes directory
# coordgen_LIBRARIES      - CoordGen's shared libraries
# coordgen_TEMPLATE_FILE  - CoordGen templates file
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

# Just in case, add parent directory above libraries to templates search hints
get_filename_component(libs_parent_dir ${coordgen_LIBRARIES} PATH)
find_file(coordgen_TEMPLATE_FILE
    NAMES templates.mae
    HINTS ${COORDGEN_DIR} ${coordgen_DIR} ${libs_parent_dir}
    PATH_SUFFIXES "share" "share/coordgen"
    DOC "templates file for coordgen"
)
message("-- coordgen templates file set as '${coordgen_TEMPLATE_FILE}'")

find_package_handle_standard_args(coordgen FOUND_VAR coordgen_FOUND
                                  REQUIRED_VARS coordgen_INCLUDE_DIRS
                                  coordgen_LIBRARIES coordgen_TEMPLATE_FILE)



