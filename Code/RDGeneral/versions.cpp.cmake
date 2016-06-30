#include <RDGeneral/versions.h>

const char * RDKit::rdkitVersion = "@RDKit_RELEASENAME@";

// The Boost version as detected at build time.
// CMake's Boost_LIB_VERSION is defined by the FindBoost.cmake module
// to be the same as the value from <boost/version.hpp>
const char * RDKit::boostVersion = "@Boost_LIB_VERSION@";

// The system/compiler on which RDKit was built as detected at build time.
const char * RDKit::rdkitBuild = "@RDKit_BUILDNAME@";
