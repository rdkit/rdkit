#include <RDGeneral/versions.h>

const char * RDKit::rdkitVersion = "@RDKit_RELEASENAME@";

// The Boost version as detected at build time.
// CMake's Boost_LIB_VERSION is defined by the FindBoost.cmake module
// to be the same as the value from <boost/version.hpp>
const char * RDKit::boostVersion = "@Boost_LIB_VERSION@";
