// $Id: versions.h 1443 2010-07-01 03:55:11Z glandrum $
//
// Copyright (c) 2010 greg Landrum
//
//   @@ All Rights Reserved  @@
//
namespace RDKit {
  const char * const rdkitVersion="2016.03.1.dev1";

  // The Boost version as detected at build time.
  // CMake's Boost_LIB_VERSION is defined by the FindBoost.cmake module
  // to be the same as the value from <boost/version.hpp>
  const char * const boostVersion="1_59";
}  
