//
// Copyright (c) 2018 greg Landrum
//
//   @@ All Rights Reserved  @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// RDKit configuration options
#cmakedefine RDK_USE_BOOST_SERIALIZATION
#cmakedefine RDK_USE_BOOST_IOSTREAMS
#cmakedefine RDK_USE_BOOST_STACKTRACE

#cmakedefine RDK_OPTIMIZE_POPCNT

#cmakedefine RDK_BUILD_THREADSAFE_SSS

#cmakedefine RDK_TEST_MULTITHREADED

#cmakedefine RDK_USE_STRICT_ROTOR_DEFINITION

#cmakedefine RDK_BUILD_DESCRIPTORS3D
#ifdef RDK_BUILD_DESCRIPTORS3D
#define RDK_HAS_EIGEN3
#endif

#cmakedefine RDK_BUILD_COORDGEN_SUPPORT

#cmakedefine RDK_BUILD_MAEPARSER_SUPPORT

#cmakedefine RDK_BUILD_AVALON_SUPPORT

#cmakedefine RDK_BUILD_INCHI_SUPPORT

#cmakedefine RDK_BUILD_SLN_SUPPORT

#cmakedefine RDK_BUILD_CAIRO_SUPPORT

#cmakedefine RDK_BUILD_FREETYPE_SUPPORT

#cmakedefine RDK_USE_URF

#cmakedefine RDK_BUILD_YAEHMOP_SUPPORT

#cmakedefine RDK_BUILD_XYZ2MOL_SUPPORT
