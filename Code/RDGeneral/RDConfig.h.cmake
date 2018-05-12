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

#cmakedefine RDK_OPTIMIZE_NATIVE
#ifdef RDK_OPTIMIZE_NATIVE
#define USE_BUILTIN_POPCOUNT
#endif

#cmakedefine RDK_BUILD_THREADSAFE_SSS
#ifdef RDK_BUILD_THREADSAFE_SSS
#define RDK_THREADSAFE_SSS
#endif

#cmakedefine RDK_TEST_MULTITHREADED

#cmakedefine RDK_USE_STRICT_ROTOR_DEFINITION

#cmakedefine RDK_BUILD_DESCRIPTORS3D
#ifdef RDK_BUILD_DESCRIPTORS3D
#define RDK_HAS_EIGEN3
#endif

#cmakedefine RDK_BUILD_COORDGEN_SUPPORT

#cmakedefine RDK_BUILD_AVALON_SUPPORT

#cmakedefine RDK_BUILD_INCHI_SUPPORT

#cmakedefine RDK_BUILD_SLN_SUPPORT

#cmakedefine RDK_BUILD_CAIRO_SUPPORT
