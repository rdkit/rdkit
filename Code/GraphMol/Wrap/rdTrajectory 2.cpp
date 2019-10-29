// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>
#include <RDBoost/python.h>

namespace python = boost::python;
void wrap_trajectory();

BOOST_PYTHON_MODULE(rdtrajectory) {
  python::scope().attr("__doc__") =
      "Module containing Trajectory and Snapshot objects\n";

  wrap_trajectory();
}
