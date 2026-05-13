// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/Wrap.h>
#include <ChemicalFeatures/FreeChemicalFeature.h>
#include <RDBoost/PySequenceHolder.h>

void wrap_freefeat();

BOOST_PYTHON_MODULE(rdChemicalFeatures) {
  python::scope().attr("__doc__") =
      "Module containing free chemical feature functionality\n\
     These are features that are not associated with molecules. They are \n\
     typically derived from pharmacophores and site-maps.\n";

  wrap_freefeat();
}
