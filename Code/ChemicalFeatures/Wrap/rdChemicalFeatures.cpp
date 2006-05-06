// $Id: rdChemicalFeatures.cpp 4944 2006-02-17 01:23:55Z glandrum $
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//     All Rights Reserved
//
//  @@ All Rights Reserved @@
//

#include <RDBoost/Wrap.h>
#include <ChemicalFeatures/FreeChemicalFeature.h>
#include <RDBoost/PySequenceHolder.h>

void wrap_freefeat();



BOOST_PYTHON_MODULE(rdChemicalFeatures)
{
  python::scope().attr("__doc__") =
    "Module containing free chemical feature functionality\n\
     These are feature that ar not derived from molecules. They are \n\
     are typically derived from pharmacophores and site-map.\n";

  wrap_freefeat();
}
