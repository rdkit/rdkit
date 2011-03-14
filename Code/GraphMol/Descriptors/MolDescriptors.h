//
//  Copyright (C) 2004-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef _RD_MOLDESCRIPTORS_H_
#define _RD_MOLDESCRIPTORS_H_

#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/Descriptors/Lipinski.h>

namespace RDKit{
  class ROMol;
  namespace Descriptors {
    /*!
      Calculates a molecule's molecular weight

      \param mol        the molecule of interest
      \param onlyHeavy  (optional) if this is true (the default is false),
          only heavy atoms will be included in the MW calculation

      \return the AMW
    */
    double calcAMW(const ROMol &mol,bool onlyHeavy=false);

  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
