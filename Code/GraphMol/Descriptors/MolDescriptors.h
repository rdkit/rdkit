//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#ifndef _RD_MOLDESCRIPTORS_H_
#define _RD_MOLDESCRIPTORS_H_

#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/AtomPairs.h>
#include <GraphMol/Descriptors/MolSurf.h>

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
    double CalcAMW(const ROMol &mol,bool onlyHeavy=false);



  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
