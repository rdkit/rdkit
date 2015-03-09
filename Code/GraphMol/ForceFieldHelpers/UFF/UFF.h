//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_UFFCONVENIENCE_H
#define RD_UFFCONVENIENCE_H
#include <ForceField/ForceField.h>

#include "Builder.h"

namespace RDKit {
  class ROMol;
  namespace UFF {
    //! Convenience function for optimizing a molecule using UFF
    /*
      \param mol        the molecule to use
      \param maxIters   the maximum number of force-field iterations
      \param vdwThresh  the threshold to be used in adding van der Waals terms
                        to the force field. Any non-bonded contact whose current
			distance is greater than \c vdwThresh * the minimum value
			for that contact will not be included.
      \param confId     the optional conformer id, if this isn't provided, the molecule's
                        default confId will be used.
      \param ignoreInterfragInteractions if true, nonbonded terms will not be added between
                                         fragments

      \return 0 if the optimization converged, 1 if more iterations are required.
    */
    std::pair<int,double> UFFOptimizeMolecule(ROMol &mol, int maxIters=1000,
                            double vdwThresh=10.0, int confId=-1,
                            bool ignoreInterfragInteractions=true ){
      ForceFields::ForceField *ff=UFF::constructForceField(mol,vdwThresh, confId,
                                                           ignoreInterfragInteractions);
      ff->initialize();
      int res=ff->minimize(maxIters);
      int e=ff->calcEnergy();
      delete ff;
      return std::make_pair(res,e);
    }
  }
}
#endif
