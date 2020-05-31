//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_FFCONTRIB_H__
#define __RD_FFCONTRIB_H__

namespace ForceFields {
class ForceField;

//! abstract base class for contributions to ForceFields
class RDKIT_FORCEFIELD_EXPORT ForceFieldContrib {
 public:
  friend class ForceField;

  ForceFieldContrib()  {};
  ForceFieldContrib(ForceFields::ForceField *owner) : dp_forceField(owner){};
  virtual ~ForceFieldContrib(){};

  //! returns our contribution to the energy of a position
  virtual double getEnergy(double *pos) const = 0;

  //! calculates our contribution to the gradients of a position
  virtual void getGrad(double *pos, double *grad) const = 0;

  //! return a copy
  virtual ForceFieldContrib *copy() const = 0;

 protected:
  ForceField *dp_forceField{nullptr};  //!< our owning ForceField
};
}  // namespace ForceFields

#endif
