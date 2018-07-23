//
//  Copyright (C) 2013 Paolo Tosco
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
#ifndef __RD_UFFPOSITIONCONSTRAINT_H__
#define __RD_UFFPOSITIONCONSTRAINT_H__
#include <iostream>
#include <ForceField/Contrib.h>
#include <Geometry/point.h>

namespace ForceFields {
namespace UFF {

//! A position constraint of the type 0.5k * deltaX^2
class RDKIT_FORCEFIELD_EXPORT PositionConstraintContrib : public ForceFieldContrib {
 public:
  PositionConstraintContrib() : d_atIdx(-1){};
  //! Constructor
  /*!
    \param owner       pointer to the owning ForceField
    \param idx         index of the atom in the ForceField's positions
    \param minDispl    minimum displacement
    \param maxDispl    maximum displacement
    \param forceConst  force constant

  */
  PositionConstraintContrib(ForceField *owner, unsigned int idx,
                            double maxDispl, double forceConst);

  ~PositionConstraintContrib() {}
  double getEnergy(double *pos) const;

  void getGrad(double *pos, double *grad) const;
  virtual PositionConstraintContrib *copy() const {
    return new PositionConstraintContrib(*this);
  };

 private:
  int d_atIdx;             //!< index of the restrained atom
  double d_maxDispl;       //!< maximum allowed displacement
  RDGeom::Point3D d_pos0;  //!< reference position
  double d_forceConstant;  //!< force constant of the bond
};
}
}
#endif
