//
//  Copyright (C) 2015 Sereina Riniker
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_TORSIONANGLEM6_H
#define RD_TORSIONANGLEM6_H

#include <ForceField/Contrib.h>
#include <boost/tuple/tuple.hpp>
#include <vector>

namespace RDGeom {
class Point3D;
}

namespace ForceFields {
class ForceField;
class ForceFieldContrib;
}

namespace ForceFields {
namespace CrystalFF {

//! the torsion term for multiplicity m = 1 - 6
class TorsionAngleContribM6 : public ForceFields::ForceFieldContrib {
 public:
  TorsionAngleContribM6()
      : d_at1Idx(-1), d_at2Idx(-1), d_at3Idx(-1), d_at4Idx(-1){};
  //! Constructor
  /*!
   The torsion is between atom1 - atom2 - atom3 - atom4
   (i.e the angle between bond atom1-atom2 and bond atom3-atom4
   while looking down bond atom2-atom3)

   \param owner       pointer to the owning ForceField
   \param idx1        index of atom1 in the ForceField's positions
   \param idx2        index of atom2 in the ForceField's positions
   \param idx3        index of atom3 in the ForceField's positions
   \param idx4        index of atom4 in the ForceField's positions
   \param V           list of 6 force constants
   \param signs       list of 6 signs (+1 or -1)
  */
  TorsionAngleContribM6(ForceFields::ForceField *owner, unsigned int idx1,
                        unsigned int idx2, unsigned int idx3, unsigned int idx4,
                        const std::vector<double> &V,
                        const std::vector<int> &signs);
  double getEnergy(double *pos) const;
  void getGrad(double *pos, double *grad) const;
  virtual TorsionAngleContribM6 *copy() const {
    return new TorsionAngleContribM6(*this);
  };

 private:
  int d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx;
  std::vector<double> d_V;
  std::vector<int> d_sign;
};

//! calculates and returns the torsional energy
double calcTorsionEnergyM6(const std::vector<double> &V,
                           const std::vector<int> &signs, const double cosPhi);
}
}
#endif
