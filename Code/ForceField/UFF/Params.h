//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_UFFPARAMS_H
#define RD_UFFPARAMS_H

#include <memory>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ForceFields {
namespace UFF {

const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;
inline bool isDoubleZero(const double x) {
  return ((x < 1.0e-10) && (x > -1.0e-10));
}
inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }

//! class to store UFF parameters for bond stretching
class RDKIT_FORCEFIELD_EXPORT UFFBond {
 public:
  double kb;
  double r0;
};

//! class to store UFF parameters for angle bending
class RDKIT_FORCEFIELD_EXPORT UFFAngle {
 public:
  double ka;
  double theta0;
};

//! class to store UFF parameters for torsions
class RDKIT_FORCEFIELD_EXPORT UFFTor {
 public:
  double V;
};

//! class to store UFF parameters for inversions
class RDKIT_FORCEFIELD_EXPORT UFFInv {
 public:
  double K;
};

//! class to store UFF parameters for van der Waals interactions
class RDKIT_FORCEFIELD_EXPORT UFFVdW {
 public:
  double x_ij;
  double D_ij;
};

//! class to store atomic parameters for the Universal Force Field
class RDKIT_FORCEFIELD_EXPORT AtomicParams {
 public:
  double r1;            //!<  valence bond radius
  double theta0;        //!< valence angle
  double x1;            //!< vdW characteristic length
  double D1;            //!< vdW atomic energy
  double zeta;          //!< vdW scaling term
  double Z1;            //!< effective charge
  double V1;            //!< sp3 torsional barrier parameter
  double U1;            //!< torsional contribution for sp2-sp3 bonds
  double GMP_Xi;        //!< GMP Electronegativity;
  double GMP_Hardness;  //!< GMP Hardness
  double GMP_Radius;    //!< GMP Radius value
};

namespace Params {
const double lambda = 0.1332;  //!< scaling factor for rBO correction
const double G = 332.06;       //!< bond force constant prefactor
const double amideBondOrder =
    1.41;  //!< special case bond order for amide C-N bonds.
};         // namespace Params

//! singleton class for retrieving UFF AtomParams
/*!
  Use the singleton like this:

  \verbatim
  ParamCollection *params=ParamCollection::getParams();
  const AtomParams *ap=params("C_3");
  \endverbatim

  If you have your own parameter data, it can be supplied as a string:
  \verbatim
  ParamCollection *params=ParamCollection::getParams(myParamData);
  const AtomParams *ap=params("C_3");
  \endverbatim
  You are responsible for making sure that the data is in the correct
  format (see Params.cpp for an example).

*/
class RDKIT_FORCEFIELD_EXPORT ParamCollection {
 public:
  //! gets a pointer to the singleton ParamCollection
  /*!
    \param paramData (optional) a string with parameter data. See
     below for more information about this argument

    \return a pointer to the singleton ParamCollection

    <b>Notes:</b>
      - do <b>not</b> delete the pointer returned here
      - if the singleton ParamCollection has already been instantiated and
        \c paramData is empty, the singleton will be returned.
      - if \c paramData is empty and the singleton ParamCollection has
        not yet been instantiated, the default UFF parameters (from Params.cpp)
        will be used.
      - if \c paramData is supplied, a new singleton will be instantiated.
        The current instantiation (if there is one) will be deleted.
  */
  static const ParamCollection *getParams(const std::string &paramData = "");
  //! Looks up the parameters for a particular UFF key and returns them.
  /*!
    \return a pointer to the AtomicParams object, NULL on failure.
  */
  const AtomicParams *operator()(const std::string &symbol) const {
    std::map<std::string, AtomicParams>::const_iterator res;
    res = d_params.find(symbol);
    if (res != d_params.end()) {
      return &((*res).second);
    }
    return nullptr;
  }

  ParamCollection(std::string paramData);

 private:
  std::map<std::string, AtomicParams> d_params;  //!< the parameter map
};
}  // namespace UFF
}  // namespace ForceFields
#endif
