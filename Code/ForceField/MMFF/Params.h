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
#ifndef __RD_MMFFPARAMS_H__
#define __RD_MMFFPARAMS_H__

#include <memory>
#include <RDGeneral/Invariant.h>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <cstdint>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// binary searches are slightly faster than std::map;
// however when I moved to binary searches I had already
// written the code for std::map, so the two methods
// can be toggled defining RDKIT_MMFF_PARAMS_USE_STD_MAP

// #define RDKIT_MMFF_PARAMS_USE_STD_MAP 1

namespace ForceFields {
namespace MMFF {

const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;
const double MDYNE_A_TO_KCAL_MOL = 143.9325;
inline bool isDoubleZero(const double x) {
  return ((x < 1.0e-10) && (x > -1.0e-10));
}
inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }

//! class to store MMFF atom type equivalence levels
class RDKIT_FORCEFIELD_EXPORT MMFFDef {
 public:
  std::uint8_t eqLevel[4];
};

//! class to store MMFF Properties
class RDKIT_FORCEFIELD_EXPORT MMFFProp {
 public:
  std::uint8_t atno;
  std::uint8_t crd;
  std::uint8_t val;
  std::uint8_t pilp;
  std::uint8_t mltb;
  std::uint8_t arom;
  std::uint8_t linh;
  std::uint8_t sbmb;
};

//! class to store MMFF Partial Bond Charge Increments
class RDKIT_FORCEFIELD_EXPORT MMFFPBCI {
 public:
  double pbci;
  double fcadj;
};

//! class to store MMFF bond-charge-increment parameters used to
//! construct MMFF partial atomic charges
class RDKIT_FORCEFIELD_EXPORT MMFFChg {
 public:
  double bci;
};

//! class to store MMFF parameters for bond stretching
class RDKIT_FORCEFIELD_EXPORT MMFFBond {
 public:
  double kb;
  double r0;
};

//! class to store parameters for Herschbach-Laurie's version
//! of Badger's rule
class RDKIT_FORCEFIELD_EXPORT MMFFHerschbachLaurie {
 public:
  double a_ij;
  double d_ij;
  double dp_ij;
};

//! class to store covalent radius and Pauling electronegativity
//! values for MMFF bond stretching empirical rule
class RDKIT_FORCEFIELD_EXPORT MMFFCovRadPauEle {
 public:
  double r0;
  double chi;
};

//! class to store MMFF parameters for angle bending
class RDKIT_FORCEFIELD_EXPORT MMFFAngle {
 public:
  double ka;
  double theta0;
};

//! class to store MMFF parameters for stretch-bending
class RDKIT_FORCEFIELD_EXPORT MMFFStbn {
 public:
  double kbaIJK;
  double kbaKJI;
};

//! class to store MMFF parameters for out-of-plane bending
class RDKIT_FORCEFIELD_EXPORT MMFFOop {
 public:
  double koop;
};

//! class to store MMFF parameters for torsions
class RDKIT_FORCEFIELD_EXPORT MMFFTor {
 public:
  double V1;
  double V2;
  double V3;
};

//! class to store MMFF parameters for non-bonded Van der Waals
class RDKIT_FORCEFIELD_EXPORT MMFFVdW {
 public:
  double alpha_i;
  double N_i;
  double A_i;
  double G_i;
  double R_star;
  std::uint8_t DA;
};

class RDKIT_FORCEFIELD_EXPORT MMFFVdWRijstarEps {
 public:
  double R_ij_starUnscaled;
  double epsilonUnscaled;
  double R_ij_star;
  double epsilon;
};

class RDKIT_FORCEFIELD_EXPORT MMFFAromCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFArom object, NULL on failure.
  */
  bool isMMFFAromatic(const unsigned int atomType) const {
    return std::find(d_params.begin(), d_params.end(), atomType) !=
           d_params.end();
  }

  MMFFAromCollection(const std::vector<std::uint8_t> *mmffArom = nullptr);
  std::vector<std::uint8_t> d_params;  //!< the aromatic type vector
};

class RDKIT_FORCEFIELD_EXPORT MMFFDefCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFDef object, NULL on failure.
  */
  const MMFFDef *operator()(const unsigned int atomType) const {
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res = d_params.find(atomType);
    return ((res != d_params.end()) ? &((*res).second) : NULL);
#else
    return ((atomType && (atomType <= d_params.size()))
                ? &d_params[atomType - 1]
                : nullptr);
#endif
  }

  MMFFDefCollection(std::string mmffDef = "");

#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int, MMFFDef> d_params;  //!< the parameter map
#else
  std::vector<MMFFDef> d_params;  //!< the parameter vector
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFPropCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFProp object, NULL on failure.
  */
  const MMFFProp *operator()(const unsigned int atomType) const {
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res = d_params.find(atomType);
    return ((res != d_params.end()) ? &((*res).second) : NULL);
#else
    auto bounds =
        std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), atomType);
    return ((bounds.first != bounds.second)
                ? &d_params[bounds.first - d_iAtomType.begin()]
                : nullptr);
#endif
  }

  MMFFPropCollection(std::string mmffProp = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int, MMFFProp> d_params;  //!< the parameter map
#else
  std::vector<MMFFProp> d_params;
  std::vector<std::uint8_t> d_iAtomType;  //!< the parameter vector
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFPBCICollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFPBCI object, NULL on failure.
  */
  const MMFFPBCI *operator()(const unsigned int atomType) const {
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res = d_params.find(atomType);
    return ((res != d_params.end()) ? &((*res).second) : NULL);
#else
    return ((atomType && (atomType <= d_params.size()))
                ? &d_params[atomType - 1]
                : nullptr);
#endif
  }

  MMFFPBCICollection(std::string mmffPBCI = "");

#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int, MMFFPBCI> d_params;  //!< the parameter map
#else
  std::vector<MMFFPBCI> d_params;  //!< the parameter vector
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFChgCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFChg object, NULL on failure.
  */
  const std::pair<int, const MMFFChg *> getMMFFChgParams(
      const unsigned int bondType, const unsigned int iAtomType,
      const unsigned int jAtomType) const {
    int sign = -1;
    const MMFFChg *mmffChgParams = nullptr;
    unsigned int canIAtomType = iAtomType;
    unsigned int canJAtomType = jAtomType;
    if (iAtomType > jAtomType) {
      canIAtomType = jAtomType;
      canJAtomType = iAtomType;
      sign = 1;
    }
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res1 = d_params[bondType].find(canIAtomType);
    if (res1 != d_params[bondType].end()) {
      const auto res2 = ((*res1).second).find(canJAtomType);
      if (res2 != ((*res1).second).end()) {
        mmffChgParams = &((*res2).second);
      }
    }
#else
    auto bounds =
        std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canIAtomType);
    if (bounds.first != bounds.second) {
      bounds = std::equal_range(
          d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
          d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()),
          canJAtomType);
      if (bounds.first != bounds.second) {
        bounds = std::equal_range(
            d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
            d_bondType.begin() + (bounds.second - d_jAtomType.begin()),
            bondType);
        if (bounds.first != bounds.second) {
          mmffChgParams = &d_params[bounds.first - d_bondType.begin()];
        }
      }
    }
#endif

    return std::make_pair(sign, mmffChgParams);
  }

  MMFFChgCollection(std::string mmffChg = "");

//!< the parameter 3D-map
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int,
           std::map<const unsigned int, std::map<const unsigned int, MMFFChg>>>
      d_params;  //!< the parameter 3D-map
#else
  std::vector<MMFFChg> d_params;          //!< the parameter vector
  std::vector<std::uint8_t> d_iAtomType;  //!< atom type vector for atom i
  std::vector<std::uint8_t> d_jAtomType;  //!< atom type vector for atom j
  std::vector<std::uint8_t> d_bondType;   //!< bond type vector for bond i-j
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFBondCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFBond object, NULL on failure.
  */
  const MMFFBond *operator()(const unsigned int bondType,
                             const unsigned int atomType,
                             const unsigned int nbrAtomType) const {
    const MMFFBond *mmffBondParams = nullptr;
    unsigned int canAtomType = atomType;
    unsigned int canNbrAtomType = nbrAtomType;
    if (atomType > nbrAtomType) {
      canAtomType = nbrAtomType;
      canNbrAtomType = atomType;
    }
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res1 = d_params.find(bondType);
    std::map<const unsigned int,
             std::map<const unsigned int, MMFFBond>>::const_iterator res2;
    std::map<const unsigned int, MMFFBond>::const_iterator res3;
    if (res1 != d_params.end()) {
      res2 = ((*res1).second).find(canAtomType);
      if (res2 != ((*res1).second).end()) {
        res3 = ((*res2).second).find(canNbrAtomType);
        if (res3 != ((*res2).second).end()) {
          mmffBondParams = &((*res3).second);
        }
      }
    }
#else
    auto bounds =
        std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canAtomType);
    if (bounds.first != bounds.second) {
      bounds = std::equal_range(
          d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
          d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()),
          canNbrAtomType);
      if (bounds.first != bounds.second) {
        bounds = std::equal_range(
            d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
            d_bondType.begin() + (bounds.second - d_jAtomType.begin()),
            bondType);
        if (bounds.first != bounds.second) {
          mmffBondParams = &d_params[bounds.first - d_bondType.begin()];
        }
      }
    }
#endif

    return mmffBondParams;
  }

  MMFFBondCollection(std::string mmffBond = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int,
           std::map<const unsigned int, std::map<const unsigned int, MMFFBond>>>
      d_params;  //!< the parameter 3D-map
#else
  std::vector<MMFFBond> d_params;         //!< the parameter vector
  std::vector<std::uint8_t> d_iAtomType;  //!< atom type vector for atom i
  std::vector<std::uint8_t> d_jAtomType;  //!< atom type vector for atom j
  std::vector<std::uint8_t> d_bondType;   //!< bond type vector for bond i-j
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFBndkCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFBndk object, NULL on failure.
  */
  const MMFFBond *operator()(const int atomicNum,
                             const int nbrAtomicNum) const {
    const MMFFBond *mmffBndkParams = nullptr;
    unsigned int canAtomicNum = atomicNum;
    unsigned int canNbrAtomicNum = nbrAtomicNum;
    if (atomicNum > nbrAtomicNum) {
      canAtomicNum = nbrAtomicNum;
      canNbrAtomicNum = atomicNum;
    }
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res1 = d_params.find(canAtomicNum);
    std::map<const unsigned int, MMFFBond>::const_iterator res2;
    if (res1 != d_params.end()) {
      res2 = ((*res1).second).find(canNbrAtomicNum);
      if (res2 != ((*res1).second).end()) {
        mmffBndkParams = &((*res2).second);
      }
    }
#else
    auto bounds = std::equal_range(d_iAtomicNum.begin(), d_iAtomicNum.end(),
                                   canAtomicNum);
    if (bounds.first != bounds.second) {
      bounds = std::equal_range(
          d_jAtomicNum.begin() + (bounds.first - d_iAtomicNum.begin()),
          d_jAtomicNum.begin() + (bounds.second - d_iAtomicNum.begin()),
          canNbrAtomicNum);
      if (bounds.first != bounds.second) {
        mmffBndkParams = &d_params[bounds.first - d_jAtomicNum.begin()];
      }
    }
#endif

    return mmffBndkParams;
  }

  MMFFBndkCollection(std::string mmffBndk = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int, std::map<const unsigned int, MMFFBond>>
      d_params;  //!< the parameter 2D-map
#else
  std::vector<MMFFBond> d_params;          //!< the parameter vector
  std::vector<std::uint8_t> d_iAtomicNum;  //!< atomic number vector for atom i
  std::vector<std::uint8_t> d_jAtomicNum;  //!< atomic number vector for atom j
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFHerschbachLaurieCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFHerschbachLaurie object, NULL on failure.
  */
  const MMFFHerschbachLaurie *operator()(const int iRow, const int jRow) const {
    const MMFFHerschbachLaurie *mmffHerschbachLaurieParams = nullptr;
    unsigned int canIRow = iRow;
    unsigned int canJRow = jRow;
    if (iRow > jRow) {
      canIRow = jRow;
      canJRow = iRow;
    }
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res1 = d_params.find(canIRow);
    std::map<const unsigned int, MMFFHerschbachLaurie>::const_iterator res2;
    if (res1 != d_params.end()) {
      res2 = ((*res1).second).find(canJRow);
      if (res2 != ((*res1).second).end()) {
        mmffHerschbachLaurieParams = &((*res2).second);
      }
    }
#else
    auto bounds = std::equal_range(d_iRow.begin(), d_iRow.end(), canIRow);
    if (bounds.first != bounds.second) {
      bounds = std::equal_range(
          d_jRow.begin() + (bounds.first - d_iRow.begin()),
          d_jRow.begin() + (bounds.second - d_iRow.begin()), canJRow);
      if (bounds.first != bounds.second) {
        mmffHerschbachLaurieParams = &d_params[bounds.first - d_jRow.begin()];
      }
    }
#endif

    return mmffHerschbachLaurieParams;
  }

  MMFFHerschbachLaurieCollection(std::string mmffHerschbachLaurie = "");

#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int,
           std::map<const unsigned int, MMFFHerschbachLaurie>>
      d_params;  //!< the parameter 2D-map
#else
  std::vector<MMFFHerschbachLaurie> d_params;  //!< the parameter vector
  std::vector<std::uint8_t> d_iRow;  //!< periodic row number vector for atom i
  std::vector<std::uint8_t> d_jRow;  //!< periodic row number vector for atom j
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFCovRadPauEleCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFCovRadPauEle object, NULL on failure.
  */
  const MMFFCovRadPauEle *operator()(const unsigned int atomicNum) const {
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res = d_params.find(atomicNum);
    return ((res != d_params.end()) ? &((*res).second) : NULL);
#else
    auto bounds =
        std::equal_range(d_atomicNum.begin(), d_atomicNum.end(), atomicNum);
    return ((bounds.first != bounds.second)
                ? &d_params[bounds.first - d_atomicNum.begin()]
                : nullptr);
#endif
  }

  MMFFCovRadPauEleCollection(std::string mmffCovRadPauEle = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int, MMFFCovRadPauEle>
      d_params;  //!< the parameter map
#else
  std::vector<MMFFCovRadPauEle> d_params;  //!< the parameter vector
  std::vector<std::uint8_t> d_atomicNum;   //!< the atomic number vector
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFAngleCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFAngle object, NULL on failure.
  */
  const MMFFAngle *operator()(const MMFFDefCollection *mmffDef,
                              const unsigned int angleType,
                              const unsigned int iAtomType,
                              const unsigned int jAtomType,
                              const unsigned int kAtomType) const {
    const MMFFAngle *mmffAngleParams = nullptr;
    unsigned int iter = 0;

// For bending of the i-j-k angle, a five-stage process based
// in the level combinations 1-1-1,2-2-2,3-2-3,4-2-4, and
// 5-2-5 is used. (MMFF.I, note 68, page 519)
// We skip 1-1-1 since Level 2 === Level 1
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    while ((iter < 4) && (!mmffAngleParams)) {
      unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iter];
      unsigned int canKAtomType = (*mmffDef)(kAtomType)->eqLevel[iter];
      if (canIAtomType > canKAtomType) {
        std::swap(canIAtomType, canKAtomType);
      }
      const auto res1 = d_params.find(angleType);
      if (res1 != d_params.end()) {
        const auto res2 = ((*res1).second).find(canIAtomType);
        if (res2 != ((*res1).second).end()) {
          const auto res3 = ((*res2).second).find(jAtomType);
          if (res3 != ((*res2).second).end()) {
            const auto res4 = ((*res3).second).find(canKAtomType);
            if (res4 != ((*res3).second).end()) {
              mmffAngleParams = &((*res4).second);
            }
          }
        }
      }
      ++iter;
    }
#else
    auto jBounds =
        std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
    if (jBounds.first != jBounds.second) {
      while ((iter < 4) && (!mmffAngleParams)) {
        unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iter];
        unsigned int canKAtomType = (*mmffDef)(kAtomType)->eqLevel[iter];
        if (canIAtomType > canKAtomType) {
          std::swap(canIAtomType, canKAtomType);
        }

        auto bounds = std::equal_range(
            d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
            d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
            canIAtomType);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range(
              d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
              d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
              canKAtomType);
          if (bounds.first != bounds.second) {
            bounds = std::equal_range(
                d_angleType.begin() + (bounds.first - d_kAtomType.begin()),
                d_angleType.begin() + (bounds.second - d_kAtomType.begin()),
                angleType);
            if (bounds.first != bounds.second) {
              mmffAngleParams = &d_params[bounds.first - d_angleType.begin()];
            }
          }
        }
        ++iter;
      }
    }
#endif

    return mmffAngleParams;
  }

  MMFFAngleCollection(std::string mmffAngle = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int,
           std::map<const unsigned int,
                    std::map<const unsigned int,
                             std::map<const unsigned int, MMFFAngle>>>>
      d_params;  //!< the parameter 4D-map
#else
  std::vector<MMFFAngle> d_params;        //!< the parameter vector
  std::vector<std::uint8_t> d_iAtomType;  //!< atom type vector for atom i
  std::vector<std::uint8_t> d_jAtomType;  //!< atom type vector for atom j
  std::vector<std::uint8_t> d_kAtomType;  //!< atom type vector for atom k
  std::vector<std::uint8_t> d_angleType;  //!< angle type vector for angle i-j-k
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFStbnCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFStbn object, NULL on failure.
  */
  const std::pair<bool, const MMFFStbn *> getMMFFStbnParams(
      const unsigned int stretchBendType, const unsigned int bondType1,
      const unsigned int bondType2, const unsigned int iAtomType,
      const unsigned int jAtomType, const unsigned int kAtomType) const {
    const MMFFStbn *mmffStbnParams = nullptr;
    bool swap = false;
    unsigned int canIAtomType = iAtomType;
    unsigned int canKAtomType = kAtomType;
    unsigned int canStretchBendType = stretchBendType;
    if (iAtomType > kAtomType) {
      canIAtomType = kAtomType;
      canKAtomType = iAtomType;
      swap = true;
    } else if (iAtomType == kAtomType) {
      swap = (bondType1 < bondType2);
    }
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res1 = d_params.find(canStretchBendType);
    if (res1 != d_params.end()) {
      const auto res2 = ((*res1).second).find(canIAtomType);
      if (res2 != ((*res1).second).end()) {
        const auto res3 = ((*res2).second).find(jAtomType);
        if (res3 != ((*res2).second).end()) {
          const auto res4 = ((*res3).second).find(canKAtomType);
          if (res4 != ((*res3).second).end()) {
            mmffStbnParams = &((*res4).second);
          }
        }
      }
    }
#else
    auto jBounds =
        std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
    if (jBounds.first != jBounds.second) {
      auto bounds = std::equal_range(
          d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
          d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
          canIAtomType);
      if (bounds.first != bounds.second) {
        bounds = std::equal_range(
            d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
            d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
            canKAtomType);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range(
              d_stretchBendType.begin() + (bounds.first - d_kAtomType.begin()),
              d_stretchBendType.begin() + (bounds.second - d_kAtomType.begin()),
              canStretchBendType);
          if (bounds.first != bounds.second) {
            mmffStbnParams =
                &d_params[bounds.first - d_stretchBendType.begin()];
          }
        }
      }
    }
#endif

    return std::make_pair(swap, mmffStbnParams);
  }

  MMFFStbnCollection(std::string mmffStbn = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int,
           std::map<const unsigned int,
                    std::map<const unsigned int,
                             std::map<const unsigned int, MMFFStbn>>>>
      d_params;  //!< the parameter 4D-map
#else
  std::vector<MMFFStbn> d_params;         //!< the parameter vector
  std::vector<std::uint8_t> d_iAtomType;  //!< atom type vector for atom i
  std::vector<std::uint8_t> d_jAtomType;  //!< atom type vector for atom j
  std::vector<std::uint8_t> d_kAtomType;  //!< atom type vector for atom k
  std::vector<std::uint8_t>
      d_stretchBendType;  //!< stretch-bend type vector for angle i-j-k
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFDfsbCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFStbn object, NULL on failure.
  */
  const std::pair<bool, const MMFFStbn *> getMMFFDfsbParams(
      const unsigned int periodicTableRow1,
      const unsigned int periodicTableRow2,
      const unsigned int periodicTableRow3) const {
    const MMFFStbn *mmffDfsbParams = nullptr;
    bool swap = false;
    unsigned int canPeriodicTableRow1 = periodicTableRow1;
    unsigned int canPeriodicTableRow3 = periodicTableRow3;
    if (periodicTableRow1 > periodicTableRow3) {
      canPeriodicTableRow1 = periodicTableRow3;
      canPeriodicTableRow3 = periodicTableRow1;
      swap = true;
    }
    const auto res1 = d_params.find(canPeriodicTableRow1);
    if (res1 != d_params.end()) {
      const auto res2 = ((*res1).second).find(periodicTableRow2);
      if (res2 != ((*res1).second).end()) {
        const auto res3 = ((*res2).second).find(canPeriodicTableRow3);
        if (res3 != ((*res2).second).end()) {
          mmffDfsbParams = &((*res3).second);
        }
      }
    }

    return std::make_pair(swap, mmffDfsbParams);
  }

  MMFFDfsbCollection(std::string mmffDfsb = "");
  std::map<const unsigned int,
           std::map<const unsigned int, std::map<const unsigned int, MMFFStbn>>>
      d_params;  //!< the parameter 3D-map
};

class RDKIT_FORCEFIELD_EXPORT MMFFOopCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFOop object, NULL on failure.
  */
  const MMFFOop *operator()(const MMFFDefCollection *mmffDef,
                            const unsigned int iAtomType,
                            const unsigned int jAtomType,
                            const unsigned int kAtomType,
                            const unsigned int lAtomType) const {
    const MMFFOop *mmffOopParams = nullptr;
    unsigned int iter = 0;
    std::vector<unsigned int> canIKLAtomType(3);
// For out-of-plane bending ijk; I , where j is the central
// atom [cf. eq. (511, the five-stage protocol 1-1-1; 1, 2-2-2; 2,
// 3-2-3;3, 4-2-4;4, 5-2-5;5 is used. The final stage provides
// wild-card defaults for all except the central atom.
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    while ((iter < 4) && (!mmffOopParams)) {
      canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
      unsigned int canJAtomType = jAtomType;
      canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
      canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
      std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
      const auto res1 = d_params.find(canIKLAtomType[0]);
      if (res1 != d_params.end()) {
        const auto res2 = ((*res1).second).find(canJAtomType);
        if (res2 != ((*res1).second).end()) {
          const auto res3 = ((*res2).second).find(canIKLAtomType[1]);
          if (res3 != ((*res2).second).end()) {
            const auto res4 = ((*res3).second).find(canIKLAtomType[2]);
            if (res4 != ((*res3).second).end()) {
              mmffOopParams = &((*res4).second);
            }
          }
        }
      }
      ++iter;
    }
#else
    auto jBounds =
        std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
    if (jBounds.first != jBounds.second) {
      while ((iter < 4) && (!mmffOopParams)) {
        canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
        canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
        canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
        std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
        auto bounds = std::equal_range(
            d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
            d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
            canIKLAtomType[0]);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range(
              d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
              d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
              canIKLAtomType[1]);
          if (bounds.first != bounds.second) {
            bounds = std::equal_range(
                d_lAtomType.begin() + (bounds.first - d_kAtomType.begin()),
                d_lAtomType.begin() + (bounds.second - d_kAtomType.begin()),
                canIKLAtomType[2]);
            if (bounds.first != bounds.second) {
              mmffOopParams = &d_params[bounds.first - d_lAtomType.begin()];
            }
          }
        }
        ++iter;
      }
    }
#endif

    return mmffOopParams;
  }

  MMFFOopCollection(const bool isMMFFs, std::string mmffOop = "");

#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int,
           std::map<const unsigned int,
                    std::map<const unsigned int,
                             std::map<const unsigned int, MMFFOop>>>>
      d_params;  //!< the parameter 4D-map
#else
  std::vector<MMFFOop> d_params;          //!< the parameter vector
  std::vector<std::uint8_t> d_iAtomType;  //!< atom type vector for atom i
  std::vector<std::uint8_t> d_jAtomType;  //!< atom type vector for atom j
  std::vector<std::uint8_t> d_kAtomType;  //!< atom type vector for atom k
  std::vector<std::uint8_t> d_lAtomType;  //!< atom type vector for atom l
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFTorCollection {
 public:
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFTor object, NULL on failure.
  */
  const std::pair<const unsigned int, const MMFFTor *> getMMFFTorParams(
      const MMFFDefCollection *mmffDef,
      const std::pair<unsigned int, unsigned int> torType,
      const unsigned int iAtomType, const unsigned int jAtomType,
      const unsigned int kAtomType, const unsigned int lAtomType) const {
    const MMFFTor *mmffTorParams = nullptr;
    unsigned int iter = 0;
    unsigned int iWildCard = 0;
    unsigned int lWildCard = 0;
    unsigned int canTorType = torType.first;
    unsigned int maxIter = 5;
// For i-j-k-2 torsion interactions, a five-stage
// process based on level combinations 1-1-1-1, 2-2-2-2,
// 3-2-2-5, 5-2-2-3, and 5-2-2-5 is used, where stages 3
// and 4 correspond to "half-default" or "half-wild-card" entries.
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
#else
#endif

    while (((iter < maxIter) && ((!mmffTorParams) || (maxIter == 4))) ||
           ((iter == 4) && (torType.first == 5) && torType.second)) {
      // The rule of setting the torsion type to the value it had
      // before being set to 5 as a last resort in case parameters
      // could not be found is not mentioned in MMFF.IV; it was
      // empirically discovered due to a number of tests in the
      // MMFF validation suite otherwise failing
      if ((maxIter == 5) && (iter == 4)) {
        maxIter = 4;
        iter = 0;
        canTorType = torType.second;
      }
      iWildCard = iter;
      lWildCard = iter;
      if (iter == 1) {
        iWildCard = 1;
        lWildCard = 3;
      } else if (iter == 2) {
        iWildCard = 3;
        lWildCard = 1;
      }
      unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iWildCard];
      unsigned int canJAtomType = jAtomType;
      unsigned int canKAtomType = kAtomType;
      unsigned int canLAtomType = (*mmffDef)(lAtomType)->eqLevel[lWildCard];
      if (canJAtomType > canKAtomType) {
        unsigned int temp = canKAtomType;
        canKAtomType = canJAtomType;
        canJAtomType = temp;
        temp = canLAtomType;
        canLAtomType = canIAtomType;
        canIAtomType = temp;
      } else if ((canJAtomType == canKAtomType) &&
                 (canIAtomType > canLAtomType)) {
        unsigned int temp = canLAtomType;
        canLAtomType = canIAtomType;
        canIAtomType = temp;
      }
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      const auto res1 = d_params.find(canTorType);
      if (res1 != d_params.end()) {
        const auto res2 = ((*res1).second).find(canIAtomType);
        if (res2 != ((*res1).second).end()) {
          const auto res3 = ((*res2).second).find(canJAtomType);
          if (res3 != ((*res2).second).end()) {
            const auto res4 = ((*res3).second).find(canKAtomType);
            if (res4 != ((*res3).second).end()) {
              const auto res5 = ((*res4).second).find(canLAtomType);
              if (res5 != ((*res4).second).end()) {
                mmffTorParams = &((*res5).second);
                if (maxIter == 4) {
                  break;
                }
              }
            }
          }
        }
      }
#else
      auto jBounds = std::equal_range(d_jAtomType.begin(), d_jAtomType.end(),
                                      canJAtomType);
      if (jBounds.first != jBounds.second) {
        auto bounds = std::equal_range(
            d_kAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
            d_kAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
            canKAtomType);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range(
              d_iAtomType.begin() + (bounds.first - d_kAtomType.begin()),
              d_iAtomType.begin() + (bounds.second - d_kAtomType.begin()),
              canIAtomType);
          if (bounds.first != bounds.second) {
            bounds = std::equal_range(
                d_lAtomType.begin() + (bounds.first - d_iAtomType.begin()),
                d_lAtomType.begin() + (bounds.second - d_iAtomType.begin()),
                canLAtomType);
            if (bounds.first != bounds.second) {
              bounds = std::equal_range(
                  d_torType.begin() + (bounds.first - d_lAtomType.begin()),
                  d_torType.begin() + (bounds.second - d_lAtomType.begin()),
                  canTorType);
              if (bounds.first != bounds.second) {
                mmffTorParams = &d_params[bounds.first - d_torType.begin()];
                if (maxIter == 4) {
                  break;
                }
              }
            }
          }
        }
      }
#endif
      ++iter;
    }

    return std::make_pair(canTorType, mmffTorParams);
  }

  MMFFTorCollection(const bool isMMFFs, std::string mmffTor = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<
      const unsigned int,
      std::map<
          const unsigned int,
          std::map<const unsigned int,
                   std::map<const unsigned int, std::map<const unsigned int,
                                                         MMFFTor>>>>>
      d_params;  //!< the parameter 5D-map
#else
  std::vector<MMFFTor> d_params;          //!< the parameter vector
  std::vector<std::uint8_t> d_iAtomType;  //!< atom type vector for atom i
  std::vector<std::uint8_t> d_jAtomType;  //!< atom type vector for atom j
  std::vector<std::uint8_t> d_kAtomType;  //!< atom type vector for atom k
  std::vector<std::uint8_t> d_lAtomType;  //!< atom type vector for atom l
  std::vector<std::uint8_t>
      d_torType;  //!< torsion type vector for angle i-j-k-l
#endif
};

class RDKIT_FORCEFIELD_EXPORT MMFFVdWCollection {
 public:
  //! gets a pointer to the singleton MMFFVdWCollection
  double power;
  double B;
  double Beta;
  double DARAD;
  double DAEPS;
  //! Looks up the parameters for a particular key and returns them.
  /*!
    \return a pointer to the MMFFVdW object, NULL on failure.
  */
  const MMFFVdW *operator()(const unsigned int atomType) const {
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
    const auto res = d_params.find(atomType);
    return (res != d_params.end() ? &((*res).second) : NULL);
#else
    auto bounds =
        std::equal_range(d_atomType.begin(), d_atomType.end(), atomType);
    return ((bounds.first != bounds.second)
                ? &d_params[bounds.first - d_atomType.begin()]
                : nullptr);
#endif
  }

  MMFFVdWCollection(std::string mmffVdW = "");
#ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
  std::map<const unsigned int, MMFFVdW> d_params;  //!< the parameter map
#else
  std::vector<MMFFVdW> d_params;         //!< the parameter vector
  std::vector<std::uint8_t> d_atomType;  //!< atom type vector
#endif
};
}  // namespace MMFF
}  // namespace ForceFields

#endif
