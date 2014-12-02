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
#ifndef __RD_MMFFPARAMS_H__
#define __RD_MMFFPARAMS_H__

#include <RDGeneral/Invariant.h>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <boost/cstdint.hpp>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

// binary searches are slightly faster than std::map;
// however when I moved to binary searches I had already
// written the code for std::map, so the two methods
// can be toggled defining RDKIT_MMFF_PARAMS_USE_STD_MAP

//#define RDKIT_MMFF_PARAMS_USE_STD_MAP 1

namespace ForceFields {
  namespace MMFF {

    const double DEG2RAD = M_PI / 180.0;
    const double RAD2DEG = 180.0 / M_PI;
    inline bool isDoubleZero(const double x) {
      return ((x < 1.0e-10) && (x > -1.0e-10));
    }
    inline void clipToOne(double &x) {
      if (x > 1.0) {
        x = 1.0;
      }
      else if (x < -1.0) {
        x = -1.0;
      }
    }

    //! class to store MMFF atom type equivalence levels
    class MMFFDef {
    public:
      boost::uint8_t eqLevel[4];
    };

    //! class to store MMFF Properties
    class MMFFProp {
    public:
      boost::uint8_t atno;
      boost::uint8_t crd;
      boost::uint8_t val;
      boost::uint8_t pilp;
      boost::uint8_t mltb;
      boost::uint8_t arom;
      boost::uint8_t linh;
      boost::uint8_t sbmb;
    };

    //! class to store MMFF Partial Bond Charge Increments
    class MMFFPBCI {
    public:
      double pbci;
      double fcadj;
    };

    //! class to store MMFF bond-charge-increment parameters used to
    //! construct MMFF partial atomic charges
    class MMFFChg {
    public:
      double bci;
    };

    //! class to store MMFF parameters for bond stretching
    class MMFFBond {
    public:
      double kb;
      double r0;
    };

    //! class to store parameters for Herschbach-Laurie's version
    //! of Badger's rule
    class MMFFHerschbachLaurie {
    public:
      double a_ij;
      double d_ij;
      double dp_ij;
    };

    //! class to store covalent radius and Pauling electronegativity
    //! values for MMFF bond stretching empirical rule
    class MMFFCovRadPauEle {
    public:
      double r0;
      double chi;
    };

    //! class to store MMFF parameters for angle bending
    class MMFFAngle {
    public:
      double ka;
      double theta0;
    };

    //! class to store MMFF parameters for stretch-bending
    class MMFFStbn {
    public:
      double kbaIJK;
      double kbaKJI;
    };

    //! class to store MMFF parameters for out-of-plane bending
    class MMFFOop {
    public:
      double koop;
    };

    //! class to store MMFF parameters for torsions
    class MMFFTor {
    public:
      double V1;
      double V2;
      double V3;
    };

    //! class to store MMFF parameters for non-bonded Van der Waals
    class MMFFVdW {
    public:
      double alpha_i;
      double N_i;
      double A_i;
      double G_i;
      double R_star;
      boost::uint8_t DA;
    };

    class MMFFAromCollection {
    public:
      //! gets a pointer to the singleton MMFFAromCollection
      /*!
	\param mmffArom (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFAromCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFAromCollection has already been instantiated and
	    \c mmffArom is empty, the singleton will be returned.
	  - if \c mmffArom is empty and the singleton MMFFAromCollection has
	    not yet been instantiated, the default MMFFArom parameters (from Params.cpp)
	    will be used.
	  - if \c mmffArom is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFAromCollection *getMMFFArom(const boost::uint8_t *aromatic_types = NULL);
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFArom object, NULL on failure.
      */
      bool isMMFFAromatic(const unsigned int atomType) const {
        return ((std::find(d_params.begin(), d_params.end(),
          atomType) != d_params.end()) ? true : false);
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFAromCollection(const boost::uint8_t mmffArom[]);
      static class MMFFAromCollection *ds_instance;    //!< the singleton
      std::vector<boost::uint8_t> d_params;  //!< the aromatic type vector
    };

    class MMFFDefCollection {
    public:
      //! gets a pointer to the singleton MMFFDefCollection
      /*!
	\param mmffDef (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFDefCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFDefCollection has already been instantiated and
	    \c mmffDef is empty, the singleton will be returned.
	  - if \c mmffDef is empty and the singleton MMFFDefCollection has
	    not yet been instantiated, the default MMFFDef parameters (from Params.cpp)
	    will be used.
	  - if \c mmffDef is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFDefCollection *getMMFFDef(const std::string &mmffDef="");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFDef object, NULL on failure.
      */
      const MMFFDef *operator()(const unsigned int atomType) const {
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, MMFFDef>::const_iterator res;
        res = d_params.find(atomType);

        return ((res != d_params.end()) ? &((*res).second) : NULL);
        #else
        return ((atomType && (atomType <= d_params.size()))
          ? &d_params[atomType - 1] : NULL);
        #endif
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFDefCollection(std::string mmffDef);
      static class MMFFDefCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, MMFFDef> d_params;  //!< the parameter map
      #else
      std::vector<MMFFDef> d_params;  //!< the parameter vector
      #endif
    };

    class MMFFPropCollection {
    public:
      //! gets a pointer to the singleton MMFFPropCollection
      /*!
	\param mmffProp (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFPropCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFPropCollection has already been instantiated and
	    \c mmffProp is empty, the singleton will be returned.
	  - if \c mmffProp is empty and the singleton MMFFPropCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffProp is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFPropCollection *getMMFFProp(const std::string &mmffProp="");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFProp object, NULL on failure.
      */
      const MMFFProp *operator()(const unsigned int atomType) const {
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, MMFFProp>::const_iterator res;
        res = d_params.find(atomType);

        return ((res != d_params.end()) ? &((*res).second) : NULL);
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds =
          std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), atomType);

        return ((bounds.first != bounds.second)
          ? &d_params[bounds.first - d_iAtomType.begin()] : NULL);
        #endif
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFPropCollection(std::string mmffProp);
      static class MMFFPropCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, MMFFProp> d_params;  //!< the parameter map
      #else
      std::vector<MMFFProp> d_params;
      std::vector<boost::uint8_t> d_iAtomType;  //!< the parameter vector
      #endif
    };

    class MMFFPBCICollection {
    public:
      //! gets a pointer to the singleton MMFFPBCICollection
      /*!
	\param mmffPBCI (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFPBCICollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFPBCICollection has already been instantiated and
	    \c mmffPBCI is empty, the singleton will be returned.
	  - if \c mmffPBCI is empty and the singleton MMFFPBCICollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffPBCI is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFPBCICollection *getMMFFPBCI(const std::string &mmffPBCI = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFPBCI object, NULL on failure.
      */
      const MMFFPBCI *operator()(const unsigned int atomType) const {
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, MMFFPBCI>::const_iterator res;
        res = d_params.find(atomType);

        return ((res != d_params.end()) ? &((*res).second) : NULL);
        #else
        return ((atomType && (atomType <= d_params.size()))
          ? &d_params[atomType - 1] : NULL);
        #endif
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFPBCICollection(std::string mmffPBCI);
      static class MMFFPBCICollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, MMFFPBCI> d_params;  //!< the parameter map
      #else
      std::vector<MMFFPBCI> d_params;  //!< the parameter vector
      #endif
    };

    class MMFFChgCollection {
    public:
      //! gets a pointer to the singleton MMFFChgCollection
      /*!
	\param mmffChg (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFChgCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFChgCollection has already been instantiated and
	    \c mmffChg is empty, the singleton will be returned.
	  - if \c mmffChg is empty and the singleton MMFFChgCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffChg is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFChgCollection *getMMFFChg(const std::string &mmffChg = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFChg object, NULL on failure.
      */
      const std::pair<int, const MMFFChg *> getMMFFChgParams(const unsigned int bondType,
        const unsigned int iAtomType, const unsigned int jAtomType) {

        int sign = -1;
        const MMFFChg *mmffChgParams = NULL;
        unsigned int canIAtomType = iAtomType;
        unsigned int canJAtomType = jAtomType;
        if (iAtomType > jAtomType) {
          canIAtomType = jAtomType;
          canJAtomType = iAtomType;
          sign = 1;
        }
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, std::map<const unsigned int, MMFFChg> >::const_iterator res1;
        std::map<const unsigned int, MMFFChg>::const_iterator res2;
        res1 = d_params[bondType].find(canIAtomType);
        if (res1 != d_params[bondType].end()) {
          res2 = ((*res1).second).find(canJAtomType);
          if (res2 != ((*res1).second).end()) {
            mmffChgParams = &((*res2).second);
          }
        }
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        
        bounds = std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canIAtomType);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range(d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
            d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()), canJAtomType);
          if (bounds.first != bounds.second) {
            bounds = std::equal_range
              (d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
              d_bondType.begin() + (bounds.second - d_jAtomType.begin()), bondType);
            if (bounds.first != bounds.second) {
              mmffChgParams = &d_params[bounds.first - d_bondType.begin()];
            }
          }
        }
        #endif

        return std::make_pair(sign, mmffChgParams);
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFChgCollection(std::string mmffChg);
      static class MMFFChgCollection *ds_instance;    //!< the singleton
      //!< the parameter 3D-map
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, std::map<const unsigned int,
        std::map<const unsigned int, MMFFChg> > > d_params;  //!< the parameter 3D-map
      #else
      std::vector<MMFFChg> d_params;  //! the parameter vector
      std::vector<boost::uint8_t> d_iAtomType;  //! atom type vector for atom i
      std::vector<boost::uint8_t> d_jAtomType;  //! atom type vector for atom j
      std::vector<boost::uint8_t> d_bondType;  //! bond type vector for bond i-j
      #endif
    };

    class MMFFBondCollection {
    public:
      //! gets a pointer to the singleton MMFFBondCollection
      /*!
	\param mmffBond (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFBondCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFBondCollection has already been instantiated and
	    \c mmffBond is empty, the singleton will be returned.
	  - if \c mmffBond is empty and the singleton MMFFBondCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffBond is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFBondCollection *getMMFFBond(const std::string &mmffBond = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFBond object, NULL on failure.
      */
      const MMFFBond *operator()(const unsigned int bondType,
        const unsigned int atomType, const unsigned int nbrAtomType) {
        
        const MMFFBond *mmffBondParams = NULL;
        unsigned int canAtomType = atomType;
        unsigned int canNbrAtomType = nbrAtomType;
        if (atomType > nbrAtomType) {
          canAtomType = nbrAtomType;
          canNbrAtomType = atomType;
        }
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, std::map<const unsigned int,
            std::map<const unsigned int, MMFFBond> > >::const_iterator res1;
        std::map<const unsigned int,
            std::map<const unsigned int, MMFFBond> >::const_iterator res2;
        std::map<const unsigned int, MMFFBond>::const_iterator res3;
        res1 = d_params.find(bondType);
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
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        bounds = std::equal_range(d_iAtomType.begin(), d_iAtomType.end(), canAtomType);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range(d_jAtomType.begin() + (bounds.first - d_iAtomType.begin()),
            d_jAtomType.begin() + (bounds.second - d_iAtomType.begin()), canNbrAtomType);
          if (bounds.first != bounds.second) {
              bounds = std::equal_range
                (d_bondType.begin() + (bounds.first - d_jAtomType.begin()),
                d_bondType.begin() + (bounds.second - d_jAtomType.begin()), bondType);
            if (bounds.first != bounds.second) {
              mmffBondParams = &d_params[bounds.first - d_bondType.begin()];
            }
          }
        }
        #endif

        return mmffBondParams;
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFBondCollection(std::string mmffBond);
      static class MMFFBondCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, std::map<const unsigned int,
        std::map<const unsigned int, MMFFBond> > > d_params;  //!< the parameter 3D-map
      #else
      std::vector<MMFFBond> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_iAtomType;  //! atom type vector for atom i
      std::vector<boost::uint8_t> d_jAtomType;  //! atom type vector for atom j
      std::vector<boost::uint8_t> d_bondType;  //! bond type vector for bond i-j
      #endif
    };

    class MMFFBndkCollection {
    public:
      //! gets a pointer to the singleton MMFFBndkCollection
      /*!
	\param mmffBndk (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFBndkCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFBndkCollection has already been instantiated and
	    \c mmffBndk is empty, the singleton will be returned.
	  - if \c mmffBndk is empty and the singleton MMFFBndkCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffBndk is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFBndkCollection *getMMFFBndk(const std::string &mmffBndk = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFBndk object, NULL on failure.
      */
      const MMFFBond *operator()(const int atomicNum, const int nbrAtomicNum) {
        
        const MMFFBond *mmffBndkParams = NULL;
        unsigned int canAtomicNum = atomicNum;
        unsigned int canNbrAtomicNum = nbrAtomicNum;
        if (atomicNum > nbrAtomicNum) {
          canAtomicNum = nbrAtomicNum;
          canNbrAtomicNum = atomicNum;
        }
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int,
          std::map<const unsigned int, MMFFBond> >::const_iterator res1;
        std::map<const unsigned int, MMFFBond>::const_iterator res2;
        res1 = d_params.find(canAtomicNum);
        if (res1 != d_params.end()) {
          res2 = ((*res1).second).find(canNbrAtomicNum);
          if (res2 != ((*res1).second).end()) {
            mmffBndkParams = &((*res2).second);
          }
        }
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        bounds = std::equal_range
          (d_iAtomicNum.begin(), d_iAtomicNum.end(), canAtomicNum);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range
            (d_jAtomicNum.begin() + (bounds.first - d_iAtomicNum.begin()),
            d_jAtomicNum.begin() + (bounds.second - d_iAtomicNum.begin()),
            canNbrAtomicNum);
          if (bounds.first != bounds.second) {
            mmffBndkParams = &d_params[bounds.first - d_jAtomicNum.begin()];
          }
        }
        #endif

        return mmffBndkParams;
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFBndkCollection(std::string mmffBndk);
      static class MMFFBndkCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, std::map<const unsigned int, MMFFBond> > d_params;  //!< the parameter 2D-map
      #else
      std::vector<MMFFBond> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_iAtomicNum;  //! atomic number vector for atom i
      std::vector<boost::uint8_t> d_jAtomicNum;  //! atomic number vector for atom j
      #endif
    };

    class MMFFHerschbachLaurieCollection {
    public:
      //! gets a pointer to the singleton MMFFHerschbachLaurieCollection
      /*!
	\param mmffHerschbachLaurie (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFHerschbachLaurieCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFHerschbachLaurieCollection has already been instantiated and
	    \c mmffHerschbachLaurie is empty, the singleton will be returned.
	  - if \c mmffHerschbachLaurie is empty and the singleton MMFFHerschbachLaurieCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffHerschbachLaurie is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFHerschbachLaurieCollection *getMMFFHerschbachLaurie(const std::string &mmffHerschbachLaurie = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFHerschbachLaurie object, NULL on failure.
      */
      const MMFFHerschbachLaurie *operator()(const int iRow, const int jRow) {
        
        const MMFFHerschbachLaurie *mmffHerschbachLaurieParams = NULL;
        unsigned int canIRow = iRow;
        unsigned int canJRow = jRow;
        if (iRow > jRow) {
          canIRow = jRow;
          canJRow = iRow;
        }
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int,
          std::map<const unsigned int, MMFFHerschbachLaurie> >::const_iterator res1;
        std::map<const unsigned int, MMFFHerschbachLaurie>::const_iterator res2;
        res1 = d_params.find(canIRow);
        if (res1 != d_params.end()) {
          res2 = ((*res1).second).find(canJRow);
          if (res2 != ((*res1).second).end()) {
            mmffHerschbachLaurieParams = &((*res2).second);
          }
        }
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        bounds = std::equal_range
          (d_iRow.begin(), d_iRow.end(), canIRow);
        if (bounds.first != bounds.second) {
          bounds = std::equal_range
            (d_jRow.begin() + (bounds.first - d_iRow.begin()),
            d_jRow.begin() + (bounds.second - d_iRow.begin()),
            canJRow);
          if (bounds.first != bounds.second) {
            mmffHerschbachLaurieParams = &d_params[bounds.first - d_jRow.begin()];
          }
        }
        #endif

        return mmffHerschbachLaurieParams;
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFHerschbachLaurieCollection(std::string mmffHerschbachLaurie);
      static class MMFFHerschbachLaurieCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, std::map<const unsigned int, MMFFHerschbachLaurie> > d_params;  //!< the parameter 2D-map
      #else
      std::vector<MMFFHerschbachLaurie> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_iRow;  //! periodic row number vector for atom i
      std::vector<boost::uint8_t> d_jRow;  //! periodic row number vector for atom j
      #endif
    };

    class MMFFCovRadPauEleCollection {
    public:
      //! gets a pointer to the singleton MMFFCovRadPauEleCollection
      /*!
	\param mmffCovRadPauEle (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFCovRadPauEleCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFCovRadPauEleCollection has already been instantiated and
	    \c mmffCovRadPauEle is empty, the singleton will be returned.
	  - if \c mmffCovRadPauEle is empty and the singleton MMFFCovRadPauEleCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffCovRadPauEle is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFCovRadPauEleCollection *getMMFFCovRadPauEle(const std::string &mmffCovRadPauEle = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFCovRadPauEle object, NULL on failure.
      */
      const MMFFCovRadPauEle *operator()(const unsigned int atomicNum) const {
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, MMFFCovRadPauEle>::const_iterator res;
        res = d_params.find(atomicNum);

        return ((res != d_params.end()) ? &((*res).second) : NULL);
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds =
          std::equal_range(d_atomicNum.begin(), d_atomicNum.end(), atomicNum);

        return ((bounds.first != bounds.second)
          ? &d_params[bounds.first - d_atomicNum.begin()] : NULL);
        #endif
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFCovRadPauEleCollection(std::string mmffCovRadPauEle);
      static class MMFFCovRadPauEleCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, MMFFCovRadPauEle> d_params;  //!< the parameter map
      #else
      std::vector<MMFFCovRadPauEle> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_atomicNum;  //!< the atomic number vector
      #endif
    };

    class MMFFAngleCollection {
    public:
      //! gets a pointer to the singleton MMFFAngleCollection
      /*!
	\param mmffAngle (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFAngleCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFAngleCollection has already been instantiated and
	    \c mmffAngle is empty, the singleton will be returned.
	  - if \c mmffAngle is empty and the singleton MMFFAngleCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffAngle is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFAngleCollection *getMMFFAngle(const std::string &mmffAngle = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFAngle object, NULL on failure.
      */
      const MMFFAngle *operator()(const unsigned int angleType,
        const unsigned int iAtomType, const unsigned int jAtomType,
        const unsigned int kAtomType ) {

        MMFFDefCollection *mmffDef = MMFFDefCollection::getMMFFDef();
        const MMFFAngle *mmffAngleParams = NULL;
        unsigned int iter = 0;

        // For bending of the i-j-k angle, a five-stage process based
        // in the level combinations 1-1-1,2-2-2,3-2-3,4-2-4, and
        // 5-2-5 is used. (MMFF.I, note 68, page 519)
        // We skip 1-1-1 since Level 2 === Level 1
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, std::map<const unsigned int, std::map<const unsigned int,
            std::map<const unsigned int, MMFFAngle> > > >::const_iterator res1;
        std::map<const unsigned int, std::map<const unsigned int, std::map<const unsigned int,
            MMFFAngle> > >::const_iterator res2;
        std::map<const unsigned int, std::map<const unsigned int, MMFFAngle> >::const_iterator res3;
        std::map<const unsigned int, MMFFAngle>::const_iterator res4;
        while ((iter < 4) && (!mmffAngleParams)) {
          unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iter];
          unsigned int canKAtomType = (*mmffDef)(kAtomType)->eqLevel[iter];
          if (canIAtomType > canKAtomType) {
            unsigned int temp = canKAtomType;
            canKAtomType = canIAtomType;
            canIAtomType = temp;
          }
          res1 = d_params.find(angleType);
          if (res1 != d_params.end()) {
            res2 = ((*res1).second).find(canIAtomType);
            if (res2 != ((*res1).second).end()) {
              res3 = ((*res2).second).find(jAtomType);
              if (res3 != ((*res2).second).end()) {
                res4 = ((*res3).second).find(canKAtomType);
                if (res4 != ((*res3).second).end()) {
                  mmffAngleParams = &((*res4).second);
                }
              }
            }
          }
          ++iter;
        }
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
            std::vector<boost::uint8_t>::const_iterator> jBounds =
          std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        if (jBounds.first != jBounds.second) {
          while ((iter < 4) && (!mmffAngleParams)) {
            unsigned int canIAtomType = (*mmffDef)(iAtomType)->eqLevel[iter];
            unsigned int canKAtomType = (*mmffDef)(kAtomType)->eqLevel[iter];
            if (canIAtomType > canKAtomType) {
              unsigned int temp = canKAtomType;
              canKAtomType = canIAtomType;
              canIAtomType = temp;
            }
            bounds = std::equal_range
              (d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
              d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
              canIAtomType);
            if (bounds.first != bounds.second) {
              bounds = std::equal_range
                (d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
                d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
                canKAtomType);
              if (bounds.first != bounds.second) {
                bounds = std::equal_range
                  (d_angleType.begin() + (bounds.first - d_kAtomType.begin()),
                  d_angleType.begin() + (bounds.second - d_kAtomType.begin()), angleType);
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
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFAngleCollection(std::string mmffAngle);
      static class MMFFAngleCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, std::map<const unsigned int,
        std::map<const unsigned int, std::map<const unsigned int, MMFFAngle> > > > d_params;  //!< the parameter 4D-map
      #else
      std::vector<MMFFAngle> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_iAtomType;  //! atom type vector for atom i
      std::vector<boost::uint8_t> d_jAtomType;  //! atom type vector for atom j
      std::vector<boost::uint8_t> d_kAtomType;  //! atom type vector for atom k
      std::vector<boost::uint8_t> d_angleType;  //! angle type vector for angle i-j-k
      #endif
    };

    class MMFFStbnCollection {
    public:
      //! gets a pointer to the singleton MMFFStbnCollection
      /*!
	\param mmffStbn (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFStbnCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFStbnCollection has already been instantiated and
	    \c mmffStbn is empty, the singleton will be returned.
	  - if \c mmffStbn is empty and the singleton MMFFStbnCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffStbn is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFStbnCollection *getMMFFStbn(const std::string &mmffStbn = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFStbn object, NULL on failure.
      */
      const std::pair<bool, const MMFFStbn *> getMMFFStbnParams
        (const unsigned int stretchBendType, const unsigned int bondType1,
        const unsigned int bondType2, const unsigned int iAtomType,
        const unsigned int jAtomType, const unsigned int kAtomType) {

        const MMFFStbn *mmffStbnParams = NULL;
        bool swap = false;
        unsigned int canIAtomType = iAtomType;
        unsigned int canKAtomType = kAtomType;
        unsigned int canStretchBendType = stretchBendType;
        if (iAtomType > kAtomType) {
          canIAtomType = kAtomType;
          canKAtomType = iAtomType;
          swap = true;
        }
        else if (iAtomType == kAtomType) {
          swap = (bondType1 < bondType2);
        }
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, std::map<const unsigned int, std::map<const unsigned int,
            std::map<const unsigned int, MMFFStbn> > > >::const_iterator res1;
        std::map<const unsigned int, std::map<const unsigned int, std::map<const unsigned int,
            MMFFStbn> > >::const_iterator res2;
        std::map<const unsigned int, std::map<const unsigned int, MMFFStbn> >::const_iterator res3;
        std::map<const unsigned int, MMFFStbn>::const_iterator res4;
        res1 = d_params.find(canStretchBendType);
        if (res1 != d_params.end()) {
          res2 = ((*res1).second).find(canIAtomType);
          if (res2 != ((*res1).second).end()) {
            res3 = ((*res2).second).find(jAtomType);
            if (res3 != ((*res2).second).end()) {
              res4 = ((*res3).second).find(canKAtomType);
              if (res4 != ((*res3).second).end()) {
                mmffStbnParams = &((*res4).second);
              }
            }
          }
        }
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
            std::vector<boost::uint8_t>::const_iterator> jBounds =
          std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        if (jBounds.first != jBounds.second) {
          bounds = std::equal_range
            (d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
            d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
            canIAtomType);
          if (bounds.first != bounds.second) {
            bounds = std::equal_range
              (d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
              d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
              canKAtomType);
            if (bounds.first != bounds.second) {
              bounds = std::equal_range
                (d_stretchBendType.begin() + (bounds.first - d_kAtomType.begin()),
                d_stretchBendType.begin() + (bounds.second - d_kAtomType.begin()),
                canStretchBendType);
              if (bounds.first != bounds.second) {
                mmffStbnParams = &d_params[bounds.first - d_stretchBendType.begin()];
              }
            }
          }
        }
        #endif
        
        return std::make_pair(swap, mmffStbnParams);
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFStbnCollection(std::string mmffStbn);
      static class MMFFStbnCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, std::map<const unsigned int,
        std::map<const unsigned int, std::map<const unsigned int, MMFFStbn> > > > d_params;  //!< the parameter 4D-map
      #else
      std::vector<MMFFStbn> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_iAtomType;  //! atom type vector for atom i
      std::vector<boost::uint8_t> d_jAtomType;  //! atom type vector for atom j
      std::vector<boost::uint8_t> d_kAtomType;  //! atom type vector for atom k
      std::vector<boost::uint8_t> d_stretchBendType;  //! stretch-bend type vector for angle i-j-k
      #endif
    };

    class MMFFDfsbCollection {
    public:
      //! gets a pointer to the singleton MMFFDfsbCollection
      /*!
	\param mmffDfsb (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFDfsbCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFDfsbCollection has already been instantiated and
	    \c mmffDfsb is empty, the singleton will be returned.
	  - if \c mmffDfsb is empty and the singleton MMFFDfsbCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffDfsb is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFDfsbCollection *getMMFFDfsb(const std::string &mmffDfsb = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFStbn object, NULL on failure.
      */
      const std::pair<bool, const MMFFStbn *> getMMFFDfsbParams(const unsigned int periodicTableRow1,
        const unsigned int periodicTableRow2, const unsigned int periodicTableRow3) {

        std::map<const unsigned int, std::map<const unsigned int, std::map<const unsigned int,
            MMFFStbn> > >::const_iterator res1;
        std::map<const unsigned int, std::map<const unsigned int, MMFFStbn> >::const_iterator res2;
        std::map<const unsigned int, MMFFStbn>::const_iterator res3;
        const MMFFStbn *mmffDfsbParams = NULL;
        bool swap = false;
        unsigned int canPeriodicTableRow1 = periodicTableRow1;
        unsigned int canPeriodicTableRow3 = periodicTableRow3;
        if (periodicTableRow1 > periodicTableRow3) {
          canPeriodicTableRow1 = periodicTableRow3;
          canPeriodicTableRow3 = periodicTableRow1;
          swap = true;
        }
        res1 = d_params.find(canPeriodicTableRow1);
        if (res1 != d_params.end()) {
          res2 = ((*res1).second).find(periodicTableRow2);
          if (res2 != ((*res1).second).end()) {
            res3 = ((*res2).second).find(canPeriodicTableRow3);
            if (res3 != ((*res2).second).end()) {
              mmffDfsbParams = &((*res3).second);
            }
          }
        }
        
        return std::make_pair(swap, mmffDfsbParams);
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFDfsbCollection(std::string mmffDfsb);
      static class MMFFDfsbCollection *ds_instance;    //!< the singleton
      std::map<const unsigned int, std::map<const unsigned int,
        std::map<const unsigned int, MMFFStbn> > > d_params;  //!< the parameter 3D-map
    };

    class MMFFOopCollection {
    public:
      //! gets a pointer to the singleton MMFFOopCollection
      /*!
	\param mmffOop (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFOopCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFOopCollection has already been instantiated and
	    \c mmffOop is empty, the singleton will be returned.
	  - if \c mmffOop is empty and the singleton MMFFOopCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffOop is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFOopCollection *getMMFFOop
        (const bool isMMFFs = false, const std::string &mmffOop = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFOop object, NULL on failure.
      */
      const MMFFOop *operator()(const unsigned int iAtomType, const unsigned int jAtomType,
        const unsigned int kAtomType, const unsigned int lAtomType) {

        MMFFDefCollection *mmffDef = MMFFDefCollection::getMMFFDef();
        const MMFFOop *mmffOopParams = NULL;
        unsigned int iter = 0;
        std::vector<unsigned int> canIKLAtomType(3);
        // For out-of-plane bending ijk; I , where j is the central
        // atom [cf. eq. (511, the five-stage protocol 1-1-1; 1, 2-2-2; 2,
        // 3-2-3;3, 4-2-4;4, 5-2-5;5 is used. The final stage provides
        // wild-card defaults for all except the central atom.
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, std::map<const unsigned int, std::map<const unsigned int,
            std::map<const unsigned int, MMFFOop> > > >::const_iterator res1;
        std::map<const unsigned int, std::map<const unsigned int, std::map<const unsigned int,
            MMFFOop> > >::const_iterator res2;
        std::map<const unsigned int, std::map<const unsigned int, MMFFOop> >::const_iterator res3;
        std::map<const unsigned int, MMFFOop>::const_iterator res4;
        while ((iter < 4) && (!mmffOopParams)) {
          canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
          unsigned int canJAtomType = jAtomType;
          canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
          canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
          std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
          res1 = d_params.find(canIKLAtomType[0]);
          if (res1 != d_params.end()) {
            res2 = ((*res1).second).find(canJAtomType);
            if (res2 != ((*res1).second).end()) {
              res3 = ((*res2).second).find(canIKLAtomType[1]);
              if (res3 != ((*res2).second).end()) {
                res4 = ((*res3).second).find(canIKLAtomType[2]);
                if (res4 != ((*res3).second).end()) {
                  mmffOopParams = &((*res4).second);
                }
              }
            }
          }
          ++iter;
        }
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> jBounds =
          std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), jAtomType);
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        if (jBounds.first != jBounds.second) {
          while ((iter < 4) && (!mmffOopParams)) {
            canIKLAtomType[0] = (*mmffDef)(iAtomType)->eqLevel[iter];
            canIKLAtomType[1] = (*mmffDef)(kAtomType)->eqLevel[iter];
            canIKLAtomType[2] = (*mmffDef)(lAtomType)->eqLevel[iter];
            std::sort(canIKLAtomType.begin(), canIKLAtomType.end());
            bounds = std::equal_range
              (d_iAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
              d_iAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
              canIKLAtomType[0]);
            if (bounds.first != bounds.second) {
              bounds = std::equal_range
                (d_kAtomType.begin() + (bounds.first - d_iAtomType.begin()),
                d_kAtomType.begin() + (bounds.second - d_iAtomType.begin()),
                canIKLAtomType[1]);
              if (bounds.first != bounds.second) {
                bounds = std::equal_range
                  (d_lAtomType.begin() + (bounds.first - d_kAtomType.begin()),
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
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFOopCollection(const bool isMMFFs, std::string mmffOop);
      static class MMFFOopCollection *ds_instance[2];    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, std::map<const unsigned int,
        std::map<const unsigned int, std::map<const unsigned int, MMFFOop> > > > d_params;  //!< the parameter 4D-map
      #else
      std::vector<MMFFOop> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_iAtomType;  //! atom type vector for atom i
      std::vector<boost::uint8_t> d_jAtomType;  //! atom type vector for atom j
      std::vector<boost::uint8_t> d_kAtomType;  //! atom type vector for atom k
      std::vector<boost::uint8_t> d_lAtomType;  //! atom type vector for atom l
      #endif
    };

    class MMFFTorCollection {
    public:
      //! gets a pointer to the singleton MMFFTorCollection
      /*!
	\param mmffTor (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFTorCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFTorCollection has already been instantiated and
	    \c mmffTor is empty, the singleton will be returned.
	  - if \c mmffTor is empty and the singleton MMFFTorCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffTor is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      static MMFFTorCollection *getMMFFTor(const bool isMMFFs, const std::string &mmffTor = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFTor object, NULL on failure.
      */
      const std::pair<const unsigned int, const MMFFTor *> getMMFFTorParams
        (const std::pair<unsigned int, unsigned int> torType,
        const unsigned int iAtomType, const unsigned int jAtomType,
        const unsigned int kAtomType, const unsigned int lAtomType)
      {
        MMFFDefCollection *mmffDef = MMFFDefCollection::getMMFFDef();
        const MMFFTor *mmffTorParams = NULL;
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
        std::map<const unsigned int,
          std::map<const unsigned int, std::map<const unsigned int,
          std::map<const unsigned int, std::map<const unsigned int,
          MMFFTor> > > > >::const_iterator res1;
        std::map<const unsigned int,
          std::map<const unsigned int, std::map<const unsigned int,
          std::map<const unsigned int, MMFFTor> > > >::const_iterator res2;
        std::map<const unsigned int,
          std::map<const unsigned int, std::map<const unsigned int,
          MMFFTor> > >::const_iterator res3;
        std::map<const unsigned int,
          std::map<const unsigned int, MMFFTor> >::const_iterator res4;
        std::map<const unsigned int, MMFFTor>::const_iterator res5;
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
            std::vector<boost::uint8_t>::const_iterator> jBounds;
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds;
        #endif
        
        while (((iter < maxIter) && ((!mmffTorParams) || (maxIter == 4)))
          || ((iter == 4) && (torType.first == 5) && torType.second)) {
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
          }
          else if (iter == 2) {
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
          }
          else if ((canJAtomType == canKAtomType)
            && (canIAtomType > canLAtomType)) {
            unsigned int temp = canLAtomType;
            canLAtomType = canIAtomType;
            canIAtomType = temp;
          }
          #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
          res1 = d_params.find(canTorType);
          if (res1 != d_params.end()) {
            res2 = ((*res1).second).find(canIAtomType);
            if (res2 != ((*res1).second).end()) {
              res3 = ((*res2).second).find(canJAtomType);
              if (res3 != ((*res2).second).end()) {
                res4 = ((*res3).second).find(canKAtomType);
                if (res4 != ((*res3).second).end()) {
                  res5 = ((*res4).second).find(canLAtomType);
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
          jBounds = std::equal_range(d_jAtomType.begin(), d_jAtomType.end(), canJAtomType);
          if (jBounds.first != jBounds.second) {
            bounds = std::equal_range
              (d_kAtomType.begin() + (jBounds.first - d_jAtomType.begin()),
              d_kAtomType.begin() + (jBounds.second - d_jAtomType.begin()),
              canKAtomType);
            if (bounds.first != bounds.second) {
              bounds = std::equal_range
                (d_iAtomType.begin() + (bounds.first - d_kAtomType.begin()),
                d_iAtomType.begin() + (bounds.second - d_kAtomType.begin()),
                canIAtomType);
              if (bounds.first != bounds.second) {
                bounds = std::equal_range
                  (d_lAtomType.begin() + (bounds.first - d_iAtomType.begin()),
                  d_lAtomType.begin() + (bounds.second - d_iAtomType.begin()),
                  canLAtomType);
                if (bounds.first != bounds.second) {
                  bounds = std::equal_range
                    (d_torType.begin() + (bounds.first - d_lAtomType.begin()),
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
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFTorCollection(const bool isMMFFs, std::string mmffTor);
      static class MMFFTorCollection *ds_instance[2];    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int,
        std::map<const unsigned int, std::map<const unsigned int,
        std::map<const unsigned int, std::map<const unsigned int,
        MMFFTor> > > > > d_params;  //!< the parameter 5D-map
      #else
      std::vector<MMFFTor> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_iAtomType;  //! atom type vector for atom i
      std::vector<boost::uint8_t> d_jAtomType;  //! atom type vector for atom j
      std::vector<boost::uint8_t> d_kAtomType;  //! atom type vector for atom k
      std::vector<boost::uint8_t> d_lAtomType;  //! atom type vector for atom l
      std::vector<boost::uint8_t> d_torType;  //! torsion type vector for angle i-j-k-l
      #endif
    };

    class MMFFVdWCollection {
    public:
      //! gets a pointer to the singleton MMFFVdWCollection
      /*!
	\param mmffVdW (optional) a string with parameter data. See
	 below for more information about this argument

	\return a pointer to the singleton MMFFVdWCollection

	<b>Notes:</b>
	  - do <b>not</b> delete the pointer returned here
	  - if the singleton MMFFVdWCollection has already been instantiated and
	    \c mmffVdW is empty, the singleton will be returned.
	  - if \c mmffVdW is empty and the singleton MMFFVdWCollection has
	    not yet been instantiated, the default parameters (from Params.cpp)
	    will be used.
	  - if \c mmffVdW is supplied, a new singleton will be instantiated.
	    The current instantiation (if there is one) will be deleted.
      */
      double power;
      double B;
      double Beta;
      double DARAD;
      double DAEPS;
      static MMFFVdWCollection *getMMFFVdW(const std::string &mmffVdW = "");
      //! Looks up the parameters for a particular key and returns them.
      /*!
	\return a pointer to the MMFFVdW object, NULL on failure.
      */
      const MMFFVdW *operator()(const unsigned int atomType) const {
        #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
        std::map<const unsigned int, MMFFVdW>::const_iterator res;
        res = d_params.find(atomType);

        return (res != d_params.end() ? &((*res).second) : NULL);
        #else
        std::pair<std::vector<boost::uint8_t>::const_iterator,
          std::vector<boost::uint8_t>::const_iterator> bounds =
          std::equal_range(d_atomType.begin(), d_atomType.end(), atomType);

        return ((bounds.first != bounds.second)
          ? &d_params[bounds.first - d_atomType.begin()] : NULL);
        #endif
      }
    private:
      //! to force this to be a singleton, the constructor must be private
      MMFFVdWCollection(std::string mmffVdW);
      static class MMFFVdWCollection *ds_instance;    //!< the singleton
      #ifdef RDKIT_MMFF_PARAMS_USE_STD_MAP
      std::map<const unsigned int, MMFFVdW> d_params;  //!< the parameter map
      #else
      std::vector<MMFFVdW> d_params;  //!< the parameter vector
      std::vector<boost::uint8_t> d_atomType;  //! atom type vector
      #endif
    };
  }
}

#endif
