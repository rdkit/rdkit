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
#ifndef _RD_MMFFATOMTYPER_H__
#define _RD_MMFFATOMTYPER_H__

#include <vector>
#include <string>
#include <iostream>
#include <ForceField/MMFF/Params.h>
#include <boost/cstdint.hpp>


namespace RDKit {
  class ROMol;
  class RWMol;
  class Atom;
  class Bond;

  namespace MMFF {
    class MMFFAtomProperties {
    public:
      MMFFAtomProperties() :
        mmffAtomType(0),
        mmffFormalCharge(0.0),
        mmffPartialCharge(0.0) {};
      ~MMFFAtomProperties() {};
      boost::uint8_t mmffAtomType;
      double mmffFormalCharge;
      double mmffPartialCharge;
    };
    
    typedef boost::shared_ptr<MMFFAtomProperties> MMFFAtomPropertiesPtr;
    enum {
      CONSTANT = 1,
      DISTANCE = 2
    };
    enum {
      MMFF_VERBOSITY_NONE = 0,
      MMFF_VERBOSITY_LOW = 1,
      MMFF_VERBOSITY_HIGH = 2
    };
    class MMFFMolProperties {
    public:
      MMFFMolProperties(unsigned int numAtoms) :
        d_mmffs(false),
        d_bondTerm(true),
        d_angleTerm(true),
        d_stretchBendTerm(true),
        d_oopTerm(true),
        d_torsionTerm(true),
        d_vdWTerm(true),
        d_eleTerm(true),
        d_dielConst(1.0),       //!< the dielectric constant
        d_dielModel(CONSTANT), //!< the dielectric model (1 = constant, 2 = distance-dependent)
        d_verbosity(MMFF_VERBOSITY_NONE),
        d_oStream(&(std::cout)),
        d_MMFFAtomPropertiesPtrVect(numAtoms) {
          for (unsigned int i = 0; i < numAtoms; ++i) {
            d_MMFFAtomPropertiesPtrVect[i]
              = MMFFAtomPropertiesPtr(new MMFFAtomProperties());
          }
        };
      ~MMFFMolProperties() {};
      const boost::uint8_t setMMFFHeavyAtomType(const Atom *atom);
      const boost::uint8_t setMMFFHydrogenType(const Atom *atom);
      const unsigned int getMMFFBondType(const Bond *bond);
      const unsigned int getMMFFAngleType(const ROMol &mol,
        const unsigned int idx1, const unsigned int idx2,
        const unsigned int idx3);
      const std::pair<unsigned int, unsigned int> getMMFFTorsionType
        (const ROMol &mol, const unsigned int idx1, const unsigned int idx2,
        const unsigned int idx3, const unsigned int idx4);
      void computeMMFFCharges(const ROMol *mol);
      const ForceFields::MMFF::MMFFTor *getMMFFTorsionEmpiricalRuleParams
        (const ROMol &mol, unsigned int idx2, unsigned int idx3);
      const ForceFields::MMFF::MMFFBond *getMMFFBondStretchEmpiricalRuleParams
        (const ROMol &mol, const Bond *bond);
      const boost::uint8_t getMMFFAtomType(const unsigned int idx)
      {
        RANGE_CHECK(0, idx, this->d_MMFFAtomPropertiesPtrVect.size() - 1);
        
        return this->d_MMFFAtomPropertiesPtrVect[idx]->mmffAtomType;
      };
      void setMMFFFormalCharge(const unsigned int idx, const double fChg)
      {
        RANGE_CHECK(0, idx, this->d_MMFFAtomPropertiesPtrVect.size() - 1);
        
        this->d_MMFFAtomPropertiesPtrVect[idx]->mmffFormalCharge = fChg;
      };
      const double getMMFFFormalCharge(const unsigned int idx)
      {
        RANGE_CHECK(0, idx, this->d_MMFFAtomPropertiesPtrVect.size() - 1);
        
        return this->d_MMFFAtomPropertiesPtrVect[idx]->mmffFormalCharge;
      };
      void setMMFFPartialCharge(const unsigned int idx, const double pChg)
      {
        RANGE_CHECK(0, idx, this->d_MMFFAtomPropertiesPtrVect.size() - 1);
        
        this->d_MMFFAtomPropertiesPtrVect[idx]->mmffPartialCharge = pChg;
      };
      const double getMMFFPartialCharge(const unsigned int idx)
      {
        RANGE_CHECK(0, idx, this->d_MMFFAtomPropertiesPtrVect.size() - 1);
        
        return this->d_MMFFAtomPropertiesPtrVect[idx]->mmffPartialCharge;
      };
      void setMMFFBondTerm(const bool state)
      {
        this->d_bondTerm = state;
      };
      const bool getMMFFBondTerm()
      {
        return this->d_bondTerm;
      };
      void setMMFFAngleTerm(const bool state)
      {
        this->d_angleTerm = state;
      };
      const bool getMMFFAngleTerm()
      {
        return this->d_angleTerm;
      };
      void setMMFFStretchBendTerm(const bool state)
      {
        this->d_stretchBendTerm = state;
      };
      const bool getMMFFStretchBendTerm()
      {
        return this->d_stretchBendTerm;
      };
      void setMMFFOopTerm(const bool state)
      {
        this->d_oopTerm = state;
      };
      const bool getMMFFOopTerm()
      {
        return this->d_oopTerm;
      };
      void setMMFFTorsionTerm(const bool state)
      {
        this->d_torsionTerm = state;
      };
      const bool getMMFFTorsionTerm()
      {
        return this->d_torsionTerm;
      };
      void setMMFFVdWTerm(const bool state)
      {
        this->d_vdWTerm = state;
      };
      const bool getMMFFVdWTerm()
      {
        return this->d_vdWTerm;
      };
      void setMMFFEleTerm(const bool state)
      {
        this->d_eleTerm = state;
      };
      const bool getMMFFEleTerm()
      {
        return this->d_eleTerm;
      };
      void setMMFFVariant(const std::string mmffVariant)
      {
        PRECONDITION((mmffVariant == "MMFF94")
          || (mmffVariant == "MMFF94s"), "bad MMFF variant");
        
        this->d_mmffs = ((mmffVariant == "MMFF94s") ? true : false);
      };
      const std::string getMMFFVariant()
      {
        return (this->d_mmffs ? "MMFF94s" : "MMFF94");
      };
      void setMMFFDielectricConstant(const double dielConst)
      {
        PRECONDITION(dielConst > 0.0, "bad dielectric constant");
        
        this->d_dielConst = dielConst;
      };
      double getMMFFDielectricConstant()
      {
        return this->d_dielConst;
      };
      void setMMFFDielectricModel(boost::uint8_t dielModel)
      {
        this->d_dielModel = dielModel;
      };
      boost::uint8_t getMMFFDielectricModel()
      {
        return this->d_dielModel;
      };
      void setMMFFVerbosity(boost::uint8_t verbosity)
      {
        this->d_verbosity = verbosity;
      };
      boost::uint8_t getMMFFVerbosity()
      {
        return this->d_verbosity;
      };
      void setMMFFOStream(std::ostream *oStream)
      {
        this->d_oStream = oStream;
      };
      std::ostream& getMMFFOStream()
      {
        return *(this->d_oStream);
      };
    private:
      bool d_mmffs;
      bool d_bondTerm;
      bool d_angleTerm;
      bool d_stretchBendTerm;
      bool d_oopTerm;
      bool d_torsionTerm;
      bool d_vdWTerm;
      bool d_eleTerm;
      double d_dielConst;       //!< the dielectric constant
      boost::uint8_t d_dielModel; //!< the dielectric model (1 = constant, 2 = distance-dependent)
      boost::uint8_t d_verbosity;
      std::ostream *d_oStream;
      std::vector<MMFFAtomPropertiesPtr> d_MMFFAtomPropertiesPtrVect;
    };
    MMFFMolProperties *setupMMFFForceField(ROMol *mol,
      std::string mmffVariant = "MMFF94", boost::uint8_t verbosity = MMFF_VERBOSITY_NONE,
      std::ostream &oStream = std::cout);
    unsigned int isAngleInRingOfSize3or4(const ROMol &mol, const unsigned int idx1,
      const unsigned int idx2, const unsigned int idx3);
    unsigned int isTorsionInRingOfSize4or5(const ROMol &mol, const unsigned int idx1,
      const unsigned int idx2, const unsigned int idx3, const unsigned int idx4);
    bool isAtomInAromaticRingOfSize(const Atom *atom, const unsigned int ringSize);
    bool isAtomNOxide(const Atom *atom);
    bool areAtomsInSameAromaticRing(const ROMol &mol,
      const unsigned int idx1, const unsigned int idx2);
    bool areAtomsInSameRingOfSize(const ROMol &mol,
      const unsigned int ringSize, const unsigned int numAtoms, ...);
    unsigned int sanitizeMMFFMol(RWMol &mol);
    void setMMFFAromaticity(RWMol *mol);
    const unsigned int getMMFFStretchBendType(const unsigned int angleType,
      const unsigned int bondType1, const unsigned int bondType2);
    const unsigned int getPeriodicTableRow(const int atomicNum);
    const ForceFields::MMFF::MMFFAngle *getMMFFAngleBendEmpiricalRuleParams
      (const ROMol &mol, const ForceFields::MMFF::MMFFAngle *oldMMFFAngleParams,
      const ForceFields::MMFF::MMFFProp *mmffPropParamsCentralAtom,
      const ForceFields::MMFF::MMFFBond *mmffBondParams1,
      const ForceFields::MMFF::MMFFBond *mmffBondParams2,
      unsigned int idx1, unsigned int idx2, unsigned int idx3);
  }
}


#endif
