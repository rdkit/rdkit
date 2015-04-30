//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_UFFATOMTYPER_H__
#define _RD_UFFATOMTYPER_H__

#include <vector>
#include <string>
#include <ForceField/UFF/Params.h>

namespace ForceFields {
  namespace UFF {
    class AtomicParams;
  }
}
namespace RDKit {
  class ROMol;
  class Atom;

  namespace UFF {
    typedef std::vector<const ForceFields::UFF::AtomicParams *> AtomicParamVect;

    std::pair<AtomicParamVect,bool> getAtomTypes(const ROMol &mol,const std::string &paramData="");
    bool getUFFBondStretchParams(const ROMol &mol, unsigned int idx1,
      unsigned int idx2, ForceFields::UFF::UFFBond &uffBondStretchParams);
    bool getUFFAngleBendParams(const ROMol &mol, unsigned int idx1,
      unsigned int idx2, unsigned int idx3, ForceFields::UFF::UFFAngle &uffAngleBendParams);
    bool getUFFTorsionParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
      unsigned int idx3, unsigned int idx4, ForceFields::UFF::UFFTor &uffTorsionParams);
    bool getUFFInversionParams(const ROMol &mol, unsigned int idx1, unsigned int idx2,
      unsigned int idx3, unsigned int idx4, ForceFields::UFF::UFFInv &uffInversionParams);
    bool getUFFVdWParams(const ROMol &mol, unsigned int idx1,
      unsigned int idx2, ForceFields::UFF::UFFVdW &uffVdWParams);

    namespace Tools {
      // these functions are primarily exposed so they can be tested.
      void addAtomChargeFlags(const Atom *atom,std::string &atomKey,
			      bool tolerateChargeMismatch=true);
      std::string getAtomLabel(const Atom *atom);
    }

  }
}


#endif
