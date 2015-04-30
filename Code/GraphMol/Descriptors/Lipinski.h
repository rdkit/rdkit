//
//  Copyright (C) 2007-2011 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file Lipinski.h

  \brief Contains Lipinski and Lipinski-like descriptors. Use MolDescriptors.h in client code.

*/
#ifndef __RD_LIPINSKI_H__
#define __RD_LIPINSKI_H__
namespace RDKit{
  class ROMol;
  namespace Descriptors {

    const std::string lipinskiHBAVersion="1.0.0";
    //! calculates the standard Lipinski HBA definition (number of Ns and Os)
    unsigned int calcLipinskiHBA(const ROMol &mol);

    const std::string lipinskiHBDVersion="2.0.0";
    //! calculates the standard Lipinski HBA definition (number of N-H and O-H bonds)
    unsigned int calcLipinskiHBD(const ROMol &mol);

    extern const std::string NumRotatableBondsVersion;
    //! calculates the number of rotatable bonds
    /*!
      \param mol           the molecule of interest
      \param strict        if set, a stricter definition of rotable bonds is used
                           this excludes amides, esters, etc.
    */
    unsigned int calcNumRotatableBonds(const ROMol &mol,bool useStrictDefinition=true);

    extern const std::string NumHBDVersion;
    //! calculates the number of H-bond donors
    unsigned int calcNumHBD(const ROMol &mol);

    extern const std::string NumHBAVersion;
    //! calculates the number of H-bond acceptors
    unsigned int calcNumHBA(const ROMol &mol);

    extern const std::string NumHeteroatomsVersion;
    //! calculates the number of heteroatoms
    unsigned int calcNumHeteroatoms(const ROMol &mol);

    extern const std::string NumAmideBondsVersion;
    //! calculates the number of amide bonds
    unsigned int calcNumAmideBonds(const ROMol &mol);

    extern const std::string FractionCSP3Version;
    //! calculates the fraction of carbons that are SP3 hybridized
    double calcFractionCSP3(const ROMol &mol);

    extern const std::string NumRingsVersion;
    //! calculates the number of SSSR rings
    unsigned int calcNumRings(const ROMol &mol);

    extern const std::string NumAromaticRingsVersion;
    //! calculates the number of aromatic SSSR rings
    unsigned int calcNumAromaticRings(const ROMol &mol);

    extern const std::string NumAliphaticRingsVersion;
    //! calculates the number of aliphatic (at least one non-aromatic bond) SSSR rings
    unsigned int calcNumAliphaticRings(const ROMol &mol);

    extern const std::string NumSaturatedRingsVersion;
    //! calculates the number of saturated SSSR rings
    unsigned int calcNumSaturatedRings(const ROMol &mol);

    extern const std::string NumHeterocyclesVersion;
    //! calculates the number of SSSR heterocycles
    unsigned int calcNumHeterocycles(const ROMol &mol);

    extern const std::string NumAromaticHeterocyclesVersion;
    //! calculates the number of aromatic SSSR heterocycles
    unsigned int calcNumAromaticHeterocycles(const ROMol &mol);

    extern const std::string NumAromaticCarbocyclesVersion;
    //! calculates the number of aromatic SSSR carbocycles
    unsigned int calcNumAromaticCarbocycles(const ROMol &mol);

    extern const std::string NumSaturatedHeterocyclesVersion;
    //! calculates the number of saturated SSSR heterocycles
    unsigned int calcNumSaturatedHeterocycles(const ROMol &mol);

    extern const std::string NumSaturatedCarbocyclesVersion;
    //! calculates the number of saturated SSSR carbocycles
    unsigned int calcNumSaturatedCarbocycles(const ROMol &mol);

    extern const std::string NumAliphaticHeterocyclesVersion;
    //! calculates the number of aliphatic (at least one non-aromatic bond) SSSR heterocycles
    unsigned int calcNumAliphaticHeterocycles(const ROMol &mol);

    extern const std::string NumAliphaticCarbocyclesVersion;
    //! calculates the number of aliphatic (at least one non-aromatic bond) SSSR carbocycles
    unsigned int calcNumAliphaticCarbocycles(const ROMol &mol);
  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
