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
    //! calculates the standard Lipinski HBA definition
    /*!  
      \param mol        the molecule of interest

      \return the number of Ns and Os in the molecule

    */
    unsigned int calcLipinskiHBA(const ROMol &mol);

    const std::string lipinskiHBDVersion="2.0.0";
    //! calculates the standard Lipinski HBA definition
    /*!  
      \param mol        the molecule of interest

      \return the number of N-H and O-H bonds in the molecule

    */
    unsigned int calcLipinskiHBD(const ROMol &mol);

    extern const std::string NumRotatableBondsVersion;
    unsigned int calcNumRotatableBonds(const ROMol &mol);
    extern const std::string NumHBDVersion;
    unsigned int calcNumHBD(const ROMol &mol);
    extern const std::string NumHBAVersion;
    unsigned int calcNumHBA(const ROMol &mol);
    extern const std::string NumHeteroatomsVersion;
    unsigned int calcNumHeteroatoms(const ROMol &mol);

    extern const std::string NumRingsVersion;
    unsigned int calcNumRings(const ROMol &mol);


  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
