//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Metal.h

        \brief Defines the MetalDisconnector class.

*/
#include <RDGeneral/export.h>
#ifndef __RD_METAL_H__
#define __RD_METAL_H__

#include <GraphMol/ROMol.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {
//! The MetalDisconnector class contains tools for disconnecting metal atoms
//! that are defined as covalently bonded to non-metals.
/*!

  <b>Notes:</b>
    -
*/

class RDKIT_MOLSTANDARDIZE_EXPORT MetalDisconnector {
 public:
  MetalDisconnector();
  MetalDisconnector(const MetalDisconnector &other);
  ~MetalDisconnector();

  ROMol *getMetalNof();  // {return metal_nof;}
  ROMol *getMetalNon();  // {return metal_non;}
  void setMetalNof(const ROMol &mol);
  void setMetalNon(const ROMol &mol);

  //! Break covalent bonds between metals and organic atoms under certain
  //! conditions.
  /*!
    <b>Notes:</b>
          The algorithm works as follows:
- Disconnect N, O, F from any metal.
- Disconnect other non-metals from transition metals + Al (but not Hg, Ga, Ge,
In, Sn, As, Tl, Pb, Bi, Po).
          - For every bond broken, adjust the charges of the begin and end atoms
accordingly.
  */
  ROMol *disconnect(const ROMol &mol);
  //! overload
  // modifies the molecule in place
  void disconnect(RWMol &mol);

 private:
  struct NonMetal {
    int cutBonds{0};
    std::vector<int> boundMetalIndices;
  };
  int chargeAdjustment(const Atom *a, int order);
  ROMOL_SPTR metal_nof;
  ROMOL_SPTR metal_non;

};  // class Metal
}  // namespace MolStandardize
}  // namespace RDKit
#endif
