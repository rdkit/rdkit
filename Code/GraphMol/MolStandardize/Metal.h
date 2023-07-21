//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
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
#ifndef RD_METAL_H
#define RD_METAL_H

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

struct RDKIT_MOLSTANDARDIZE_EXPORT MetalDisconnectorOptions {
  bool splitGrignards = false;  // Whether to split Grignard-type complexes.
  bool splitAromaticC = false;  // Whether to split metal-aromatic C bonds.
  bool adjustCharges = true;    // Whether to adjust charges on ligand atoms.
  bool removeHapticDummies =
      false;  // Whether to remove the dummy atoms representing haptic bonds.
              // Such dummies are bonded to the metal with a bond
              // that has the _MolFileBondEndPts prop set.
};

class RDKIT_MOLSTANDARDIZE_EXPORT MetalDisconnector {
 public:
  MetalDisconnector(
      const MetalDisconnectorOptions &options = MetalDisconnectorOptions());
  MetalDisconnector(const MetalDisconnector &other);
  ~MetalDisconnector();

  ROMol *getMetalNof();  // {return dp_metal_nof;}
  ROMol *getMetalNon();  // {return dp_metal_non;}
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
  /// modifies the molecule in place
  void disconnect(RWMol &mol);
  void disconnectInPlace(RWMol &mol) { disconnect(mol); };

 private:
  struct NonMetal {
    int cutBonds{0};
    std::vector<int> boundMetalIndices;
  };
  int chargeAdjustment(const Atom *a, int order);
  ROMOL_SPTR dp_metal_nof;
  ROMOL_SPTR dp_metal_non;
  ROMOL_SPTR dp_metalDummy;

  const MetalDisconnectorOptions d_options;

  void adjust_charges(RDKit::RWMol &mol, std::map<int, NonMetal> &nonMetals,
                      std::map<int, int> &metalChargeExcess);
  // Remove any dummy atoms that are bonded to a metal and have the ENDPTS
  // prop.  These are assumed to marking a haptic bond from the aotms in
  // ENDPTS to the metal, e.g. in ferrocene.
  void remove_haptic_dummies(RDKit::RWMol &mol);

};  // class Metal

}  // namespace MolStandardize
}  // namespace RDKit
#endif
