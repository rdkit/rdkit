//
//  Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file MacroBondInfo.h

  \brief Defines bond-level macro bond information

*/
#include <RDGeneral/export.h>
#ifndef RD_MACROBONDINFO_H
#define RD_MACROBONDINFO_H

#include <cstddef>
#include <memory>
#include <vector>

namespace RDKit {

//! Captures bond-level information for macro bonds.
/*!
  Macro bond information is owned by a Bond. Since macro molecules may have
  multiple macro bonds between the same pair of atoms, this stores one entry
  for each macro-bond connection represented by the owning graph bond.
  Use this to retrieve macro bond types; the owning graph bond's bond type
  is not the macro bond type.
*/
class RDKIT_GRAPHMOL_EXPORT MacroBondInfo {
 public:
  //! One macro-bond connection in the owning graph bond's orientation.
  struct RDKIT_GRAPHMOL_EXPORT BondInfo {
    //! The attachment point on the begin atom.
    int beginAttachPt{-1};

    //! The attachment point on the end atom.
    int endAttachPt{-1};

    //! Stored as the numeric Bond::BondType value to avoid a Bond.h dependency.
    unsigned int bondType{0};
  };

  virtual ~MacroBondInfo() = default;

  MacroBondInfo() = default;

  //! Construct macro bond information with one macro-bond connection.
  /*!
    \param beginAttachPt the attachment point on the begin atom
    \param endAttachPt   the attachment point on the end atom
    \param bondType      the numeric Bond::BondType value
  */
  MacroBondInfo(int beginAttachPt, int endAttachPt, unsigned int bondType) {
    addBond(beginAttachPt, endAttachPt, bondType);
  }
  MacroBondInfo(const MacroBondInfo &other) = default;

  //! Adds a macro-bond connection.
  /*!
    \param beginAttachPt the attachment point on the begin atom
    \param endAttachPt   the attachment point on the end atom
    \param bondType      the numeric Bond::BondType value
  */
  void addBond(int beginAttachPt, int endAttachPt, unsigned int bondType) {
    d_bonds.push_back({beginAttachPt, endAttachPt, bondType});
  }

  //! Returns the number of macro-bond connections.
  std::size_t getNumBonds() const { return d_bonds.size(); }

  //! Returns a macro-bond connection by index.
  /*!
    \param idx the index of the macro-bond connection

    \return the requested macro-bond connection
  */
  const BondInfo &getBond(unsigned int idx) const { return d_bonds.at(idx); }

  //! Returns all macro-bond connections.
  const std::vector<BondInfo> &getBonds() const { return d_bonds; }

  //! Returns a copy of this macro bond information.
  virtual std::unique_ptr<MacroBondInfo> copy() const {
    return std::make_unique<MacroBondInfo>(*this);
  }

 private:
  std::vector<BondInfo> d_bonds;
};
}  // namespace RDKit

#endif
