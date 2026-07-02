//
//  Copyright (C) 2026 Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file MacroAtomInfo.h

  \brief Defines atom-level macro atom information

*/
#include <RDGeneral/export.h>
#ifndef RD_MACROATOMINFO_H
#define RD_MACROATOMINFO_H

#include <memory>
#include <string>
#include <utility>

namespace RDKit {

//! Captures atom-level information for macro atoms.
/*!
  Macro atom information is owned by an Atom and stores the macro atom's
  display symbol and monomer class.
*/
class RDKIT_GRAPHMOL_EXPORT MacroAtomInfo {
 public:
  virtual ~MacroAtomInfo() = default;

  MacroAtomInfo() = default;

  //! Construct macro atom information.
  /*!
    \param symbol       the symbol used to identify the monomer
    \param monomerClass the class of monomer the macro atom represents
  */
  MacroAtomInfo(std::string symbol, std::string monomerClass = "")
      : d_symbol(std::move(symbol)), d_monomerClass(std::move(monomerClass)) {}
  MacroAtomInfo(const MacroAtomInfo &other) = default;

  //! Returns the macro atom symbol.
  const std::string &getSymbol() const { return d_symbol; }

  //! Sets the macro atom symbol.
  /*!
    \param symbol the symbol used to identify the monomer
  */
  void setSymbol(const std::string &symbol) { d_symbol = symbol; }

  //! Returns the macro atom monomer class.
  const std::string &getMonomerClass() const { return d_monomerClass; }

  //! Sets the macro atom monomer class.
  /*!
    \param monomerClass the class of monomer the macro atom represents
  */
  void setMonomerClass(const std::string &monomerClass) {
    d_monomerClass = monomerClass;
  }

  //! Returns a copy of this macro atom information.
  virtual std::unique_ptr<MacroAtomInfo> copy() const {
    return std::make_unique<MacroAtomInfo>(*this);
  }

 private:
  std::string d_symbol{""};
  std::string d_monomerClass{""};
};
}  // namespace RDKit

#endif
