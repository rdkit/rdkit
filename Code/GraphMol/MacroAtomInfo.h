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
#ifndef RD_MACROATOMINFO_H
#define RD_MACROATOMINFO_H

#include <RDGeneral/export.h>

#include <memory>
#include <string>
#include <utility>

namespace RDKit {

//! Format-neutral classes of monomer that a macro atom can represent.
/*!
  These are broad semantic categories, not the class vocabularies used by
  SCSR, HELM, PDB, or another external format. Format parsers and writers are
  responsible for mapping their own classifications to and from these values.
  Such mappings may require additional format-specific context and should fail
  rather than guess when the available information is ambiguous.

  A regular enum is used here because MonomerClass appears in public MacroMol
  method signatures. RDKit's BETTER_ENUM macro can expand to different C++
  types depending on per-source-file preprocessor settings, which can cause
  link errors in shared-library and unity builds.
*/
enum MonomerClass : int {
  //! An amino-acid monomer, without encoding chirality, modification, or
  //! crosslink status.
  AminoAcid,
  //! A nucleic-acid monomer, without distinguishing DNA, RNA, base, sugar, or
  //! phosphate classifications.
  NucleicAcid,
  //! A non-biopolymer chemical monomer.
  Chemical,
  //! A recognized monomer that does not belong to another neutral category.
  Other
};

//! Converts a monomer class to its canonical format-neutral RDKit name.
/*!
  The returned name is not a mapping to an external format's class vocabulary.
*/
RDKIT_GRAPHMOL_EXPORT const char *monomerClassToString(
    MonomerClass monomerClass);

//! Converts a canonical format-neutral RDKit name to its monomer class.
/*!
  \throws ValueErrorException if monomerClass is not a canonical name. External
          format adapters must validate and map their own class vocabulary
          instead of passing arbitrary values to this function.
*/
RDKIT_GRAPHMOL_EXPORT MonomerClass
monomerClassFromString(const std::string &monomerClass);

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
  MacroAtomInfo(std::string symbol,
                MonomerClass monomerClass = MonomerClass::Other)
      : d_symbol(std::move(symbol)), d_monomerClass(monomerClass) {}
  MacroAtomInfo(const MacroAtomInfo &other) = default;

  //! Returns the macro atom symbol.
  const std::string &getSymbol() const { return d_symbol; }

  //! Sets the macro atom symbol.
  /*!
    \param symbol the symbol used to identify the monomer
  */
  void setSymbol(const std::string &symbol) { d_symbol = symbol; }

  //! Returns the macro atom monomer class.
  MonomerClass getMonomerClass() const { return d_monomerClass; }

  //! Sets the macro atom monomer class.
  /*!
    \param monomerClass the class of monomer the macro atom represents
  */
  void setMonomerClass(MonomerClass monomerClass) {
    d_monomerClass = monomerClass;
  }

  //! Returns a copy of this macro atom information.
  virtual std::unique_ptr<MacroAtomInfo> copy() const {
    return std::make_unique<MacroAtomInfo>(*this);
  }

 private:
  std::string d_symbol{""};
  MonomerClass d_monomerClass{MonomerClass::Other};
};
}  // namespace RDKit

#endif
