//
//  Copyright (C) 2026 Tad Hurst, Schrödinger and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RD_MACROMOLTEMPLATE_H
#define RD_MACROMOLTEMPLATE_H

#include <RDGeneral/export.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SubstanceGroup.h>
#include <GraphMol/MacroAtomInfo.h>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace RDKit {
//! A template molecule for a macromolecule monomer.
/*!
  Wraps an RWMol and annotates it with SUP SGroups: a single main group that
  defines the monomer's core atoms and monomer class, plus zero or more
  leaving-group SGroups describing the atoms removed when the monomer is
  connected to its neighbors.
*/
class RDKIT_GRAPHMOL_EXPORT MacroMolTemplate {
 public:
  MacroMolTemplate() = default;
  explicit MacroMolTemplate(const RWMol &other) : d_mol(other) {}

  //! Returns a copy of the underlying template molecule.
  RWMol getMol() const;

  //! Sets the main SUP SGroup for this template.
  /*!
    \param atomIdxs      atom indices in the main monomer group
    \param monomerClass  monomer class used as the main SGroup CLASS value
  */
  unsigned int setMainGroup(const std::vector<unsigned int> &atomIdxs,
                            MonomerClass monomerClass);

  //! Adds a leaving-group SUP SGroup and records its attachment point.
  /*!
    \param atomIdxs        atom indices in the leaving group
    \param attachAtomIdx   atom in the main group used for the attachment
    \param leavingAtomIdx  atom in the leaving group used for the attachment
    \param attachPt        attachment point label, conventionally 1, 2, etc.
  */
  unsigned int addLeavingGroup(const std::vector<unsigned int> &atomIdxs,
                               unsigned int attachAtomIdx,
                               unsigned int leavingAtomIdx, const int attachPt);

  //! Returns the main SUP SGroup, or nullptr if one has not been set.
  const SubstanceGroup *getMainSgroup() const;

  //! Returns all leaving-group SUP SGroups.
  std::vector<const SubstanceGroup *> getLeavingGroups() const;

 private:
  SubstanceGroup *getMainSgroup();

  RWMol d_mol;
};

struct RDKIT_GRAPHMOL_EXPORT MacroMolEntry {
  MonomerClass monomerClass = MonomerClass::OTHER;  // e.g., AA, NA, CHEM
  std::string templateName;  // Name of the template, e.g., "ALA"
  std::string symbol;        // e.g., "A" for Alanine
  std::string original_data;  // Original definition (SMILES, SDF, etc.)
  // Parsed, annotated template molecule
  std::shared_ptr<MacroMolTemplate> molTemplate;
};

class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLibrary {
 public:
  void addEntry(const std::shared_ptr<MacroMolEntry> &macroMolEntry);
  //! Returns the entries ordered by descending main-group size.
  std::vector<std::shared_ptr<MacroMolEntry>> entries() const;
  //! Returns the entry for the given monomer class and template name, or
  //! nullptr if no such entry has been added.
  std::shared_ptr<MacroMolEntry> getByTemplateName(
      MonomerClass monomerClass, const std::string &templateName) const;
  //! Returns the entry for the given monomer class and symbol, or nullptr if
  //! no such entry has been added.
  std::shared_ptr<MacroMolEntry> getBySymbol(
      MonomerClass monomerClass, const std::string &symbol) const;

 private:
  using MacroMolTemplateKey = std::pair<MonomerClass, std::string>;
  using MacroMolTemplateEntry = std::pair<size_t, std::shared_ptr<MacroMolEntry>>;

  //! Orders entries by descending main-group size. Entries with equal size
  //! keep their insertion order, which std::multiset guarantees for
  //! equivalent elements.
  struct MainGroupSizeGreater {
    bool operator()(const MacroMolTemplateEntry &lhs,
                    const MacroMolTemplateEntry &rhs) const {
      return lhs.first > rhs.first;
    }
  };

  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolEntry>> byTemplateName;
  std::map<MacroMolTemplateKey, std::shared_ptr<MacroMolEntry>> bySymbol;
  std::multiset<MacroMolTemplateEntry, MainGroupSizeGreater> orderedEntries;
};

}  // namespace RDKit

#endif
