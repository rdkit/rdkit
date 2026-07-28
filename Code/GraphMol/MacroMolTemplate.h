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
class RDKIT_GRAPHMOL_EXPORT MacroMolTemplate : public RWMol {
 public:
  //! Constructs a template from a molecule and its immutable lookup metadata.
  MacroMolTemplate(const RWMol &mol, MonomerClass monomerClass,
                   std::string templateName, std::string symbol,
                   std::string originalData)
      : RWMol(mol),
        d_monomerClass(monomerClass),
        d_templateName(std::move(templateName)),
        d_symbol(std::move(symbol)),
        d_originalData(std::move(originalData)) {}

  //! Returns the monomer class used for lookup and the main SGroup CLASS.
  MonomerClass getMonomerClass() const { return d_monomerClass; }
  //! Returns the template name, e.g. "ALA".
  const std::string &getTemplateName() const { return d_templateName; }
  //! Returns the template symbol, e.g. "A" for alanine.
  const std::string &getSymbol() const { return d_symbol; }
  //! Returns the original template definition (SMILES, SDF, etc.).
  const std::string &getOriginalData() const { return d_originalData; }

  //! Sets the main SUP SGroup for this template.
  /*!
    The SGroup CLASS value is derived from this template's monomer class.

    \param atomIdxs  atom indices in the main monomer group
  */
  unsigned int setMainGroup(const std::vector<unsigned int> &atomIdxs);

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

  MonomerClass d_monomerClass;
  std::string d_templateName;
  std::string d_symbol;
  std::string d_originalData;
};

class RDKIT_GRAPHMOL_EXPORT MacroMolTemplateLibrary {
 public:
  //! Adds a completed template and takes ownership of it.
  /*!
    The template must contain exactly one main SUP SGroup. Registered
    templates are exposed as const objects so that their lookup keys and
    main-group-size ordering cannot be invalidated.
  */
  void addTemplate(std::unique_ptr<MacroMolTemplate> macroMolTemplate);

  //! Returns the templates ordered by descending main-group size.
  /*!
    Templates with equal main-group size remain in insertion order.
  */
  const std::vector<const MacroMolTemplate *> &entries() const;

  //! Returns the template for the given monomer class and template name, or
  //! nullptr if no such template has been added.
  const MacroMolTemplate *getByTemplateName(
      MonomerClass monomerClass, const std::string &templateName) const;

  //! Returns the template for the given monomer class and symbol, or nullptr
  //! if no such template has been added.
  const MacroMolTemplate *getBySymbol(
      MonomerClass monomerClass, const std::string &symbol) const;

 private:
  using MacroMolTemplateKey = std::pair<MonomerClass, std::string>;

  std::map<MacroMolTemplateKey, const MacroMolTemplate *> byTemplateName;
  std::map<MacroMolTemplateKey, const MacroMolTemplate *> bySymbol;
  std::vector<std::unique_ptr<const MacroMolTemplate>> ownedTemplates;
  std::vector<const MacroMolTemplate *> orderedEntries;
};

}  // namespace RDKit

#endif
