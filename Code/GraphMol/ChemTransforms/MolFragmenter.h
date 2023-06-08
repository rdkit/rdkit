//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_MOLFRAGMENTER_H__
#define _RD_MOLFRAGMENTER_H__

#include <istream>
#include <GraphMol/ROMol.h>

namespace RDKit {
namespace MolFragmenter {
struct RDKIT_CHEMTRANSFORMS_EXPORT FragmenterBondType {
  unsigned int atom1Label, atom2Label;
  unsigned int atom1Type, atom2Type;
  Bond::BondType bondType;
  ROMOL_SPTR query;
};

//! \brief Fragments a molecule by breaking a set of bonds
//!
/*!

  \param mol            - the molecule to be modified
  \param bondIndices    - indices of the bonds to be broken

  optional:
  \param addDummies  - toggles addition of dummy atoms to indicate where
  bonds were broken
  \param dummyLabels - used to provide the labels to be used for the dummies.
  the first element in each pair is the label for the dummy
  that replaces the bond's beginAtom, the second is for the
  dummy that replaces the bond's endAtom. If not provided, the
  dummies are labeled with atom indices.
  \param bondTypes - used to provide the bond type to use between the
  fragments and the dummy atoms. If not provided, defaults to single.
  \param nCutsPerAtom - used to return the number of bonds that were
   cut at each atom. Should be nAtoms long.

  \return a new ROMol with the modifications
  The client is responsible for deleting this molecule.

*/
RDKIT_CHEMTRANSFORMS_EXPORT ROMol *fragmentOnBonds(
    const ROMol &mol, const std::vector<unsigned int> &bondIndices,
    bool addDummies = true,
    const std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels =
        nullptr,
    const std::vector<Bond::BondType> *bondTypes = nullptr,
    std::vector<unsigned int> *nCutsPerAtom = nullptr);
//! \overload
RDKIT_CHEMTRANSFORMS_EXPORT ROMol *fragmentOnBonds(
    const ROMol &mol, const std::vector<FragmenterBondType> &bondPatterns,
    const std::map<unsigned int, ROMOL_SPTR> *atomEnvirons = nullptr,
    std::vector<unsigned int> *nCutsPerAtom = nullptr);
RDKIT_CHEMTRANSFORMS_EXPORT void fragmentOnSomeBonds(
    const ROMol &mol, const std::vector<unsigned int> &bondIndices,
    std::vector<ROMOL_SPTR> &resMols, unsigned int maxToCut = 1,
    bool addDummies = true,
    const std::vector<std::pair<unsigned int, unsigned int>> *dummyLabels =
        nullptr,
    const std::vector<Bond::BondType> *bondTypes = nullptr,
    std::vector<std::vector<unsigned int>> *nCutsPerAtom = nullptr);

//! \brief Fragments a molecule by breaking all BRICS bonds
/*!
  \return a new ROMol with the modifications
  The client is responsible for deleting this molecule.

*/
RDKIT_CHEMTRANSFORMS_EXPORT ROMol *fragmentOnBRICSBonds(const ROMol &mol);

RDKIT_CHEMTRANSFORMS_EXPORT void constructFragmenterAtomTypes(
    std::istream *inStream, std::map<unsigned int, std::string> &defs,
    const std::string &comment = "//", bool validate = true,
    std::map<unsigned int, ROMOL_SPTR> *environs = nullptr);
RDKIT_CHEMTRANSFORMS_EXPORT void constructFragmenterAtomTypes(
    const std::string &str, std::map<unsigned int, std::string> &defs,
    const std::string &comment = "//", bool validate = true,
    std::map<unsigned int, ROMOL_SPTR> *environs = nullptr);
RDKIT_CHEMTRANSFORMS_EXPORT void constructBRICSAtomTypes(
    std::map<unsigned int, std::string> &defs,
    std::map<unsigned int, ROMOL_SPTR> *environs = nullptr);
RDKIT_CHEMTRANSFORMS_EXPORT void constructFragmenterBondTypes(
    std::istream *inStream,
    const std::map<unsigned int, std::string> &atomTypes,
    std::vector<FragmenterBondType> &defs, const std::string &comment = "//",
    bool validate = true, bool labelByConnector = true);
RDKIT_CHEMTRANSFORMS_EXPORT void constructFragmenterBondTypes(
    const std::string &str,
    const std::map<unsigned int, std::string> &atomTypes,
    std::vector<FragmenterBondType> &defs, const std::string &comment = "//",
    bool validate = true, bool labelByConnector = true);
RDKIT_CHEMTRANSFORMS_EXPORT void constructBRICSBondTypes(
    std::vector<FragmenterBondType> &defs);
}  // namespace MolFragmenter

// n.b. AtomProperty must resolve to an unsigned integer value on an atom
// property
enum class MolzipLabel {
  AtomMapNumber,
  Isotope,
  FragmentOnBonds,
  AtomType,
  AtomProperty
};

struct RDKIT_CHEMTRANSFORMS_EXPORT MolzipParams {
  MolzipLabel label = MolzipLabel::AtomMapNumber;
  std::vector<std::string> atomSymbols;
  std::string atomProperty;
  bool enforceValenceRules = true;
  bool generateCoordinates = false;
};

RDKIT_CHEMTRANSFORMS_EXPORT std::unique_ptr<ROMol> molzip(
    const ROMol &a, const ROMol &b,
    const MolzipParams &params = MolzipParams());

RDKIT_CHEMTRANSFORMS_EXPORT std::unique_ptr<ROMol> molzip(
    const ROMol &a, const MolzipParams &params = MolzipParams());

//! \brief Creates a molecule from an R group decomposition
/*!
 *
 * @param decomposition - A list of molecules that comprises an R group
 * decomposition.  The core must be the first molecule in the list. If
 * generateCoordinates is set in the parameters then aligned depiction
 * coordinates will be set on the returned molecule and the input decomposition
 *
 * optional:
 * @param params - molzip parameters
 *
 * @return the zipped molecule
 */
RDKIT_CHEMTRANSFORMS_EXPORT std::unique_ptr<ROMol> molzip(
    std::vector<ROMOL_SPTR> &decomposition,
    const MolzipParams &params = MolzipParams());
}  // namespace RDKit
#endif
