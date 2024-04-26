//
//  Copyright (C) 2002-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SMILESWRITE_H_012020
#define RD_SMILESWRITE_H_012020

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <limits>

namespace RDKit {
class Atom;
class Bond;
class ROMol;

struct RDKIT_SMILESPARSE_EXPORT SmilesWriteParams {
  bool doIsomericSmiles =
      true;              /**< include stereochemistry and isotope information */
  bool doKekule = false; /**< kekulize the molecule before generating the SMILES
                            and output single/double bonds. NOTE that the output
                            is not canonical and that this will thrown an
                            exception if the molecule cannot be kekulized. */
  bool canonical = true; /**< generate canonical SMILES */
  bool allBondsExplicit = false; /**< include symbols for all bonds */
  bool allHsExplicit = false;    /**< provide hydrogen counts for every atom */
  bool doRandom = false; /**< randomize the output order. The resulting SMILES
                            is not canonical */
  int rootedAtAtom = -1; /**< make sure the SMILES starts at the specified
                             atom. The resulting SMILES is not canonical */
  bool includeDativeBonds =
      true; /**< include the RDKit extension for dative bonds. Otherwise dative
               bonds will be written as single bonds*/
};

namespace SmilesWrite {

#define CXSMILESFIELDS_ENUM_ITEMS                        \
  CXSMILESFIELDS_ENUM_ITEM(CX_NONE, 0)                   \
  CXSMILESFIELDS_ENUM_ITEM(CX_ATOM_LABELS, 1 << 0)       \
  CXSMILESFIELDS_ENUM_ITEM(CX_MOLFILE_VALUES, 1 << 1)    \
  CXSMILESFIELDS_ENUM_ITEM(CX_COORDS, 1 << 2)            \
  CXSMILESFIELDS_ENUM_ITEM(CX_RADICALS, 1 << 3)          \
  CXSMILESFIELDS_ENUM_ITEM(CX_ATOM_PROPS, 1 << 4)        \
  CXSMILESFIELDS_ENUM_ITEM(CX_LINKNODES, 1 << 5)         \
  CXSMILESFIELDS_ENUM_ITEM(CX_ENHANCEDSTEREO, 1 << 6)    \
  CXSMILESFIELDS_ENUM_ITEM(CX_SGROUPS, 1 << 7)           \
  CXSMILESFIELDS_ENUM_ITEM(CX_POLYMER, 1 << 8)           \
  CXSMILESFIELDS_ENUM_ITEM(CX_BOND_CFG, 1 << 9)          \
  CXSMILESFIELDS_ENUM_ITEM(CX_BOND_ATROPISOMER, 1 << 10) \
  CXSMILESFIELDS_ENUM_ITEM(CX_COORDINATE_BONDS, 1 << 11) \
  CXSMILESFIELDS_ENUM_ITEM(CX_ALL, 0x7fffffff)           \
  CXSMILESFIELDS_ENUM_ITEM(CX_ALL_BUT_COORDS, CX_ALL ^ CX_COORDS)

#define CXSMILESFIELDS_ENUM_ITEM(k, v) k = (v),
enum CXSmilesFields : uint32_t { CXSMILESFIELDS_ENUM_ITEMS };
#undef CXSMILESFIELDS_ENUM_ITEM
#define CXSMILESFIELDS_STD_MAP_ITEM(k) {#k, SmilesWrite::CXSmilesFields::k},
#define CXSMILESFIELDS_ENUM_ITEM(k, v) CXSMILESFIELDS_STD_MAP_ITEM(k)
#define CXSMILESFIELDS_ITEMS_MAP                       \
  std::map<std::string, SmilesWrite::CXSmilesFields> { \
    CXSMILESFIELDS_ENUM_ITEMS                          \
  }

//! \brief returns the cxsmiles data for a molecule
RDKIT_SMILESPARSE_EXPORT std::string getCXExtensions(
    const ROMol &mol, std::uint32_t flags = CXSmilesFields::CX_ALL);

//! \brief returns true if the atom number is in the SMILES organic subset
RDKIT_SMILESPARSE_EXPORT bool inOrganicSubset(int atomicNumber);

//! \brief returns the SMILES for an atom
/*!
  \param atom : the atom to work with
  \param ps : the parameters controlling the SMILES generation
*/
RDKIT_SMILESPARSE_EXPORT std::string GetAtomSmiles(const Atom *atom,
                                                   const SmilesWriteParams &ps);

//! \brief returns the SMILES for an atom
/*!
  \param atom : the atom to work with
  \param doKekule : we're doing kekulized smiles (e.g. don't use
    lower case for the atom label)
  \param bondIn : the bond we came into the atom on (unused)
  \param allHsExplicit : if true, hydrogen counts will be provided for every
  atom.
  \param isomericSmiles : if true, isomeric SMILES will be generated
*/
inline std::string GetAtomSmiles(const Atom *atom, bool doKekule = false,
                                 const Bond * = nullptr,
                                 bool allHsExplicit = false,
                                 bool isomericSmiles = true) {
  // RDUNUSED_PARAM(bondIn);
  SmilesWriteParams ps;
  ps.doIsomericSmiles = isomericSmiles;
  ps.doKekule = doKekule;
  ps.allHsExplicit = allHsExplicit;
  return GetAtomSmiles(atom, ps);
};

//! \brief returns the SMILES for a bond
/*!
  \param bond : the bond to work with
  \param ps : the parameters controlling the SMILES generation
  \param atomToLeftIdx : the index of the atom preceding \c bond
    in the SMILES
*/
RDKIT_SMILESPARSE_EXPORT std::string GetBondSmiles(const Bond *bond,
                                                   const SmilesWriteParams &ps,
                                                   int atomToLeftIdx = -1);
//! \brief returns the SMILES for a bond
/*!
  \param bond : the bond to work with
  \param atomToLeftIdx : the index of the atom preceding \c bond
    in the SMILES
  \param doKekule : we're doing kekulized smiles (e.g. write out
    bond orders for aromatic bonds)
  \param allBondsExplicit : if true, symbols will be included for all bonds.
*/
inline std::string GetBondSmiles(const Bond *bond, int atomToLeftIdx = -1,
                                 bool doKekule = false,
                                 bool allBondsExplicit = false) {
  SmilesWriteParams ps;
  ps.doKekule = doKekule;
  ps.allBondsExplicit = allBondsExplicit;
  ps.doIsomericSmiles = false;
  return GetBondSmiles(bond, ps, atomToLeftIdx);
};

namespace detail {
RDKIT_SMILESPARSE_EXPORT std::string MolToSmiles(
    const ROMol &mol, const SmilesWriteParams &params, bool doingCXSmiles);
}

}  // namespace SmilesWrite

//! \brief returns canonical SMILES for a molecule
RDKIT_SMILESPARSE_EXPORT std::string MolToSmiles(
    const ROMol &mol, const SmilesWriteParams &params);

//! \brief returns canonical SMILES for a molecule
/*!
  \param mol : the molecule in question.
  \param doIsomericSmiles : include stereochemistry and isotope information
      in the SMILES

  \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds) NOTE that
      this will throw an exception if the molecule cannot be kekulized.

  \param rootedAtAtom : make sure the SMILES starts at the specified atom.
      The resulting SMILES is not, of course, canonical.
  \param canonical : if false, no attempt will be made to canonicalize the
  SMILES
  \param allBondsExplicit : if true, symbols will be included for all bonds.
  \param allHsExplicit : if true, hydrogen counts will be provided for every
  atom.
 */
inline std::string MolToSmiles(const ROMol &mol, bool doIsomericSmiles = true,
                               bool doKekule = false, int rootedAtAtom = -1,
                               bool canonical = true,
                               bool allBondsExplicit = false,
                               bool allHsExplicit = false,
                               bool doRandom = false) {
  SmilesWriteParams ps;
  ps.doIsomericSmiles = doIsomericSmiles;
  ps.doKekule = doKekule;
  ps.rootedAtAtom = rootedAtAtom;
  ps.canonical = canonical;
  ps.allBondsExplicit = allBondsExplicit;
  ps.allHsExplicit = allHsExplicit;
  ps.doRandom = doRandom;
  return MolToSmiles(mol, ps);
};

//! \brief returns a vector of random SMILES for a molecule (may contain
//! duplicates)
/*!
  \param mol : the molecule in question.
  \param numSmiles : the number of SMILES to return
  \param randomSeed : if >0, will be used to seed the random number generator
  \param doIsomericSmiles : include stereochemistry and isotope information
      in the SMILES
  \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds)
  \param allBondsExplicit : if true, symbols will be included for all bonds.
  \param allHsExplicit : if true, hydrogen counts will be provided for every
  atom.
 */
RDKIT_SMILESPARSE_EXPORT std::vector<std::string> MolToRandomSmilesVect(
    const ROMol &mol, unsigned int numSmiles, unsigned int randomSeed = 0,
    bool doIsomericSmiles = true, bool doKekule = false,
    bool allBondsExplicit = false, bool allHsExplicit = false);

//! \brief returns canonical SMILES for part of a molecule
RDKIT_SMILESPARSE_EXPORT std::string MolFragmentToSmiles(
    const ROMol &mol, const SmilesWriteParams &params,
    const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr,
    const std::vector<std::string> *atomSymbols = nullptr,
    const std::vector<std::string> *bondSymbols = nullptr);

//! \brief returns canonical SMILES for part of a molecule
/*!
  \param mol : the molecule in question.
  \param atomsToUse : indices of the atoms in the fragment
  \param bondsToUse : indices of the bonds in the fragment. If this is not
  provided,
                      all bonds between the atoms in atomsToUse will be included
  \param atomSymbols : symbols to use for the atoms in the output SMILES
  \param bondSymbols : symbols to use for the bonds in the output SMILES
  \param doIsomericSmiles : include stereochemistry and isotope information
      in the SMILES
  \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds)
  \param rootedAtAtom : make sure the SMILES starts at the specified atom.
      The resulting SMILES is not, of course, canonical.
  \param canonical : if false, no attempt will be made to canonicalize the
  SMILES
  \param allBondsExplicit : if true, symbols will be included for all bonds.
  \param allHsExplicit : if true, hydrogen counts will be provided for every
  atom.
  \param doRandom : generate a randomized smiles string by randomly choosing
                    the priority to follow in the DFS traversal. [default false]

  \b NOTE: the bondSymbols are *not* currently used in the canonicalization.

 */
inline std::string MolFragmentToSmiles(
    const ROMol &mol, const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr,
    const std::vector<std::string> *atomSymbols = nullptr,
    const std::vector<std::string> *bondSymbols = nullptr,
    bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1,
    bool canonical = true, bool allBondsExplicit = false,
    bool allHsExplicit = false) {
  SmilesWriteParams ps;
  ps.doIsomericSmiles = doIsomericSmiles;
  ps.doKekule = doKekule;
  ps.rootedAtAtom = rootedAtAtom;
  ps.canonical = canonical;
  ps.allBondsExplicit = allBondsExplicit;
  ps.allHsExplicit = allHsExplicit;
  return MolFragmentToSmiles(mol, ps, atomsToUse, bondsToUse, atomSymbols,
                             bondSymbols);
}

#define RESTOREBONDDIROPTION_ENUM_ITEMS                          \
  RESTOREBONDDIROPTION_ENUM_ITEM(RestoreBondDirOptionTrue,       \
                                 0) /*!< DO restore bond dirs */ \
  RESTOREBONDDIROPTION_ENUM_ITEM(RestoreBondDirOptionClear,      \
                                 1) /*!< clear all bond dir information */

#define RESTOREBONDDIROPTION_ENUM_ITEM(k, v) k = v,
enum RestoreBondDirOption { RESTOREBONDDIROPTION_ENUM_ITEMS };
#undef RESTOREBONDDIROPTION_ENUM_ITEM
#define RESTOREBONDDIROPTION_STD_MAP_ITEM(k) {#k, k},
#define RESTOREBONDDIROPTION_ENUM_ITEM(k, v) \
  RESTOREBONDDIROPTION_STD_MAP_ITEM(k)
#define RESTOREBONDDIROPTION_ITEMS_MAP          \
  std::map<std::string, RestoreBondDirOption> { \
    RESTOREBONDDIROPTION_ENUM_ITEMS             \
  }

//! \brief returns canonical CXSMILES for a molecule
RDKIT_SMILESPARSE_EXPORT std::string MolToCXSmiles(
    const ROMol &mol, const SmilesWriteParams &ps,
    std::uint32_t flags = SmilesWrite::CXSmilesFields::CX_ALL,
    RestoreBondDirOption restoreBondDirs = RestoreBondDirOptionClear);

//! \brief returns canonical CXSMILES for a molecule
/*!
  \param mol : the molecule in question.
  \param doIsomericSmiles : include stereochemistry and isotope information
      in the SMILES
  \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds)
  \param rootedAtAtom : make sure the SMILES starts at the specified atom.
      The resulting SMILES is not, of course, canonical.
  \param canonical : if false, no attempt will be made to canonicalize the
  SMILES
  \param allBondsExplicit : if true, symbols will be included for all bonds.
  \param allHsExplicit : if true, hydrogen counts will be provided for every
  \param doRandom : generate a randomized smiles string by randomly choosing
                    the priority to follow in the DFS traversal. [default false]
  atom.
 */
inline std::string MolToCXSmiles(const ROMol &mol, bool doIsomericSmiles = true,
                                 bool doKekule = false, int rootedAtAtom = -1,
                                 bool canonical = true,
                                 bool allBondsExplicit = false,
                                 bool allHsExplicit = false,
                                 bool doRandom = false) {
  SmilesWriteParams ps;
  ps.doIsomericSmiles = doIsomericSmiles;
  ps.doKekule = doKekule;
  ps.rootedAtAtom = rootedAtAtom;
  ps.canonical = canonical;
  ps.allBondsExplicit = allBondsExplicit;
  ps.allHsExplicit = allHsExplicit;
  ps.doRandom = doRandom;
  return MolToCXSmiles(mol, ps, SmilesWrite::CXSmilesFields::CX_ALL);
};

//! \brief returns canonical CXSMILES for part of a molecule
RDKIT_SMILESPARSE_EXPORT std::string MolFragmentToCXSmiles(
    const ROMol &mol, const SmilesWriteParams &params,
    const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr,
    const std::vector<std::string> *atomSymbols = nullptr,
    const std::vector<std::string> *bondSymbols = nullptr);

//! \brief returns canonical CXSMILES for part of a molecule
/*!
  \param mol : the molecule in question.
  \param atomsToUse : indices of the atoms in the fragment
  \param bondsToUse : indices of the bonds in the fragment. If this is not
  provided,
                      all bonds between the atoms in atomsToUse will be included
  \param atomSymbols : symbols to use for the atoms in the output SMILES
  \param bondSymbols : symbols to use for the bonds in the output SMILES
  \param doIsomericSmiles : include stereochemistry and isotope information
      in the SMILES
  \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds)
  \param rootedAtAtom : make sure the SMILES starts at the specified atom.
      The resulting SMILES is not, of course, canonical.
  \param canonical : if false, no attempt will be made to canonicalize the
  SMILES
  \param allBondsExplicit : if true, symbols will be included for all bonds.
  \param allHsExplicit : if true, hydrogen counts will be provided for every
  atom.

  \b NOTE: the bondSymbols are *not* currently used in the canonicalization.

 */
inline std::string MolFragmentToCXSmiles(
    const ROMol &mol, const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr,
    const std::vector<std::string> *atomSymbols = nullptr,
    const std::vector<std::string> *bondSymbols = nullptr,
    bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1,
    bool canonical = true, bool allBondsExplicit = false,
    bool allHsExplicit = false) {
  SmilesWriteParams ps;
  ps.doIsomericSmiles = doIsomericSmiles;
  ps.doKekule = doKekule;
  ps.rootedAtAtom = rootedAtAtom;
  ps.canonical = canonical;
  ps.allBondsExplicit = allBondsExplicit;
  ps.allHsExplicit = allHsExplicit;
  return MolFragmentToCXSmiles(mol, ps, atomsToUse, bondsToUse, atomSymbols,
                               bondSymbols);
}

void updateSmilesWriteParamsFromJSON(SmilesWriteParams &params,
                                     const std::string &details_json);
void updateSmilesWriteParamsFromJSON(SmilesWriteParams &params,
                                     const char *details_json);
void updateCXSmilesFieldsFromJSON(SmilesWrite::CXSmilesFields &cxSmilesFields,
                                  RestoreBondDirOption &restoreBondDirs,
                                  const std::string &details_json);
void updateCXSmilesFieldsFromJSON(SmilesWrite::CXSmilesFields &cxSmilesFields,
                                  RestoreBondDirOption &restoreBondDirs,
                                  const char *details_json);

}  // namespace RDKit
#endif
