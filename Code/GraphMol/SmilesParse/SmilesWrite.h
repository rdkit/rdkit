//
//  Copyright (C) 2002-2020 Greg Landrum and Rational Discovery LLC
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

namespace RDKit {
class Atom;
class Bond;
class ROMol;
namespace SmilesWrite {

//! \brief returns the cxsmiles data for a molecule
RDKIT_SMILESPARSE_EXPORT std::string getCXExtensions(const ROMol &mol);

//! \brief returns true if the atom number is in the SMILES organic subset
RDKIT_SMILESPARSE_EXPORT bool inOrganicSubset(int atomicNumber);

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
RDKIT_SMILESPARSE_EXPORT std::string GetAtomSmiles(const Atom *atom,
                                                   bool doKekule = false,
                                                   const Bond *bondIn = nullptr,
                                                   bool allHsExplicit = false,
                                                   bool isomericSmiles = true);

//! \brief returns the SMILES for a bond
/*!
  \param bond : the bond to work with
  \param atomToLeftIdx : the index of the atom preceding \c bond
    in the SMILES
  \param doKekule : we're doing kekulized smiles (e.g. write out
    bond orders for aromatic bonds)
  \param allBondsExplicit : if true, symbols will be included for all bonds.
*/
RDKIT_SMILESPARSE_EXPORT std::string GetBondSmiles(
    const Bond *bond, int atomToLeftIdx = -1, bool doKekule = false,
    bool allBondsExplicit = false);
}  // namespace SmilesWrite

//! \brief returns canonical SMILES for a molecule
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
  atom.
 */
RDKIT_SMILESPARSE_EXPORT std::string MolToSmiles(
    const ROMol &mol, bool doIsomericSmiles = true, bool doKekule = false,
    int rootedAtAtom = -1, bool canonical = true, bool allBondsExplicit = false,
    bool allHsExplicit = false, bool doRandom = false);

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
RDKIT_SMILESPARSE_EXPORT std::string MolFragmentToSmiles(
    const ROMol &mol, const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr,
    const std::vector<std::string> *atomSymbols = nullptr,
    const std::vector<std::string> *bondSymbols = nullptr,
    bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1,
    bool canonical = true, bool allBondsExplicit = false,
    bool allHsExplicit = false);

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
  atom.
 */
RDKIT_SMILESPARSE_EXPORT std::string MolToCXSmiles(
    const ROMol &mol, bool doIsomericSmiles = true, bool doKekule = false,
    int rootedAtAtom = -1, bool canonical = true, bool allBondsExplicit = false,
    bool allHsExplicit = false, bool doRandom = false);

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
  \param doRandom : generate a randomized smiles string by randomly choosing
                    the priority to follow in the DFS traversal. [default false]

  \b NOTE: the bondSymbols are *not* currently used in the canonicalization.

 */
RDKIT_SMILESPARSE_EXPORT std::string MolFragmentToCXSmiles(
    const ROMol &mol, const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr,
    const std::vector<std::string> *atomSymbols = nullptr,
    const std::vector<std::string> *bondSymbols = nullptr,
    bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1,
    bool canonical = true, bool allBondsExplicit = false,
    bool allHsExplicit = false);

}  // namespace RDKit
#endif
