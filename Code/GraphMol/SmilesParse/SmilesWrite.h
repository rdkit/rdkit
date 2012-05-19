//
//  Copyright (C) 2002-2012 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SMILESWRITE_H
#define _RD_SMILESWRITE_H

#include <string>
#include <vector>

namespace RDKit{
  class Atom;
  class Bond;
  class ROMol;
  namespace SmilesWrite {
    //! \brief returns true if the atom number is in the SMILES organic subset
    bool inOrganicSubset(int atomicNumber);

    //! \brief returns the SMILES for an atom
    /*!
      \param atom : the atom to work with
      \param doKekule : we're doing kekulized smiles (e.g. don't use
        lower case for the atom label)
      \param bondIn : the bond we came into the atom on (used for
        chirality calculation
    */
    std::string GetAtomSmiles(const Atom *atom,bool doKekule=false,
                              const Bond *bondIn=0);

    //! \brief returns the SMILES for a bond
    /*!
      \param bond : the bond to work with
      \param atomToLeftIdx : the index of the atom preceding \c bond
        in the SMILES
      \param doKekule : we're doing kekulized smiles (e.g. write out
        bond orders for aromatic bonds)
      \param allBondsExplicit : if true, symbols will be included for all bonds.
    */
    std::string GetBondSmiles(const Bond *bond,int atomToLeftIdx=-1,
                              bool doKekule=false,bool allBondsExplicit=false);
  } 
  
  //! \brief returns canonical SMILES for a molecule
  /*!
    \param mol : the molecule in question. 
    \param doIsomericSmiles : include stereochemistry and isotope information
        in the SMILES
    \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds)
    \param rootedAtAtom : make sure the SMILES starts at the specified atom.
        The resulting SMILES is not, of course, canonical.
    \param canonical : if false, no attempt will be made to canonicalize the SMILES
    \param allBondsExplicit : if true, symbols will be included for all bonds.
   */
  std::string MolToSmiles(const ROMol &mol,bool doIsomericSmiles=false,
			  bool doKekule=false,int rootedAtAtom=-1,
                          bool canonical=true,
                          bool allBondsExplicit=false);

  //! \brief returns canonical SMILES for part of a molecule
  /*!
    \param mol : the molecule in question.
    \param atomsToUse : indices of the atoms in the fragment
    \param bondsToUse : indices of the bonds in the fragment. If this is not provided,
                        all bonds between the atoms in atomsToUse will be included
    \param atomSymbols : symbols to use for the atoms in the output SMILES
    \param bondSymbols : sybmols to use for the bonds in the output SMILES
    \param doIsomericSmiles : include stereochemistry and isotope information
        in the SMILES
    \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds)
    \param rootedAtAtom : make sure the SMILES starts at the specified atom.
        The resulting SMILES is not, of course, canonical.
    \param canonical : if false, no attempt will be made to canonicalize the SMILES
    \param allBondsExplicit : if true, symbols will be included for all bonds.
   */
  std::string MolFragmentToSmiles(const ROMol &mol,
                                  const std::vector<int> &atomsToUse,
                                  const std::vector<int> *bondsToUse=0,
                                  const std::vector<std::string> *atomSymbols=0,
                                  const std::vector<std::string> *bondSymbols=0,
                                  bool doIsomericSmiles=false,
                                  bool doKekule=false,int rootedAtAtom=-1,
                                  bool canonical=true,
                                  bool allBondsExplicit=false);

}
#endif
