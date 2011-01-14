//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
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
    */
    std::string GetBondSmiles(const Bond *bond,int atomToLeftIdx=-1,
                              bool doKekule=false);
  } 
  
  //! \brief returns canonical SMILES for a molecule
  /*!
    \param mol : the molecule in question. NOTE that the molecule may
        be modified as part of the canonicalization process.
    \param doIsomericSmiles : include stereochemistry and isotope information
        in the SMILES
    \param doKekule : do Kekule smiles (i.e. don't use aromatic bonds)
    \param rootedAtAtom : make sure the SMILES starts at the specified atom.
        The resulting SMILES is not, of course, canonical.
    \param canonical : if false, no attempt will be made to canonicalize the SMILES
   */
  std::string MolToSmiles(ROMol &mol,bool doIsomericSmiles=false,
			  bool doKekule=false,int rootedAtAtom=-1,
                          bool canonical=true);

  //! \brief returns SMILES for a piece of a molecule
  /*!
    NOTE: the SMILES generated is, at the moment, not canonical
    
    \param mol : the molecule in question. 
    \param rootedAtAtom : the atom at the center of the fragment
    \param radius : number of atoms to move out from the central atom
    \param doIsomericSmiles : include stereochemistry and isotope information
        in the SMILES
   */
  std::string GetFragmentSmiles(ROMol &mol,unsigned int rootedAtAtom,unsigned int radius,
                                bool doIsomericSmiles=false);

}
#endif
