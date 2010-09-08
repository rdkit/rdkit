//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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
   */
  std::string MolToSmiles(ROMol &mol,bool doIsomericSmiles=false,
			  bool doKekule=false,int rootedAtAtom=-1);
}
#endif
