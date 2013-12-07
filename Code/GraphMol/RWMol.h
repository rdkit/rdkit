//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file RWMol.h

  \brief Defines the editable molecule class \c RWMol

*/  

#ifndef __RD_RWMOL_H__
#define __RD_RWMOL_H__

// our stuff
#include "ROMol.h"
#include "RingInfo.h"

namespace RDKit{

  //! RWMol is a molecule class that is intended to be edited
  /*!
      See documentation for ROMol for general remarks

   */
  class RWMol : public ROMol {
  public:
  
    RWMol() { d_partialBonds.clear(); }

    //! copy constructor with a twist
    /*!
      \param other     the molecule to be copied
      \param quickCopy (optional) if this is true, the resulting ROMol will not
           copy any of the properties or bookmarks and conformers from \c other.  This can
           make the copy substantially faster (thus the name).
    */
    RWMol(const ROMol &other,bool quickCopy=false) {d_partialBonds.clear(); initFromOther(other,quickCopy);};
    RWMol &operator=(const RWMol &);

    //! insert the atoms and bonds from \c other into this molecule
    void insertMol( const ROMol &other);
  

    //! \name Atoms
    //@{

    //! adds an empty Atom to our collection
    /*!
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom

      \return the new number of atoms

    */
    unsigned int addAtom(bool updateLabel=true);

    //! adds an Atom to our collection
    /*!
      \param atom          pointer to the Atom to add
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom
      \param takeOwnership (optional) if this is true, we take ownership of \c atom
                           instead of copying it.

      \return the new number of atoms
    */
    unsigned int addAtom(Atom *atom,bool updateLabel=true,bool takeOwnership=false){
      return ROMol::addAtom(atom,updateLabel,takeOwnership);
    };

    //! adds an Atom to our collection
    /*!
      \param atom          pointer to the Atom to add
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom


      \return the new number of atoms

      <b>Note:</b> since this is using a smart pointer, we don't need to worry about
      issues of ownership.

    */
    unsigned int addAtom(ATOM_SPTR atom,bool updateLabel=true){
      return ROMol::addAtom(atom,updateLabel);
    };

    //! replaces a particular Atom
    /*!
      \param idx          the index of the Atom to replace
      \param atom         the new atom, which will be copied.
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom

    */
    void replaceAtom(unsigned int idx,Atom *atom,bool updateLabel=false);
    //! returns a pointer to the highest-numbered Atom
    Atom *getLastAtom() { return getAtomWithIdx(getNumAtoms()-1); };
    //! returns a pointer to the "active" Atom
    /*!
       If we have an \c activeAtom, it will be returned,
       otherwise the results of getLastAtom() will be returned.
     */
    Atom *getActiveAtom();
    //! sets our \c activeAtom
    void setActiveAtom(Atom *atom);
    //! \overload
    void setActiveAtom(unsigned int idx);
    //! removes an Atom from the molecule
    void removeAtom(unsigned int idx);
    //! \overload
    void removeAtom(Atom *atom);

    //@}


    //! \name Bonds
    //@{

    //! adds a Bond between the indicated Atoms
    /*!
       \return the number of Bonds
    */
    unsigned int addBond(unsigned int beginAtomIdx,unsigned int endAtomIdx,
			 Bond::BondType order=Bond::UNSPECIFIED);
    //! \overload
    unsigned int addBond(ATOM_SPTR beginAtom,ATOM_SPTR endAtom,
			 Bond::BondType order=Bond::UNSPECIFIED);
    //! \overload
    unsigned int addBond(Atom *beginAtom, Atom *endAtom,
			 Bond::BondType order=Bond::UNSPECIFIED);


    //! adds a Bond to our collection
    /*!
      \param bond          pointer to the Bond to add
      \param takeOwnership (optional) if this is true, we take ownership of \c bond
                           instead of copying it.

      \return the new number of bonds
    */
    unsigned int addBond(Bond *bond,bool takeOwnership=false){
      return ROMol::addBond(bond,takeOwnership);
    };
    //! adds a Bond to our collection
    /*!
      \param bsp        smart pointer to the Bond to add

      \return the new number of bonds

      <b>Note:</b> since this is using a smart pointer, we don't need to worry about
      issues of ownership.
    */
    unsigned int addBond(BOND_SPTR bsp){
      return ROMol::addBond(bsp);
    };


    //! starts a Bond and sets its beginAtomIdx
    /*!
      \return a pointer to the new bond

      The caller should set a bookmark to the returned Bond in order
      to be able to later complete it:

      \verbatim
        Bond *pBond = mol->createPartialBond(1);
	mol->setBondBookmark(pBond,666);
	... do some other stuff ...
	mol->finishPartialBond(2,666,Bond::SINGLE);
	mol->clearBondBookmark(666,pBond);
      \endverbatim

      or, if we want to set the \c BondType initially:
      \verbatim
        Bond *pBond = mol->createPartialBond(1,Bond::DOUBLE);
	mol->setBondBookmark(pBond,666);
	... do some other stuff ...
	mol->finishPartialBond(2,666);
	mol->clearBondBookmark(666,pBond);
      \endverbatim

      the call to finishPartialBond() will take priority if you set the
      \c BondType in both calls.
      
    */
    Bond *createPartialBond(unsigned int beginAtomIdx,
			    Bond::BondType order=Bond::UNSPECIFIED);
    //! finishes a partially constructed bond
    /*!
      \return the final number of Bonds

      See the documentation for createPartialBond() for more details
    */
    unsigned int finishPartialBond(unsigned int endAtomIdx,int bondBookmark,
				   Bond::BondType order=Bond::UNSPECIFIED);

    //! removes a bond from the molecule
    void removeBond(unsigned int beginAtomIdx, unsigned int endAtomIdx);
    //@}

    //! removes all atoms, bonds, properties, bookmarks, etc.
    void clear() {
      d_atomBookmarks.clear();
      d_bondBookmarks.clear();
      d_graph.clear();
      d_confs.clear();
      if(dp_props){
        dp_props->reset();
        STR_VECT computed;
        dp_props->setVal(detail::computedPropName, computed);
      }
      if(dp_ringInfo) dp_ringInfo->reset();
    };


  private:
    std::vector<BOND_SPTR> d_partialBonds;
    void destroy();

  };

  typedef boost::shared_ptr<RWMol>    RWMOL_SPTR;
  typedef std::vector< RWMOL_SPTR > RWMOL_SPTR_VECT;

}; // end of RDKit namespace

#endif
