//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#ifndef __RD_RWMOL_H__
#define __RD_RWMOL_H__

// our stuff
#include "ROMol.h"

namespace RDKit{

  //! RWMol is a molecule class that is intended to be edited
  /*!
      See documentation for ROMol for general remarks

   */
  class RWMol : public ROMol {
  public:
  
    using ROMol::addAtom;
    using ROMol::addBond;

    //! insert the atoms and bonds from \c other into this molecule
    void insertMol( const ROMol &other);
  

    // --------------------------------------------
    //
    //  Atoms
    //
    // --------------------------------------------

    //! adds an empty Atom to our collection
    /*!
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom

      \return the new number of atoms

    */
    unsigned int addAtom(bool updateLabel=true);
    //! replaces a particular Atom
    /*!
      \param idx          the index of the Atom to replace
      \param atom         the new atom, which will be copied.
      \param updateLabel   (optional) if this is true, the new Atom will be
                           our \c activeAtom

    */
    void replaceAtom(unsigned int idx,Atom *atom,bool updateLabel=false);
    //! returns a pointer to the highest-numbered Atom
    GRAPH_NODE_TYPE getLastAtom() { return getAtomWithIdx(getNumAtoms()-1); };
    //! returns a pointer to the "active" Atom
    /*!
       If we have an \c activeAtom, it will be returned,
       otherwise the results of getLastAtom() will be returned.
     */
    GRAPH_NODE_TYPE getActiveAtom();
    //! sets our \c activeAtom
    void setActiveAtom(Atom *atom);
    //! \overload
    void setActiveAtom(unsigned int idx);
    //! removes an Atom from the molecule
    void removeAtom(unsigned int idx);
    //! \overload
    void removeAtom(Atom *atom);


    // --------------------------------------------
    //
    //  Bonds
    //
    // --------------------------------------------

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

  private:
    std::vector<BOND_SPTR> d_partialBonds;
    void destroy();
    ROMol &operator=(const ROMol &); // disable assignment

  };

  typedef boost::shared_ptr<RWMol>    RWMOL_SPTR;
  typedef std::vector< RWMOL_SPTR > RWMOL_SPTR_VECT;

}; // end of RDKit namespace

#endif
