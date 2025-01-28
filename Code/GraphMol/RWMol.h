//
//  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
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

#include <RDGeneral/export.h>

#ifndef RD_RWMOL_H
#define RD_RWMOL_H

// our stuff
#include "ROMol.h"
#include "RingInfo.h"

namespace RDKit {

//! RWMol is a molecule class that is intended to be edited
/*!
    See documentation for ROMol for general remarks

 */
class RDKIT_GRAPHMOL_EXPORT RWMol : public ROMol {
 public:
  RWMol() : ROMol() {}
  //! copy constructor with a twist
  /*!
    \param other     the molecule to be copied
    \param quickCopy (optional) if this is true, the resulting ROMol will not
         copy any of the properties or bookmarks and conformers from \c other.
    This can
         make the copy substantially faster (thus the name).
    \param confId if this is >=0, the resulting ROMol will contain only
         the specified conformer from \c other.
  */
  RWMol(const ROMol &other, bool quickCopy = false, int confId = -1)
      : ROMol(other, quickCopy, confId) {}
  RWMol(const RWMol &other) : ROMol(other) {}
  RWMol &operator=(const RWMol &);
  RWMol(RWMol &&other) noexcept : ROMol(std::move(other)) {}
  RWMol &operator=(RWMol &&other) noexcept {
    ROMol::operator=(std::move(other));
    return *this;
  }

  //! insert the atoms and bonds from \c other into this molecule
  void insertMol(const ROMol &other);

  //! \name Atoms
  //! @{

  //! adds an empty Atom to our collection
  /*!
    \param updateLabel   (optional) if this is true, the new Atom will be
                         our \c activeAtom

    \return the index of the added atom

  */
  unsigned int addAtom(bool updateLabel = true);

  //! adds an Atom to our collection
  /*!
    \param atom          pointer to the Atom to add
    \param updateLabel   (optional) if this is true, the new Atom will be
                         our \c activeAtom
    \param takeOwnership (optional) if this is true, we take ownership of \c
    atom
                         instead of copying it.

    \return the index of the added atom
  */
  unsigned int addAtom(Atom *atom, bool updateLabel = true,
                       bool takeOwnership = false) {
    return ROMol::addAtom(atom, updateLabel, takeOwnership);
  }

  //! adds an Atom to our collection

  //! replaces a particular Atom
  /*!
    \param idx          the index of the Atom to replace
    \param atom         the new atom, which will be copied.
    \param updateLabel   (optional) if this is true, the new Atom will be
                         our \c activeAtom
    \param preserveProps if true preserve the original atom property data

  */
  void replaceAtom(unsigned int idx, Atom *atom, bool updateLabel = false,
                   bool preserveProps = false);
  //! returns a pointer to the highest-numbered Atom
  Atom *getLastAtom() { return getAtomWithIdx(getNumAtoms() - 1); }
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
  //! \overload
  void removeAtom(Atom *atom, bool clearProps);

  //! @}

  //! \name Bonds
  //! @{

  //! adds a Bond between the indicated Atoms
  /*!
     \return the number of Bonds
  */
  unsigned int addBond(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                       Bond::BondType order = Bond::UNSPECIFIED);
  //! \overload
  unsigned int addBond(Atom *beginAtom, Atom *endAtom,
                       Bond::BondType order = Bond::UNSPECIFIED);

  //! adds a Bond to our collection
  /*!
    \param bond          pointer to the Bond to add
    \param takeOwnership (optional) if this is true, we take ownership of \c
    bond
                         instead of copying it.

    \return the new number of bonds
  */
  unsigned int addBond(Bond *bond, bool takeOwnership = false) {
    return ROMol::addBond(bond, takeOwnership);
  }

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
                          Bond::BondType order = Bond::UNSPECIFIED);
  //! finishes a partially constructed bond
  /*!
    \return the final number of Bonds

    See the documentation for createPartialBond() for more details
  */
  unsigned int finishPartialBond(unsigned int endAtomIdx, int bondBookmark,
                                 Bond::BondType order = Bond::UNSPECIFIED);

  //! removes a bond from the molecule
  void removeBond(unsigned int beginAtomIdx, unsigned int endAtomIdx);

  //! replaces a particular Bond
  /*!
    \param idx          the index of the Bond to replace
    \param bond         the new bond, which will be copied.
    \param preserveProps if true preserve the original bond property data
    \param keepSGroups if true, keep Substance groups referencing the bond

  */
  void replaceBond(unsigned int idx, Bond *bond, bool preserveProps = false,
                   bool keepSGroups = true);

  //! @}

  //! removes all atoms, bonds, properties, bookmarks, etc.
  void clear() {
    destroy();
    d_confs.clear();
    ROMol::initMol();  // make sure we have a "fresh" ready to go copy
    numBonds = 0;
  }

  void beginBatchEdit();
  void rollbackBatchEdit() {
    dp_delAtoms.reset();
    dp_delBonds.reset();
  }
  void commitBatchEdit();

 private:
  void batchRemoveBonds();
  void batchRemoveAtoms();
};

typedef boost::shared_ptr<RWMol> RWMOL_SPTR;
typedef std::vector<RWMOL_SPTR> RWMOL_SPTR_VECT;

};  // namespace RDKit

#endif
