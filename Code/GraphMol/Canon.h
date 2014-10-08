//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_CANON_H_
#define _RD_CANON_H_

#include <boost/tuple/tuple.hpp>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {
  class ROMol;
  class Atom;
  class Bond;
};

namespace Canon {
  const int MAX_NATOMS=5000; //!< used in the canonical traversal code
  const int MAX_CYCLES=1000;   //!< used in the canonical traversal code
  const int MAX_BONDTYPE=32; //!< used in the canonical traversal code

  //! used in traversals of the molecule
  typedef enum {
    WHITE_NODE=0,  //! not visited
    GREY_NODE,     //! visited, but not finished
    BLACK_NODE,    //! visited and finished
  } AtomColors; 

  //! used to indicate types of entries in the molecular stack:
  typedef enum {
    MOL_STACK_ATOM=0,       //!< an Atom
    MOL_STACK_BOND,         //!< a Bond
    MOL_STACK_RING,         //!< a ring closure
    MOL_STACK_BRANCH_OPEN,  //!< beginning of a branch
    MOL_STACK_BRANCH_CLOSE, //!< end of a branch
  } MolStackTypes;

  //! used to store components in the molecular stack
  typedef union{
    RDKit::Atom *atom;
    RDKit::Bond *bond;
  } MolStackUnion;

  //! these are the actual elements in the molecular stack
  class MolStackElem {
  public:
    //! construct an Atom node
    explicit MolStackElem(RDKit::Atom *at) {
      type = MOL_STACK_ATOM;
      obj.atom = at;
    };
    //! construct a bond node
    /*!

       \param bond  pointer to the Bond being added
       \param idx   index of the Atom traversed before this Bond
         (beginAtom in the canonical traversal order)
    */
    explicit MolStackElem(RDKit::Bond *bond,int idx) {
      type = MOL_STACK_BOND;
      obj.bond = bond;
      number = idx;
    };
    //! construct for a ring closure
    explicit MolStackElem(int idx) {
      type = MOL_STACK_RING;
      number = idx;
    };
    //! construct for a branch opening or closing
    explicit MolStackElem(const char *chr,int idx) {
      switch(chr[0]){
      case '(':
	type = MOL_STACK_BRANCH_OPEN;
	break;
      case ')':
	type = MOL_STACK_BRANCH_CLOSE;
	break;
      default:
	break;
      }
      number=idx;
    }
    MolStackTypes type; //!< stores the type of node
    MolStackUnion obj;  //!< holds our pointer (if appropriate)
    int number;         //!< stores our number (relevant for bonds and ring closures)
  };
  typedef std::vector<MolStackElem> MolStack;


  //! used to represent possible branches from an atom
  typedef boost::tuple<int,int,RDKit::Bond *> PossibleType;

  //! constructs the canonical traversal order for a molecular fragment
  /*!

    \param mol       the ROMol we're working on
    \param atomIdx   the index of the atom to start the traversal from
    \param colors    the traversal status of each atom in \c mol
    \param ranks     the assigned rank of each atom in \c mol
    \param molStack  the current traversal stack (used to return the results)

    <b>Notes</b>
      - \c mol will, in general, be modified by this operation as bond directions
        and the like are changed to fit the canonical traversal order

   */
  void canonicalizeFragment(RDKit::ROMol &mol,int atomIdx,
			    std::vector<AtomColors> &colors,
			    std::vector<int> &ranks,
			    MolStack &molStack,
                            const boost::dynamic_bitset<> *bondsInPlay=0,
                            const std::vector<std::string> *bondSymbols=0
                            );

};

#endif
