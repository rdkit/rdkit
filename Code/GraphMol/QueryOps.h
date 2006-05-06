//
//  Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

//! \file QueryOps.h
/*!
    \brief Includes a bunch of functionality for handling Atom and Bond queries.
*/
#ifndef _RD_QUERY_OPS_H
#define _RD_QUERY_OPS_H

#include <GraphMol/RDKitBase.h>
#include <Query/QueryObjects.h>

namespace RDKit{
  typedef Queries::Query<bool,Atom const *,true> ATOM_BOOL_QUERY;
  typedef Queries::Query<bool,Bond const *,true> BOND_BOOL_QUERY;

  typedef Queries::AndQuery<int,Atom const *,true> ATOM_AND_QUERY;
  typedef Queries::AndQuery<int,Bond const *,true> BOND_AND_QUERY;

  typedef Queries::OrQuery<int,Atom const *,true> ATOM_OR_QUERY;
  typedef Queries::OrQuery<int,Bond const *,true> BOND_OR_QUERY;

  typedef Queries::XOrQuery<int,Atom const *,true> ATOM_XOR_QUERY;
  typedef Queries::XOrQuery<int,Bond const *,true> BOND_XOR_QUERY;

  typedef Queries::EqualityQuery<int,Atom const *,true> ATOM_EQUALS_QUERY;
  typedef Queries::EqualityQuery<int,Bond const *,true> BOND_EQUALS_QUERY;

  typedef Queries::Query<int,Bond const *,true> BOND_NULL_QUERY;
  typedef Queries::Query<int,Atom const *,true> ATOM_NULL_QUERY;

  //! returns a Query for matching atomic number
  ATOM_EQUALS_QUERY *makeAtomNumEqualsQuery(int what);
  //! returns a Query for matching implicit valence
  ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what);
  //! returns a Query for matching total valence
  ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what);
  //! returns a Query for matching explicit valence
  ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what);
  //! returns a Query for matching atomic degree
  ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what);
  //! returns a Query for matching hydrogen count
  ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what);
  //! returns a Query for matching the \c isAromatic flag
  ATOM_EQUALS_QUERY *makeAtomAromaticQuery();
  //! returns a Query for matching aliphatic atoms
  ATOM_EQUALS_QUERY *makeAtomAliphaticQuery();
  //! returns a Query for matching atoms with a particular mass (for isotopes)
  ATOM_EQUALS_QUERY *makeAtomMassQuery(int what);
  //! returns a Query for matching formal charge
  ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what);
  //! returns a Query for matching hybridization
  ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what);

  //! returns a Query for matching ring atoms
  ATOM_EQUALS_QUERY *makeAtomInRingQuery();
  //! returns a Query for matching atoms in a particular number of rings
  ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what);
  //! returns a Query for matching atoms in rings of a particular size
  ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt);

  //! returns a Query for matching bond orders
  BOND_EQUALS_QUERY *makeBondOrderEqualsQuery(Bond::BondType what);
  //! returns a Query for matching bond directions
  BOND_EQUALS_QUERY *makeBondDirEqualsQuery(Bond::BondDir what);
  //! returns a Query for matching ring bonds
  BOND_EQUALS_QUERY *makeBondIsInRingQuery();
  //! returns a Query for matching bonds in rings of a particular size
  BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int what);
  //! returns a Query for matching bonds in a particular number of rings
  BOND_EQUALS_QUERY *makeBondInNRingsQuery(int tgt);

  //! returns a Query for matching any bond
  BOND_NULL_QUERY *makeBondNullQuery();
  //! returns a Query for matching any atom
  ATOM_NULL_QUERY *makeAtomNullQuery();


  //! allows use of recursive structure queries (e.g. recursive SMARTS)
  class RecursiveStructureQuery : public Queries::SetQuery<int,Atom const *,true> {
  public:
    RecursiveStructureQuery() : Queries::SetQuery<int,Atom const *,true>() {
      setDataFunc(getAtIdx);
      setDescription("RecursiveStructure");
    };
    //! initialize from an ROMol pointer
    /*!
      <b>Notes</b>
        - this takes over ownership of the pointer
    */
    RecursiveStructureQuery(ROMol const *query) : Queries::SetQuery<int,Atom const *,true>() {
      setQueryMol(query);
      setDataFunc(getAtIdx);
      setDescription("RecursiveStructure");
    };
    //! returns the index of an atom
    static int getAtIdx(Atom const *at) { return at->getIdx(); };

    //! sets the molecule we'll use recursively
    /*!
      <b>Notes</b>
        - this takes over ownership of the pointer
    */
    void setQueryMol(ROMol const *query) {
      dp_queryMol.reset(query);
    }
    //! returns a pointer to our query molecule
    ROMol const * getQueryMol() const { return dp_queryMol.get(); };

    //! returns a copy of this query
    Queries::Query<int,Atom const *,true> *
    copy() const {
      RecursiveStructureQuery *res =
	new RecursiveStructureQuery();
      res->dp_queryMol = dp_queryMol;

      std::set<int>::const_iterator i;
      for(i=d_set.begin();i!=d_set.end();i++){
	res->insert(*i);
      }
      res->setNegation(getNegation());
      res->d_description = d_description;
      return res;
    }
  
  private:
    boost::shared_ptr<const ROMol>dp_queryMol;
  
  };


  
};


#endif
