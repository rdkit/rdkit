//
//  Copyright (C) 2003-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//! \file QueryOps.h
/*!
    \brief Includes a bunch of functionality for handling Atom and Bond queries.
*/
#ifndef _RD_QUERY_OPS_H
#define _RD_QUERY_OPS_H

#include <GraphMol/RDKitBase.h>
#include <Query/QueryObjects.h>

#ifdef RDK_THREADSAFE_SSS
#include <boost/thread/mutex.hpp>
#endif

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

  typedef Queries::GreaterQuery<int,Atom const *,true> ATOM_GREATER_QUERY;
  typedef Queries::GreaterQuery<int,Bond const *,true> BOND_GREATER_QUERY;

  typedef Queries::GreaterEqualQuery<int,Atom const *,true> ATOM_GREATEREQUAL_QUERY;
  typedef Queries::GreaterEqualQuery<int,Bond const *,true> BOND_GREATEREQUAL_QUERY;

  typedef Queries::LessQuery<int,Atom const *,true> ATOM_LESS_QUERY;
  typedef Queries::LessQuery<int,Bond const *,true> BOND_LESS_QUERY;

  typedef Queries::LessEqualQuery<int,Atom const *,true> ATOM_LESSEQUAL_QUERY;
  typedef Queries::LessEqualQuery<int,Bond const *,true> BOND_LESSEQUAL_QUERY;

  typedef Queries::RangeQuery<int,Atom const *,true> ATOM_RANGE_QUERY;
  typedef Queries::RangeQuery<int,Bond const *,true> BOND_RANGE_QUERY;

  typedef Queries::SetQuery<int,Atom const *,true> ATOM_SET_QUERY;
  typedef Queries::SetQuery<int,Bond const *,true> BOND_SET_QUERY;

  typedef Queries::Query<int,Bond const *,true> BOND_NULL_QUERY;
  typedef Queries::Query<int,Atom const *,true> ATOM_NULL_QUERY;

  // -------------------------------------------------
  // common atom queries

  static int queryAtomAromatic(Atom const * at) { return at->getIsAromatic(); };
  static int queryAtomAliphatic(Atom const * at) { return !(at->getIsAromatic()); };
  static int queryAtomExplicitDegree(Atom const * at) { return at->getDegree(); };
  static int queryAtomTotalDegree(Atom const * at) { return at->getTotalDegree(); };
  static int queryAtomHeavyAtomDegree(Atom const * at) { return at->getTotalDegree()-at->getTotalNumHs(true); };
  static int queryAtomHCount(Atom const * at) { return at->getTotalNumHs(true); };
  static int queryAtomImplicitHCount(Atom const * at) { return at->getTotalNumHs(false); };
  static int queryAtomHasImplicitH(Atom const * at) { return int(at->getTotalNumHs(false)>0); };
  static int queryAtomImplicitValence(Atom const * at) { return at->getImplicitValence(); };
  static int queryAtomExplicitValence(Atom const * at) { return at->getExplicitValence() - at->getNumExplicitHs(); };
  static int queryAtomTotalValence(Atom const * at) { return at->getExplicitValence()+at->getImplicitValence(); };
  static int queryAtomUnsaturated(Atom const * at) { return static_cast<int>(at->getDegree())<at->getExplicitValence(); };
  static int queryAtomNum(Atom const * at) { return at->getAtomicNum(); };
  static int massIntegerConversionFactor=1000;
  static int queryAtomMass(Atom const * at) {
    return static_cast<int>(round(massIntegerConversionFactor*at->getMass()));
  };
  static int queryAtomIsotope(Atom const * at) {
    return static_cast<int>(at->getIsotope());
  };
  static int queryAtomFormalCharge(Atom const * at) { 
      return static_cast<int>(at->getFormalCharge()); 
  };
  static int queryAtomHybridization(Atom const * at) { return at->getHybridization(); };
  unsigned int queryAtomBondProduct(Atom const * at);
  unsigned int queryAtomAllBondProduct(Atom const * at);
    
    
  // -------------------------------------------------
  // common bond queries

  static int queryBondOrder(Bond const * bond) { return static_cast<int>(bond->getBondType()); };
  static int queryBondDir(Bond const * bond) { return static_cast<int>(bond->getBondDir()); };
  static int queryIsBondInNRings(Bond const * at) {
    return at->getOwningMol().getRingInfo()->numBondRings(at->getIdx());
  };

  // -------------------------------------------------
  // ring queries 

  static int queryIsAtomInNRings(Atom const * at) {
    return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx());
  };
  static int queryIsAtomInRing(Atom const * at) {
    return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx())!=0;
  };
  static int queryAtomHasRingBond(Atom const * at) {
    ROMol::OBOND_ITER_PAIR atomBonds=at->getOwningMol().getAtomBonds(at);
    while(atomBonds.first != atomBonds.second){
      unsigned int bondIdx=at->getOwningMol().getTopology()[*atomBonds.first]->getIdx();
      if(at->getOwningMol().getRingInfo()->numBondRings(bondIdx)) {
        return 1;
      }
      ++atomBonds.first;  
    }
    return 0;
  };
  static int queryIsBondInRing(Bond const * bond) {
    return bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx())!=0;
  };
  static int queryAtomMinRingSize(Atom const *at){
    return at->getOwningMol().getRingInfo()->minAtomRingSize(at->getIdx());
  };
  static int queryBondMinRingSize(Bond const *bond){
    return bond->getOwningMol().getRingInfo()->minBondRingSize(bond->getIdx());
  };

  static int queryAtomRingBondCount(Atom const *at) {
    // EFF: cache this result
    int res=0;
    ROMol::OBOND_ITER_PAIR atomBonds=at->getOwningMol().getAtomBonds(at);
    while(atomBonds.first != atomBonds.second){
      unsigned int bondIdx=at->getOwningMol().getTopology()[*atomBonds.first]->getIdx();
      if(at->getOwningMol().getRingInfo()->numBondRings(bondIdx)) {
        res++;
      }
      ++atomBonds.first;  
    }
    return res;
  }

  template <int tgt>
  int queryAtomIsInRingOfSize(Atom const *at) {
    if(at->getOwningMol().getRingInfo()->isAtomInRingOfSize(at->getIdx(),tgt)){
      return tgt;
    } else {
      return 0;
    }
  };
  template <int tgt>
  int queryBondIsInRingOfSize(Bond const *bond) {
    if(bond->getOwningMol().getRingInfo()->isBondInRingOfSize(bond->getIdx(),tgt)){
      return tgt;
    } else {
      return 0;
    }
  };


  template <class T>
  T *makeAtomSimpleQuery(int what,int func(Atom const *),const std::string &description="Atom Simple"){
    T *res = new T;
    res->setVal(what);
    res->setDataFunc(func);
    res->setDescription(description);
    return res;
  }


  //! returns a Query for matching atomic number
  template <class T>
  T *makeAtomNumQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomNum,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomNumQuery(int what);

  //! returns a Query for matching implicit valence
  template <class T>
  T *makeAtomImplicitValenceQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomImplicitValence,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what);

  //! returns a Query for matching explicit valence
  template <class T>
  T *makeAtomExplicitValenceQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomExplicitValence,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomExplicitValenceQuery(int what);

  //! returns a Query for matching total valence
  template <class T>
  T *makeAtomTotalValenceQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomTotalValence,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what);

  //! returns a Query for matching explicit degree
  template <class T>
  T *makeAtomExplicitDegreeQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomExplicitDegree,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what);

  //! returns a Query for matching atomic degree
  template <class T>
  T *makeAtomTotalDegreeQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomTotalDegree,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what);

  //! returns a Query for matching hydrogen count
  template <class T>
  T *makeAtomHCountQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomHCount,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what);

  //! returns a Query for matching ring atoms
  template <class T>
  T *makeAtomHasImplicitHQuery(const std::string &descr){
    return makeAtomSimpleQuery<T>(true,queryAtomHasImplicitH,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomHasImplicitHQuery();


  //! returns a Query for matching implicit hydrogen count
  template <class T>
  T *makeAtomImplicitHCountQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomImplicitHCount,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomImplicitHCountQuery(int what);

  //! returns a Query for matching the \c isAromatic flag
  template <class T>
  T *makeAtomAromaticQuery(const std::string &descr){
    return makeAtomSimpleQuery<T>(true,queryAtomAromatic,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomAromaticQuery();

  //! returns a Query for matching aliphatic atoms
  template <class T>
  T *makeAtomAliphaticQuery(const std::string &descr){
    return makeAtomSimpleQuery<T>(true,queryAtomAliphatic,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomAliphaticQuery();

  //! returns a Query for matching atoms with a particular mass
  template <class T>
  T *makeAtomMassQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(massIntegerConversionFactor*what,queryAtomMass,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomMassQuery(int what);

  //! returns a Query for matching atoms with a particular isotope
  template <class T>
  T *makeAtomIsotopeQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomIsotope,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what);

  //! returns a Query for matching formal charge
  template <class T>
  T *makeAtomFormalChargeQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomFormalCharge,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what);

  //! returns a Query for matching hybridization
  template <class T>
  T *makeAtomHybridizationQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomHybridization,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what);

  //! returns a Query for matching atoms with unsaturation:
  template <class T>
  T *makeAtomUnsaturatedQuery(const std::string &descr){
    return makeAtomSimpleQuery<T>(true,queryAtomUnsaturated,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery();


  //! returns a Query for matching ring atoms
  template <class T>
  T *makeAtomInRingQuery(const std::string &descr){
    return makeAtomSimpleQuery<T>(true,queryIsAtomInRing,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomInRingQuery();

  //! returns a Query for matching atoms in a particular number of rings
  template <class T>
  T *makeAtomInNRingsQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryIsAtomInNRings,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what);

  //! returns a Query for matching atoms in rings of a particular size
  ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt);

  //! returns a Query for matching an atom's minimum ring size
  template <class T>
  T *makeAtomMinRingSizeQuery(int tgt,const std::string &descr){
    return makeAtomSimpleQuery<T>(tgt,queryAtomMinRingSize,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomMinRingSizeQuery(int tgt);

  //! returns a Query for matching atoms with a particular number of ring bonds
  template <class T>
  T *makeAtomRingBondCountQuery(int what,const std::string &descr){
    return makeAtomSimpleQuery<T>(what,queryAtomRingBondCount,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what);

  //! returns a Query for matching atoms with a particular number of ring bonds
  template <class T>
  T *makeAtomHasRingBondQuery(const std::string &descr){
    return makeAtomSimpleQuery<T>(1,queryAtomHasRingBond,descr);
  }
  //! \overload
  ATOM_EQUALS_QUERY *makeAtomHasRingBondQuery();

  //! returns a Query for matching bond orders
  BOND_EQUALS_QUERY *makeBondOrderEqualsQuery(Bond::BondType what);
  //! returns a Query for matching bond directions
  BOND_EQUALS_QUERY *makeBondDirEqualsQuery(Bond::BondDir what);
  //! returns a Query for matching ring bonds
  BOND_EQUALS_QUERY *makeBondIsInRingQuery();
  //! returns a Query for matching bonds in rings of a particular size
  BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int what);
  //! returns a Query for matching a bond's minimum ring size
  BOND_EQUALS_QUERY *makeBondMinRingSizeQuery(int what);
  //! returns a Query for matching bonds in a particular number of rings
  BOND_EQUALS_QUERY *makeBondInNRingsQuery(int tgt);

  //! returns a Query for matching any bond
  BOND_NULL_QUERY *makeBondNullQuery();
  //! returns a Query for matching any atom
  ATOM_NULL_QUERY *makeAtomNullQuery();

  static int queryAtomRingMembership(Atom const *at) {
    return static_cast<int>(at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx()));
  }
  // I'm pretty sure that this typedef shouldn't be necessary,
  // but VC++ generates a warning about const Atom const * in
  // the definition of Match, then complains about an override
  // that differs only by const/volatile (c4301), then generates
  // incorrect code if we don't do this... so let's do it.
  typedef Atom const *ConstAtomPtr;
  
  class AtomRingQuery : public Queries::EqualityQuery<int, ConstAtomPtr,true> {
  public:
    AtomRingQuery() : Queries::EqualityQuery<int,ConstAtomPtr,true>(-1) {
      // default is to just do a number of rings query:
      this->setDescription("AtomInNRings");
      this->setDataFunc(queryAtomRingMembership);
    };
    explicit AtomRingQuery(int v) : Queries::EqualityQuery<int,ConstAtomPtr,true>(v) {
      // default is to just do a number of rings query:
      this->setDescription("AtomInNRings");
      this->setDataFunc(queryAtomRingMembership);
    };

    virtual bool Match(const ConstAtomPtr what) const {
      int v = this->TypeConvert(what,Queries::Int2Type<true>());
      bool res;
      if(this->d_val<0){
        res = v!=0;
      } else {
        res=!Queries::queryCmp(v,this->d_val,this->d_tol);
      }
      if(this->getNegation()){
        res=!res;
      }
      return res;
    }

    //! returns a copy of this query
    Queries::Query<int,ConstAtomPtr,true> *
    copy() const {
      AtomRingQuery *res = new AtomRingQuery(this->d_val);
      res->setNegation(getNegation());
      res->setTol(this->getTol());
      res->d_description = this->d_description;
      res->d_dataFunc = this->d_dataFunc;
      return res;
    }
  };
  
  //! allows use of recursive structure queries (e.g. recursive SMARTS)
  class RecursiveStructureQuery : public Queries::SetQuery<int,Atom const *,true> {
  public:
    RecursiveStructureQuery() :
      Queries::SetQuery<int,Atom const *,true>(),
      d_serialNumber(0)
    {
      setDataFunc(getAtIdx);
      setDescription("RecursiveStructure");
    };
    //! initialize from an ROMol pointer
    /*!
      <b>Notes</b>
        - this takes over ownership of the pointer
    */
    RecursiveStructureQuery(ROMol const *query,unsigned int serialNumber=0) :
      Queries::SetQuery<int,Atom const *,true>(),
      d_serialNumber(serialNumber)
    {
      setQueryMol(query);
      setDataFunc(getAtIdx);
      setDescription("RecursiveStructure");
    };
    //! returns the index of an atom
    static int getAtIdx(Atom const *at) {
      PRECONDITION(at,"bad atom argument");
      return at->getIdx();
    };

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
      res->dp_queryMol.reset(new ROMol(*dp_queryMol,true));

      std::set<int>::const_iterator i;
      for(i=d_set.begin();i!=d_set.end();i++){
	res->insert(*i);
      }
      res->setNegation(getNegation());
      res->d_description = d_description;
      res->d_serialNumber=d_serialNumber;
      return res;
    }
    unsigned int getSerialNumber() const { return d_serialNumber; };
  
#ifdef RDK_THREADSAFE_SSS
    boost::mutex d_mutex;
#endif
  private:
    boost::shared_ptr<const ROMol>dp_queryMol;
    unsigned int d_serialNumber;
  };

  template <typename T>
  int nullDataFun(T arg) { return 1; }
  template <typename T>
  bool nullQueryFun(T arg) { return true; } 
  
};


#endif
