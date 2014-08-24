// $Id$
//
//  Copyright (C) 2001-2013 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/MonomerInfo.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/types.h>
#include <Query/QueryObjects.h>
#include <map>
#include <boost/cstdint.hpp>
using boost::int32_t;
using boost::uint32_t;
namespace RDKit{

  const int32_t MolPickler::versionMajor=7;
  const int32_t MolPickler::versionMinor=2;
  const int32_t MolPickler::versionPatch=0;
  const int32_t MolPickler::endianId=0xDEADBEEF;

  void streamWrite(std::ostream &ss,const std::string &what){
    unsigned int l=what.length();
    ss.write((const char *)&l,sizeof(l));
    ss.write(what.c_str(),sizeof(char)*l);
  };  

  void streamWrite(std::ostream &ss,MolPickler::Tags tag){
    unsigned char tmp=static_cast<unsigned char>(tag);
    streamWrite(ss,tmp);
  }
  template <typename T>
  void streamWrite(std::ostream &ss,MolPickler::Tags tag,const T &what){
    streamWrite(ss,tag);
    streamWrite(ss,what);
  };  

  template<class T>
  void streamRead(std::istream &ss,T &obj,int version){
    streamRead(ss,obj);
  }    

  void streamRead(std::istream &ss,std::string &what,int version){
    unsigned int l;
    ss.read((char *)&l,sizeof(l));
    char *buff=new char[l+1];
    ss.read(buff,sizeof(char)*l);
    buff[l]=0;
    what=buff;
    delete [] buff;
  };  
  void streamRead(std::istream &ss,MolPickler::Tags &tag,int version){
    if(version<7000){
      int32_t tmp;
      streamRead(ss,tmp,version);
      tag=static_cast<MolPickler::Tags>(tmp);
    } else {
      unsigned char tmp;
      streamRead(ss,tmp,version);
      tag=static_cast<MolPickler::Tags>(tmp);
    }
  }

  namespace {
    using namespace Queries;
    template <class T>
    void pickleQuery(std::ostream &ss,const Query<int,T const *,true> *query) {
      PRECONDITION(query,"no query");
      streamWrite(ss,query->getDescription());
      if(query->getNegation()) streamWrite(ss,MolPickler::QUERY_ISNEGATED);
      int32_t queryVal;
      //if (typeid(*query)==typeid(ATOM_BOOL_QUERY)){
      //  streamWrite(ss,QUERY_BOOL);
      if (typeid(*query)==typeid(AndQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_AND);
      } else if (typeid(*query)==typeid(OrQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_OR);
      } else if (typeid(*query)==typeid(XOrQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_XOR);
      } else if (typeid(*query)==typeid(EqualityQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_EQUALS);
        queryVal=static_cast<const EqualityQuery<int,T const *,true>*>(query)->getVal();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        queryVal=static_cast<const EqualityQuery<int,T const *,true>*>(query)->getTol();
        streamWrite(ss,queryVal);
      } else if (typeid(*query)==typeid(GreaterQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_GREATER);
        queryVal=static_cast<const GreaterQuery<int,T const *,true>*>(query)->getVal();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        queryVal=static_cast<const GreaterQuery<int,T const *,true>*>(query)->getTol();
        streamWrite(ss,queryVal);
      } else if (typeid(*query)==typeid(GreaterEqualQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_GREATEREQUAL);
        queryVal=static_cast<const GreaterEqualQuery<int,T const *,true>*>(query)->getVal();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        queryVal=static_cast<const GreaterEqualQuery<int,T const *,true>*>(query)->getTol();
        streamWrite(ss,queryVal);
      } else if (typeid(*query)==typeid(LessQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_LESS);
        queryVal=static_cast<const LessQuery<int,T const *,true>*>(query)->getVal();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        queryVal=static_cast<const LessQuery<int,T const *,true>*>(query)->getTol();
        streamWrite(ss,queryVal);
      } else if (typeid(*query)==typeid(LessEqualQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_LESSEQUAL);
        queryVal=static_cast<const LessEqualQuery<int,T const *,true>*>(query)->getVal();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        queryVal=static_cast<const LessEqualQuery<int,T const *,true>*>(query)->getTol();
        streamWrite(ss,queryVal);
      } else if (typeid(*query)==typeid(RangeQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_RANGE);
        queryVal=static_cast<const RangeQuery<int,T const *,true>*>(query)->getLower();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        queryVal=static_cast<const RangeQuery<int,T const *,true>*>(query)->getUpper();
        streamWrite(ss,queryVal);
        queryVal=static_cast<const RangeQuery<int,T const *,true>*>(query)->getTol();
        streamWrite(ss,queryVal);
        char ends;
        bool lowerOpen,upperOpen;
        boost::tie(lowerOpen,upperOpen)=static_cast<const RangeQuery<int,T const *,true>*>(query)->getEndsOpen();
        ends=0|(lowerOpen<<1)|upperOpen;
        streamWrite(ss,ends);
      } else if (typeid(*query)==typeid(SetQuery<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_SET);
        queryVal=static_cast<const SetQuery<int,T const *,true>*>(query)->size();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        typename SetQuery<int,T const *,true>::CONTAINER_TYPE::const_iterator cit;
        for(cit=static_cast<const SetQuery<int,T const *,true>*>(query)->beginSet();
            cit!=static_cast<const SetQuery<int,T const *,true>*>(query)->endSet();
            ++cit){
          queryVal=*cit;
          streamWrite(ss,queryVal);
        }
      } else if (typeid(*query)==typeid(AtomRingQuery)){
        streamWrite(ss,MolPickler::QUERY_ATOMRING);
        queryVal=static_cast<const EqualityQuery<int,T const *,true>*>(query)->getVal();
        streamWrite(ss,MolPickler::QUERY_VALUE,queryVal);
        queryVal=static_cast<const EqualityQuery<int,T const *,true>*>(query)->getTol();
        streamWrite(ss,queryVal);
      } else if (typeid(*query)==typeid(RecursiveStructureQuery)){
        streamWrite(ss,MolPickler::QUERY_RECURSIVE);
        streamWrite(ss,MolPickler::QUERY_VALUE);
        MolPickler::pickleMol(((const RecursiveStructureQuery *)query)->getQueryMol(),ss);
      } else if (typeid(*query)==typeid(Query<int,T const *,true>)){
        streamWrite(ss,MolPickler::QUERY_NULL);
      } else {
        throw MolPicklerException("do not know how to pickle part of the query.");
      }

      // now the children:
      streamWrite(ss,MolPickler::QUERY_NUMCHILDREN,
                  static_cast<unsigned char>(query->endChildren()-query->beginChildren()));
      typename Query<int,T const *,true>::CHILD_VECT_CI cit;
      for(cit=query->beginChildren();cit!=query->endChildren();++cit){
        pickleQuery(ss,cit->get());
      }
    }

    void finalizeQueryFromDescription(Query<int,Atom const *,true> *query,
                                      Atom const *owner){
      std::string descr=query->getDescription();
      Query<int,Atom const *,true> *tmpQuery;
      if(descr=="AtomRingBondCount"){
        query->setDataFunc(queryAtomRingBondCount);
      } else if(descr=="AtomHasRingBond"){
        query->setDataFunc(queryAtomHasRingBond);
      } else if(descr=="AtomRingSize"){
        tmpQuery=makeAtomInRingOfSizeQuery(static_cast<ATOM_EQUALS_QUERY *>(query)->getVal());
        query->setDataFunc(tmpQuery->getDataFunc());
        delete tmpQuery;
      } else if(descr=="AtomMinRingSize"){
        query->setDataFunc(queryAtomMinRingSize);
      } else if(descr=="AtomRingBondCount"){
        query->setDataFunc(queryAtomRingBondCount);
      } else if(descr=="AtomImplicitValence"){
        query->setDataFunc(queryAtomImplicitValence);
      } else if(descr=="AtomTotalValence"){
        query->setDataFunc(queryAtomTotalValence);
      } else if(descr=="AtomAtomicNum"){
        query->setDataFunc(queryAtomNum);
      } else if(descr=="AtomExplicitDegree"){
        query->setDataFunc(queryAtomExplicitDegree);
      } else if(descr=="AtomTotalDegree"){
        query->setDataFunc(queryAtomTotalDegree);
      } else if(descr=="AtomHCount"){
        query->setDataFunc(queryAtomHCount);
      } else if(descr=="AtomImplicitHCount"){
        query->setDataFunc(queryAtomImplicitHCount);
      } else if(descr=="AtomHasImplicitH"){
        query->setDataFunc(queryAtomHasImplicitH);
      } else if(descr=="AtomIsAromatic"){
        query->setDataFunc(queryAtomAromatic);
      } else if(descr=="AtomIsAliphatic"){
        query->setDataFunc(queryAtomAliphatic);
      } else if(descr=="AtomUnsaturated"){
        query->setDataFunc(queryAtomUnsaturated);
      } else if(descr=="AtomMass"){
        query->setDataFunc(queryAtomMass);
      } else if(descr=="AtomIsotope"){
        query->setDataFunc(queryAtomIsotope);
      } else if(descr=="AtomFormalCharge"){
        query->setDataFunc(queryAtomFormalCharge);
      } else if(descr=="AtomHybridization"){
        query->setDataFunc(queryAtomHybridization);
      } else if(descr=="AtomInRing"){
        query->setDataFunc(queryIsAtomInRing);
      } else if(descr=="AtomInNRings"){
        query->setDataFunc(queryIsAtomInNRings);
      } else if(descr=="AtomNull"){
        query->setDataFunc(nullDataFun);
        query->setMatchFunc(nullQueryFun);
      } else if(descr=="AtomInNRings"||descr=="RecursiveStructure"){
        // don't need to do anything here because the classes
        // automatically have everything set
      } else if(descr=="AtomAnd"||descr=="AtomOr"||descr=="AtomXor"){
        // don't need to do anything here because the classes
        // automatically have everything set
      } else {
        throw MolPicklerException("Do not know how to finalize query: '"+descr+"'");
      }
    }

    void finalizeQueryFromDescription(Query<int,Bond const *,true> *query,
                                      Bond const *owner){
      std::string descr=query->getDescription();
      Query<int,Bond const *,true> *tmpQuery;
      if(descr=="BondRingSize"){
        tmpQuery=makeBondInRingOfSizeQuery(static_cast<BOND_EQUALS_QUERY *>(query)->getVal());
        query->setDataFunc(tmpQuery->getDataFunc());
        delete tmpQuery;
      } else if(descr=="BondMinRingSize"){
        query->setDataFunc(queryBondMinRingSize);
      } else if(descr=="BondOrder"){
        query->setDataFunc(queryBondOrder);
      } else if(descr=="BondDir"){
        query->setDataFunc(queryBondDir);
      } else if(descr=="BondInRing"){
        query->setDataFunc(queryIsBondInRing);
      } else if(descr=="BondInNRings"){
        query->setDataFunc(queryIsBondInNRings);
      } else if(descr=="BondNull"){
        query->setDataFunc(nullDataFun);
        query->setMatchFunc(nullQueryFun);
      } else if(descr=="BondAnd"||descr=="BondOr"||descr=="BondXor"){
        // don't need to do anything here because the classes
        // automatically have everything set
      } else {
        throw MolPicklerException("Do not know how to finalize query: '"+descr+"'");
      }
      
    }
    
    template <class T>
    Query<int,T const *,true> *buildBaseQuery(std::istream &ss,
                                              T const *owner,MolPickler::Tags tag,
                                              int version) {
      PRECONDITION(owner,"no query");
      std::string descr;
      Query<int,T const *,true> *res=0;
      int32_t val;
      int32_t nMembers;
      switch(tag){
      case MolPickler::QUERY_AND:
        res = new AndQuery<int,T const *,true>();
        break;
      case MolPickler::QUERY_OR:
        res = new OrQuery<int,T const *,true>();
        break;
      case MolPickler::QUERY_XOR:
        res = new XOrQuery<int,T const *,true>();
        break;
      case MolPickler::QUERY_EQUALS:
        res = new EqualityQuery<int,T const *,true>();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,val,version);
        static_cast<EqualityQuery<int,T const *,true> *>(res)->setVal(val);
        streamRead(ss,val,version);
        static_cast<EqualityQuery<int,T const *,true> *>(res)->setTol(val);
        break;
      case MolPickler::QUERY_GREATER:
        res = new GreaterQuery<int,T const *,true>();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,val,version);
        static_cast<GreaterQuery<int,T const *,true> *>(res)->setVal(val);
        streamRead(ss,val,version);
        static_cast<GreaterQuery<int,T const *,true> *>(res)->setTol(val);
        break;
      case MolPickler::QUERY_GREATEREQUAL:
        res = new GreaterEqualQuery<int,T const *,true>();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,val,version);
        static_cast<GreaterEqualQuery<int,T const *,true> *>(res)->setVal(val);
        streamRead(ss,val,version);
        static_cast<GreaterEqualQuery<int,T const *,true> *>(res)->setTol(val);
        break;
      case MolPickler::QUERY_LESS:
        res = new LessQuery<int,T const *,true>();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,val,version);
        static_cast<LessQuery<int,T const *,true> *>(res)->setVal(val);
        streamRead(ss,val,version);
        static_cast<LessQuery<int,T const *,true> *>(res)->setTol(val);
        break;
      case MolPickler::QUERY_LESSEQUAL:
        res = new LessEqualQuery<int,T const *,true>();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,val,version);
        static_cast<LessEqualQuery<int,T const *,true> *>(res)->setVal(val);
        streamRead(ss,val,version);
        static_cast<LessEqualQuery<int,T const *,true> *>(res)->setTol(val);
        break;
      case MolPickler::QUERY_RANGE:
        res = new RangeQuery<int,T const *,true>();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,val,version);
        static_cast<RangeQuery<int,T const *,true> *>(res)->setLower(val);
        streamRead(ss,val,version);
        static_cast<RangeQuery<int,T const *,true> *>(res)->setUpper(val);
        streamRead(ss,val,version);
        static_cast<RangeQuery<int,T const *,true> *>(res)->setTol(val);
        break;
      case MolPickler::QUERY_SET:
        res = new SetQuery<int,T const *,true>();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,nMembers);
        while(nMembers>0){
          streamRead(ss,val,version);
          static_cast<SetQuery<int,T const *,true> *>(res)->insert(val);
          --nMembers;
        }
        break;
      case MolPickler::QUERY_NULL:
        res = new Query<int,T const *,true>();
        break;
      default:
        throw MolPicklerException("unknown query-type tag encountered");
      }

      POSTCONDITION(res,"no match found");
      return res;
    }
    
    Query<int,Atom const *,true> *unpickleQuery(std::istream &ss,Atom const *owner,
                                                int version) {
      PRECONDITION(owner,"no query");
      std::string descr;
      bool isNegated=false;
      Query<int,Atom const *,true> *res;      
      streamRead(ss,descr,version);
      MolPickler::Tags tag;
      streamRead(ss,tag,version);
      if(tag==MolPickler::QUERY_ISNEGATED){
        isNegated=true;
        streamRead(ss,tag,version);
      }
      int32_t val;
      ROMol *tmpMol;
      switch(tag){
      case MolPickler::QUERY_ATOMRING:
        res=new AtomRingQuery();
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        streamRead(ss,val,version);
        static_cast<EqualityQuery<int,Atom const *,true> *>(res)->setVal(val);
        streamRead(ss,val,version);
        static_cast<EqualityQuery<int,Atom const *,true> *>(res)->setTol(val);
        break;
      case MolPickler::QUERY_RECURSIVE:
        streamRead(ss,tag,version);
        if(tag != MolPickler::QUERY_VALUE){
          throw MolPicklerException("Bad pickle format: QUERY_VALUE tag not found.");
        }
        tmpMol=new ROMol();
        MolPickler::molFromPickle(ss,tmpMol);
        res=new RecursiveStructureQuery(tmpMol);
        break;
      default:
        res = buildBaseQuery(ss,owner,tag,version);
        break;
      }
      CHECK_INVARIANT(res,"no query!");
      
      res->setNegation(isNegated);
      res->setDescription(descr);
      
      finalizeQueryFromDescription(res,owner);

      // read in the children:
      streamRead(ss,tag,version);
      if(tag != MolPickler::QUERY_NUMCHILDREN){
        throw MolPicklerException("Bad pickle format: QUERY_NUMCHILDREN tag not found.");
      }
      unsigned char numChildren;
      streamRead(ss,numChildren,version);
      while(numChildren>0){
        Query<int,Atom const*,true> *child=unpickleQuery(ss,owner,version);
        res->addChild(Query<int,Atom const *,true>::CHILD_TYPE(child));
        --numChildren;
      }
      return res;
    }
    Query<int,Bond const *,true> *unpickleQuery(std::istream &ss,Bond const *owner,
                                                int version) {
      PRECONDITION(owner,"no query");
      std::string descr;
      bool isNegated=false;
      Query<int,Bond const *,true> *res;      
      streamRead(ss,descr,version);
      MolPickler::Tags tag;
      streamRead(ss,tag,version);
      if(tag==MolPickler::QUERY_ISNEGATED){
        isNegated=true;
        streamRead(ss,tag,version);
      }
      res = buildBaseQuery(ss,owner,tag,version);
      CHECK_INVARIANT(res,"no query!");
      
      res->setNegation(isNegated);
      res->setDescription(descr);
      
      finalizeQueryFromDescription(res,owner);

      // read in the children:
      streamRead(ss,tag,version);
      if(tag != MolPickler::QUERY_NUMCHILDREN){
        throw MolPicklerException("Bad pickle format: QUERY_NUMCHILDREN tag not found.");
      }
      unsigned char numChildren;
      streamRead(ss,numChildren,version);
      while(numChildren>0){
        Query<int,Bond const*,true> *child=unpickleQuery(ss,owner,version);
        res->addChild(Query<int,Bond const *,true>::CHILD_TYPE(child));
        --numChildren;
      }
      return res;
    }

    void pickleAtomPDBResidueInfo(std::ostream &ss,const AtomPDBResidueInfo *info){
      PRECONDITION(info,"no info");
      if(info->getSerialNumber())
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_SERIALNUMBER,info->getSerialNumber());
      if(info->getAltLoc()!="")
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_ALTLOC,info->getAltLoc());
      if(info->getResidueName()!="")
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_RESIDUENAME,info->getResidueName());
      if(info->getResidueNumber())
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_RESIDUENUMBER,info->getResidueNumber());
      if(info->getChainId()!="")
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_CHAINID,info->getChainId());
      if(info->getInsertionCode()!="")
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_INSERTIONCODE,info->getInsertionCode());
      if(info->getOccupancy())
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_OCCUPANCY,info->getOccupancy());
      if(info->getTempFactor())
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_TEMPFACTOR,info->getTempFactor());
      if(info->getIsHeteroAtom())
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_ISHETEROATOM,static_cast<char>(info->getIsHeteroAtom()));
      if(info->getSecondaryStructure())
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE,info->getSecondaryStructure());
      if(info->getSegmentNumber())
        streamWrite(ss,MolPickler::ATOM_PDB_RESIDUE_SEGMENTNUMBER,info->getSegmentNumber());
    }

    void unpickleAtomPDBResidueInfo(std::istream &ss,AtomPDBResidueInfo *info,
                                        int version){
      PRECONDITION(info,"no info");
      std::string sval;
      double dval;
      char cval;
      unsigned int uival;
      int ival;
      MolPickler::Tags tag=MolPickler::BEGIN_ATOM_MONOMER;
      while(tag!=MolPickler::END_ATOM_MONOMER){
        streamRead(ss,tag,version);
        switch(tag){
        case MolPickler::ATOM_PDB_RESIDUE_SERIALNUMBER:
          streamRead(ss,ival,version);
          info->setSerialNumber(ival);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_ALTLOC:
          streamRead(ss,sval,version);
          info->setAltLoc(sval);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_RESIDUENAME:
          streamRead(ss,sval,version);
          info->setResidueName(sval);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_RESIDUENUMBER:
          streamRead(ss,ival,version);
          info->setResidueNumber(ival);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_CHAINID:
          streamRead(ss,sval,version);
          info->setChainId(sval);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_INSERTIONCODE:
          streamRead(ss,sval,version);
          info->setInsertionCode(sval);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_OCCUPANCY:
          streamRead(ss,dval,version);
          info->setOccupancy(dval);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_TEMPFACTOR:
          streamRead(ss,dval,version);
          info->setTempFactor(dval);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_ISHETEROATOM:
          streamRead(ss,cval,version);
          info->setIsHeteroAtom(cval);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE:
          streamRead(ss,uival,version);
          info->setSecondaryStructure(uival);
          break;
        case MolPickler::ATOM_PDB_RESIDUE_SEGMENTNUMBER:
          streamRead(ss,uival,version);
          info->setSegmentNumber(uival);
          break;
        case MolPickler::END_ATOM_MONOMER:
          break;
        default:
          throw MolPicklerException("unrecognized tag while parsing atom peptide residue info");
        }
      }
    }

    
    void pickleAtomMonomerInfo(std::ostream &ss,const AtomMonomerInfo *info) {
      PRECONDITION(info,"no info");
      streamWrite(ss,info->getName());
      streamWrite(ss,static_cast<unsigned int>(info->getMonomerType()));
      switch(info->getMonomerType()){
      case AtomMonomerInfo::UNKNOWN:
      case AtomMonomerInfo::OTHER:
        break;
      case AtomMonomerInfo::PDBRESIDUE:
        pickleAtomPDBResidueInfo(ss,static_cast<const AtomPDBResidueInfo *>(info));
        break;
      default:
        throw MolPicklerException("unrecognized MonomerType");
      }
    }
    AtomMonomerInfo *unpickleAtomMonomerInfo(std::istream &ss,int version) {
      MolPickler::Tags tag;
      std::string nm;
      streamRead(ss,nm,version);
      unsigned int typ;
      streamRead(ss,typ,version);

      AtomMonomerInfo *res;
      switch(typ){
      case AtomMonomerInfo::UNKNOWN:
      case AtomMonomerInfo::OTHER:
        res = new AtomMonomerInfo(RDKit::AtomMonomerInfo::AtomMonomerType(typ),nm);
        streamRead(ss,tag,version);
        if(tag!=MolPickler::END_ATOM_MONOMER)
          throw MolPicklerException("did not find expected end of atom monomer info");          
        break;
      case AtomMonomerInfo::PDBRESIDUE:
        res = static_cast<AtomMonomerInfo *>(new AtomPDBResidueInfo(nm));
        unpickleAtomPDBResidueInfo(ss,static_cast<AtomPDBResidueInfo *>(res),version);
        break;
      default:
        throw MolPicklerException("unrecognized MonomerType");
      }
      return res;
    }


  } // end of anonymous namespace


  void MolPickler::pickleMol(const ROMol *mol,std::ostream &ss){
    PRECONDITION(mol,"empty molecule");
    streamWrite(ss,endianId);
    streamWrite(ss,static_cast<int>(VERSION));
    streamWrite(ss,versionMajor);
    streamWrite(ss,versionMinor);
    streamWrite(ss,versionPatch);
#ifndef OLD_PICKLE
    if(mol->getNumAtoms()>255){
      _pickle<int32_t>(mol,ss);
    } else {
      _pickle<unsigned char>(mol,ss);
    }
#else
    _pickleV1(mol,ss);
#endif    
  }
  void MolPickler::pickleMol(const ROMol *mol,std::string &res){
    PRECONDITION(mol,"empty molecule");
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    MolPickler::pickleMol(mol,ss);
    res = ss.str();
  }

  // NOTE: if the mol passed in here already has atoms and bonds, they will
  // be left intact.  The side effect is that ALL atom and bond bookmarks
  // will be blown out by the end of this process.
  void MolPickler::molFromPickle(std::istream &ss,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    int32_t tmpInt;

    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();
  
    streamRead(ss,tmpInt);
    if(tmpInt!=endianId){
      throw MolPicklerException("Bad pickle format: bad endian ID or invalid file format");
    }
  
    streamRead(ss,tmpInt);
    if(static_cast<Tags>(tmpInt)!=VERSION){
      throw MolPicklerException("Bad pickle format: no version tag");
    }
    int32_t majorVersion,minorVersion,patchVersion;
    streamRead(ss,majorVersion);
    streamRead(ss,minorVersion);
    streamRead(ss,patchVersion);
    if(majorVersion>versionMajor||(majorVersion==versionMajor&&minorVersion>versionMinor)){
      BOOST_LOG(rdWarningLog)<<"Depickling from a version number ("<<majorVersion<<"." << minorVersion<<")" << "that is higher than our version ("<<versionMajor<<"."<<versionMinor<<").\nThis probably won't work."<<std::endl;
    }
    majorVersion=1000*majorVersion+minorVersion*10+patchVersion;
    if(majorVersion==1){
      _depickleV1(ss,mol);
    } else {
      int32_t numAtoms;
      streamRead(ss,numAtoms,majorVersion);
      if(numAtoms>255){
	_depickle<int32_t>(ss,mol,majorVersion,numAtoms);
      } else {
	_depickle<unsigned char>(ss,mol,majorVersion,numAtoms);
      }
    }
    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();
    if(majorVersion<4000){
      // FIX for issue 220 - probably better to change the pickle format later
      MolOps::assignStereochemistry(*mol,true);
    }
  }
  void MolPickler::molFromPickle(const std::string &pickle,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    std::stringstream ss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
    ss.write(pickle.c_str(),pickle.length());
    MolPickler::molFromPickle(ss,mol);
  }


  //--------------------------------------
  //
  //            Molecules
  //
  //--------------------------------------
  template <typename T>
  void MolPickler::_pickle(const ROMol *mol,std::ostream &ss){
    PRECONDITION(mol,"empty molecule");
    int32_t tmpInt;
    bool includeAtomCoords=true;
    std::map<int,int> atomIdxMap;

    tmpInt = static_cast<int32_t>(mol->getNumAtoms());
    streamWrite(ss,tmpInt);
    tmpInt = static_cast<int32_t>(mol->getNumBonds());
    streamWrite(ss,tmpInt);

    char flag = 0;
    if(includeAtomCoords) flag |= 0x1<<7;
    streamWrite(ss,flag);

    // -------------------
    //
    // Write Atoms
    //
    // -------------------
    streamWrite(ss,BEGINATOM);
    ROMol::ConstAtomIterator atIt;
    int nWritten=0;
    for(atIt=mol->beginAtoms();atIt!=mol->endAtoms();++atIt){
      _pickleAtom<T>(ss,*atIt);
      atomIdxMap[(*atIt)->getIdx()] = nWritten;
      nWritten++;
    }
  
    // -------------------
    //
    // Write Bonds
    //
    // -------------------
    streamWrite(ss,BEGINBOND);
    for(unsigned int i=0;i<mol->getNumBonds();i++){
      _pickleBond<T>(ss,mol->getBondWithIdx(i),atomIdxMap);
    }

    // -------------------
    //
    // Write Rings (if present)
    //
    // -------------------
    const RingInfo *ringInfo=mol->getRingInfo();
    if(ringInfo && ringInfo->isInitialized()){
      streamWrite(ss,BEGINSSSR);
      _pickleSSSR<T>(ss,ringInfo,atomIdxMap);
    }
    
    // pickle the conformations if necessary
    
    if (includeAtomCoords) {
      streamWrite(ss,BEGINCONFS);
      tmpInt = static_cast<int32_t>(mol->getNumConformers());
      streamWrite(ss,tmpInt);
      
      ROMol::ConstConformerIterator ci;
      for (ci = mol->beginConformers(); ci != mol->endConformers(); ++ci) {
        const Conformer *conf = ci->get();
        _pickleConformer<T>(ss, conf);
      }
    }
    streamWrite(ss,ENDMOL);
  }

  template <typename T>
  void MolPickler::_depickle(std::istream &ss,ROMol *mol, int version,int numAtoms){
    PRECONDITION(mol,"empty molecule");
    bool directMap= mol->getNumAtoms()==0;
    Tags tag;
    int32_t tmpInt;
    //int numAtoms,numBonds;
    int numBonds;
    bool haveQuery=false;
    
    streamRead(ss,tmpInt,version);
    numBonds = tmpInt;

    // did we include coordinates
    bool includeCoords=false;
    if (version >= 3000) {
      char flag;
      streamRead(ss,flag,version);
      if (flag & 0x1<<7) includeCoords = true;
    }
    // -------------------
    //
    // Read Atoms
    //
    // -------------------
    streamRead(ss,tag,version);
    if(tag != BEGINATOM){
      throw MolPicklerException("Bad pickle format: BEGINATOM tag not found.");
    }
    Conformer *conf = 0;
    if ((version >= 2000 && version<3000) && includeCoords) {
      // there can only one conformation - since the poositions were stored on
      // the atoms themselves in this version
      conf = new Conformer(numAtoms);
      mol->addConformer(conf, true);
    }
    for(int i=0;i<numAtoms;i++){
      RDGeom::Point3D pos;
      Atom *atom=_addAtomFromPickle<T>(ss,mol,pos,version,directMap);
      if ((version >= 2000 && version<3000) && includeCoords) {
        // this is a older pickle so we go the pos
        conf->setAtomPos(i, pos);
      }
      if(!directMap){
	mol->setAtomBookmark(atom,i);
      }
      if(atom->hasQuery()){
        haveQuery=true;
      }
    }

    // -------------------
    //
    // Read Bonds
    //
    // -------------------
    streamRead(ss,tag,version);
    if(tag != BEGINBOND){
      throw MolPicklerException("Bad pickle format: BEGINBOND tag not found.");
    }
    for(int i=0;i<numBonds;i++){
      Bond *bond=_addBondFromPickle<T>(ss,mol,version,directMap);
      if(!directMap){
        mol->setBondBookmark(bond,i);
      }
    }


    // -------------------
    //
    // Read Rings (if needed)
    //
    // -------------------
    streamRead(ss,tag,version);
    if(tag == BEGINSSSR){
      _addRingInfoFromPickle<T>(ss,mol,version,directMap);
      streamRead(ss,tag,version);
    }
    
    if (tag == BEGINCONFS) {
      // read in the conformation
      streamRead(ss, tmpInt,version);
      int i;
      for (i = 0; i < tmpInt; i++) {
        Conformer *conf = _conformerFromPickle<T>(ss,version);
        mol->addConformer(conf);
      }
      streamRead(ss,tag,version);
    }

    if(tag != ENDMOL){
      throw MolPicklerException("Bad pickle format: ENDMOL tag not found.");
    }

    if(haveQuery){
      // we didn't read any property info for atoms with associated
      // queries. update their property caches
      // (was sf.net Issue 3316407)
      for(ROMol::AtomIterator atIt=mol->beginAtoms();
          atIt!=mol->endAtoms();atIt++){
        Atom *atom = *atIt;
        if(atom->hasQuery()){
          atom->updatePropertyCache(false);
        }
      }
    }
  }


  //--------------------------------------
  //
  //            Atoms
  //
  //--------------------------------------

  namespace {
    bool getAtomMapNumber(const Atom *atom,int &mapNum){
      PRECONDITION(atom,"bad atom");
      if(!atom->hasProp("molAtomMapNumber")) return false;
      bool res=true;
      int tmpInt;
      try{
        atom->getProp("molAtomMapNumber",tmpInt);
      } catch (boost::bad_any_cast &exc) {
        const std::string &tmpSVal=atom->getProp<std::string>("molAtomMapNumber");
        try{
          tmpInt = boost::lexical_cast<int>(tmpSVal);
        } catch(boost::bad_lexical_cast &lexc) {
          res=false;
        }
      }
      if(res) mapNum=tmpInt;
      return res;
    }
  }
  
  // T refers to the type of the atom indices written
  template <typename T>
  void MolPickler::_pickleAtom(std::ostream &ss,const Atom *atom) {
    PRECONDITION(atom,"empty atom");
    char tmpChar;
    signed char tmpSchar;
    int tmpInt;
    char flags;

    tmpChar = atom->getAtomicNum()%128;
    streamWrite(ss,tmpChar);

    flags = 0;
    if(atom->getIsAromatic()) flags |= 0x1<<6;
    if(atom->getNoImplicit()) flags |= 0x1<<5;
    if(atom->hasQuery()) flags |= 0x1<<4;
    if(getAtomMapNumber(atom,tmpInt)) flags |= 0x1<<3;
    if(atom->hasProp("dummyLabel")) flags |= 0x1<<2;
    if(atom->getMonomerInfo()) flags |= 0x1<<1;

    streamWrite(ss,flags);
    
    if(!atom->hasQuery()){
      std::stringstream tss(std::ios_base::binary|std::ios_base::out|std::ios_base::in);
      int32_t propFlags=0;
      // tmpFloat=atom->getMass()-PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
      // if(fabs(tmpFloat)>.0001){
      //   propFlags |= 1;
      //   streamWrite(tss,tmpFloat);
      // }
      tmpSchar=static_cast<signed char>(atom->getFormalCharge());
      if(tmpSchar!=0){
        propFlags |= 1<<1;
        streamWrite(tss,tmpSchar);
      }
      tmpChar = static_cast<char>(atom->getChiralTag());
      if(tmpChar!=0){
        propFlags |= 1<<2;
        streamWrite(tss,tmpChar);
      }
      tmpChar = static_cast<char>(atom->getHybridization());
      if(tmpChar!=static_cast<char>(Atom::SP3)){
        propFlags |= 1<<3;
        streamWrite(tss,tmpChar);
      }

      tmpChar = static_cast<char>(atom->getNumExplicitHs());
      if(tmpChar!=0){
        propFlags |= 1<<4;
        streamWrite(tss,tmpChar);
      }
      if(atom->d_explicitValence>0){
        tmpChar = static_cast<char>(atom->d_explicitValence);
        propFlags |= 1<<5;
        streamWrite(tss,tmpChar);
      }
      if(atom->d_implicitValence>0){
        tmpChar = static_cast<char>(atom->d_implicitValence);
        propFlags |= 1<<6;
        streamWrite(tss,tmpChar);
      }
      tmpChar = static_cast<char>(atom->getNumRadicalElectrons());
      if(tmpChar!=0){
        propFlags |= 1<<7;
        streamWrite(tss,tmpChar);
      }

      unsigned int tmpuint=atom->getIsotope();
      if(tmpuint>0){
        propFlags |= 1<<8;
        streamWrite(tss,tmpuint);
      }

      streamWrite(ss,propFlags);
      ss.write(tss.str().c_str(),tss.str().size());
    } else {
      streamWrite(ss,BEGINQUERY);
      pickleQuery(ss,static_cast<const QueryAtom*>(atom)->getQuery());
      streamWrite(ss,ENDQUERY);
    }
    if(getAtomMapNumber(atom,tmpInt)){
      tmpChar=static_cast<char>(tmpInt%256);
      streamWrite(ss,ATOM_MAPNUMBER,tmpChar);
    }
    if(atom->hasProp("dummyLabel")){
      streamWrite(ss,ATOM_DUMMYLABEL,atom->getProp<std::string>("dummyLabel"));
    }
    if(atom->getMonomerInfo()){
      streamWrite(ss,BEGIN_ATOM_MONOMER);
      pickleAtomMonomerInfo(ss,atom->getMonomerInfo());
      streamWrite(ss,END_ATOM_MONOMER);      
    }
  }

  template <typename T>
  void MolPickler::_pickleConformer(std::ostream &ss,const Conformer *conf) {
    PRECONDITION(conf,"empty conformer");
    char tmpChr = static_cast<int>(conf->is3D());
    streamWrite(ss,tmpChr);
    int32_t tmpInt = static_cast<int32_t>(conf->getId());
    streamWrite(ss,tmpInt);
    T tmpT = static_cast<T>(conf->getNumAtoms());
    streamWrite(ss,tmpT);
    const RDGeom::POINT3D_VECT &pts = conf->getPositions();
    for (RDGeom::POINT3D_VECT_CI pti = pts.begin(); pti != pts.end(); pti++) {
      float tmpFloat;
      tmpFloat = static_cast<float>(pti->x);
      streamWrite(ss,tmpFloat);
      tmpFloat = static_cast<float>(pti->y);
      streamWrite(ss,tmpFloat);
      tmpFloat = static_cast<float>(pti->z);
      streamWrite(ss,tmpFloat);
    } 
  }
    
  template <typename T> 
  Conformer *MolPickler::_conformerFromPickle(std::istream &ss,int version) {
    float tmpFloat;
    bool is3D=true;
    if(version>4000){
      char tmpChr;
      streamRead(ss, tmpChr,version);
      is3D=static_cast<bool>(tmpChr);
    }
    int tmpInt;
    streamRead(ss, tmpInt,version);
    unsigned int cid = static_cast<unsigned int>(tmpInt);
    T tmpT;
    streamRead(ss, tmpT,version);
    unsigned int numAtoms = static_cast<unsigned int>(tmpT);
    Conformer *conf = new Conformer(numAtoms);
    conf->setId(cid);
    conf->set3D(is3D);
    for (unsigned int i = 0; i < numAtoms; i++) {
      streamRead(ss, tmpFloat,version);
      conf->getAtomPos(i).x = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat,version);
      conf->getAtomPos(i).y = static_cast<double>(tmpFloat);
      streamRead(ss, tmpFloat,version);
      conf->getAtomPos(i).z = static_cast<double>(tmpFloat);
    }
    return conf;
  }
    
  template <typename T>
  Atom *MolPickler::_addAtomFromPickle(std::istream &ss,ROMol *mol,
                                       RDGeom::Point3D &pos,
				       int version,bool directMap){
    PRECONDITION(mol,"empty molecule");
    float x,y,z;
    char tmpChar;
    signed char tmpSchar;
    char flags;
    Tags tag;
    Atom *atom=0;
    int atomicNum=0;

    streamRead(ss,tmpChar,version);
    atomicNum=tmpChar;

    bool hasQuery=false;
    streamRead(ss,flags,version);
    if(version>5000){
      hasQuery=flags&0x1<<4;
    }
    if(!hasQuery){
      atom = new Atom(atomicNum);
    } else {
      atom = new QueryAtom();
      if(atomicNum){
        // can't set this in the constructor because that builds a
        // query and we're going to take care of that later:
        atom->setAtomicNum(atomicNum);
      }
    }
    atom->setIsAromatic(flags & 0x1<<6);
    atom->setNoImplicit(flags & 0x1<<5);

    bool hasAtomMap=0,hasDummyLabel=0;
    if(version>=6020){
      hasAtomMap=flags & 0x1<<3;
      hasDummyLabel=flags & 0x1<<2;
    }
    bool hasMonomerInfo=0;
    if(version>=7020){
      hasMonomerInfo=flags & 0x1<<1;
    }

    // are coordinates present?
    if(flags & 0x1<<7){
      streamRead(ss,x,version);
      pos.x = static_cast<double>(x);
      streamRead(ss,y,version);
      pos.y = static_cast<double>(y);
      streamRead(ss,z,version);
      pos.z = static_cast<double>(z);
    }
    
    if(version<=5000 || !hasQuery){
      if(version<7000){
        if(version<6030){
          streamRead(ss,tmpSchar,version);
          atom->setMass(PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum())+
                        static_cast<int>(tmpSchar));
        } else {
          float tmpFloat;
          streamRead(ss,tmpFloat,version);
          atom->setMass(tmpFloat);
        }

        streamRead(ss,tmpSchar,version);
        atom->setFormalCharge(static_cast<int>(tmpSchar));

        streamRead(ss,tmpChar,version);
        atom->setChiralTag(static_cast<Atom::ChiralType>(tmpChar));
        streamRead(ss,tmpChar,version);
        atom->setHybridization(static_cast<Atom::HybridizationType>(tmpChar));
        streamRead(ss,tmpChar,version);
        atom->setNumExplicitHs(static_cast<int>(tmpChar));
        streamRead(ss,tmpChar,version);
        atom->d_explicitValence = tmpChar;
        streamRead(ss,tmpChar,version);
        atom->d_implicitValence = tmpChar;
        if(version>6000){
          streamRead(ss,tmpChar,version);
          atom->d_numRadicalElectrons = static_cast<unsigned int>(tmpChar);
        }
      } else {
        int propFlags;
        streamRead(ss,propFlags,version);
        atom->setMass(PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
        if(propFlags&1){
          float tmpFloat;
          streamRead(ss,tmpFloat,version);
          int iso=static_cast<int>(floor(tmpFloat+atom->getMass()+.0001));
          atom->setIsotope(iso);
        }

        if(propFlags&(1<<1)){
          streamRead(ss,tmpSchar,version);
        } else {
          tmpSchar=0;
        }
        atom->setFormalCharge(static_cast<int>(tmpSchar));

        if(propFlags&(1<<2)){
          streamRead(ss,tmpChar,version);
        } else {
          tmpChar=0;
        }
        atom->setChiralTag(static_cast<Atom::ChiralType>(tmpChar));

        if(propFlags&(1<<3)){
          streamRead(ss,tmpChar,version);
        } else {
          tmpChar=Atom::SP3;
        }
        atom->setHybridization(static_cast<Atom::HybridizationType>(tmpChar));

        if(propFlags&(1<<4)){
          streamRead(ss,tmpChar,version);
        } else {
          tmpChar=0;
        }
        atom->setNumExplicitHs(tmpChar);

        if(propFlags&(1<<5)){
          streamRead(ss,tmpChar,version);
        } else {
          tmpChar=0;
        }
        atom->d_explicitValence=tmpChar;
        
        if(propFlags&(1<<6)){
          streamRead(ss,tmpChar,version);
        } else {
          tmpChar=0;
        }
        atom->d_implicitValence=tmpChar;
        if(propFlags&(1<<7)){
          streamRead(ss,tmpChar,version);
        } else {
          tmpChar=0;
        }          
        atom->d_numRadicalElectrons=static_cast<unsigned int>(tmpChar);

        atom->d_isotope=0;
        if(propFlags&(1<<8)){
          unsigned int tmpuint;
          streamRead(ss,tmpuint,version);
          atom->setIsotope(tmpuint);
        }
      }
      
    } else if(version>5000){
      // we have a query:
      streamRead(ss,tag,version);
      if(tag != BEGINQUERY){
        throw MolPicklerException("Bad pickle format: BEGINQUERY tag not found.");
      }
      static_cast<QueryAtom *>(atom)->setQuery(unpickleQuery(ss,atom,version));
      streamRead(ss,tag,version);
      if(tag != ENDQUERY){
        throw MolPicklerException("Bad pickle format: ENDQUERY tag not found.");
      }
    
      atom->setMass(PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
      atom->setNumExplicitHs(0);

    }

    if(version>5000){
      if(version<6020){
        unsigned int sPos=ss.tellg();
        Tags tag;
        streamRead(ss,tag,version);
        if(tag==ATOM_MAPNUMBER){
          int tmpInt;
          streamRead(ss,tmpChar,version);
          tmpInt=tmpChar;
          atom->setProp("molAtomMapNumber",tmpInt);
        } else {
          ss.seekg(sPos);
        }
      } else {
        if(hasAtomMap) {
          Tags tag;
          streamRead(ss,tag,version);
          if(tag != ATOM_MAPNUMBER){
            throw MolPicklerException("Bad pickle format: ATOM_MAPNUMBER tag not found.");
          }
          int tmpInt;
          streamRead(ss,tmpChar,version);
          tmpInt=tmpChar;
          atom->setProp("molAtomMapNumber",tmpInt);
        }
        if(hasDummyLabel){
          streamRead(ss,tag,version);
          if(tag != ATOM_DUMMYLABEL){
            throw MolPicklerException("Bad pickle format: ATOM_DUMMYLABEL tag not found.");
          }
          std::string tmpStr;
          streamRead(ss,tmpStr,version);
          atom->setProp("dummyLabel",tmpStr);
        }
      }
    }
    if(version>=7020){
      if(hasMonomerInfo){
        streamRead(ss,tag,version);
        if(tag != BEGIN_ATOM_MONOMER){
          throw MolPicklerException("Bad pickle format: BEGIN_ATOM_MONOMER tag not found.");
        }
        atom->setMonomerInfo(unpickleAtomMonomerInfo(ss,version));
      }
    }
    mol->addAtom(atom,false,true);
    return atom;
  }

  //--------------------------------------
  //
  //            Bonds
  //
  //--------------------------------------

  template <typename T>
  void MolPickler::_pickleBond(std::ostream &ss,const Bond *bond,
			       std::map<int,int> &atomIdxMap){
    PRECONDITION(bond,"empty bond");
    T tmpT;
    char tmpChar;
    char flags;

    tmpT = static_cast<T>(atomIdxMap[bond->getBeginAtomIdx()]);
    streamWrite(ss,tmpT);
    tmpT = static_cast<T>(atomIdxMap[bond->getEndAtomIdx()]);
    streamWrite(ss,tmpT);

    flags = 0;
    if(bond->getIsAromatic())	flags |= 0x1<<6;
    if(bond->getIsConjugated()) flags |= 0x1<<5;
    if(bond->hasQuery()) flags |= 0x1<<4;
    if(bond->getBondType()!=Bond::SINGLE) flags |= 0x1<<3;
    if(bond->getBondDir()!=Bond::NONE) flags |= 0x1<<2;
    if(bond->getStereo()!=Bond::STEREONONE) flags |= 0x1<<1;
    streamWrite(ss,flags);
    
    if(bond->getBondType()!=Bond::SINGLE){
      tmpChar = static_cast<char>(bond->getBondType());
      streamWrite(ss,tmpChar);
    }
    if(bond->getBondDir()!=Bond::NONE){
      tmpChar = static_cast<char>(bond->getBondDir());
      streamWrite(ss,tmpChar);
    }

    // write info about the stereochemistry:
    if(bond->getStereo()!=Bond::STEREONONE){
      tmpChar = static_cast<char>(bond->getStereo());
      streamWrite(ss,tmpChar);
      const INT_VECT &stereoAts=bond->getStereoAtoms();
      tmpChar = stereoAts.size();
      streamWrite(ss,tmpChar);
      for(INT_VECT_CI idxIt=stereoAts.begin();idxIt!=stereoAts.end();++idxIt){
        tmpT = static_cast<T>(*idxIt);
        streamWrite(ss,tmpT);
      }
    }
    if(bond->hasQuery()){
      streamWrite(ss,BEGINQUERY);
      pickleQuery(ss,static_cast<const QueryBond*>(bond)->getQuery());
      streamWrite(ss,ENDQUERY);
    }
  }

  template <typename T>
  Bond *MolPickler::_addBondFromPickle(std::istream &ss,ROMol *mol,int version,
				       bool directMap){
    PRECONDITION(mol,"empty molecule");
    char tmpChar;
    char flags;
    int begIdx,endIdx;
    T tmpT;

    Bond *bond=NULL;
    streamRead(ss,tmpT,version);
    if(directMap){
      begIdx=tmpT;
    } else {
      begIdx=mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx();
    }
    streamRead(ss,tmpT,version);
    if(directMap){
      endIdx=tmpT;

    } else {
      endIdx=mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx();
    }
    streamRead(ss,flags,version);
    bool hasQuery=flags&0x1<<4;

    if(version<=5000 || (version<=7000 && !hasQuery) || version>7000){
      bond = new Bond();
      bond->setIsAromatic(flags & 0x1<<6);
      bond->setIsConjugated(flags & 0x1<<5);

      if(version<7000){
        streamRead(ss,tmpChar,version);
        bond->setBondType(static_cast<Bond::BondType>(tmpChar));
        streamRead(ss,tmpChar,version);
        bond->setBondDir(static_cast<Bond::BondDir>(tmpChar));

        if(version>3000){
          streamRead(ss,tmpChar,version);
          Bond::BondStereo stereo=static_cast<Bond::BondStereo>(tmpChar);
          bond->setStereo(stereo);
          if(stereo!=Bond::STEREONONE){
            streamRead(ss,tmpChar,version);
            for(char i=0;i<tmpChar;++i){
              streamRead(ss,tmpT,version);
              bond->getStereoAtoms().push_back(static_cast<int>(tmpT));
            }
          }
        }
      } else {
        if(flags & (0x1<<3)){
          streamRead(ss,tmpChar,version);
          bond->setBondType(static_cast<Bond::BondType>(tmpChar));
        } else {
          bond->setBondType(Bond::SINGLE);
        }
            
        if(flags & (0x1<<2)){
          streamRead(ss,tmpChar,version);
          bond->setBondDir(static_cast<Bond::BondDir>(tmpChar));
        } else {
          bond->setBondDir(Bond::NONE);
        }

        if(flags & (0x1<<1)){
          streamRead(ss,tmpChar,version);
          Bond::BondStereo stereo=static_cast<Bond::BondStereo>(tmpChar);
          bond->setStereo(stereo);
          streamRead(ss,tmpChar,version);
          for(char i=0;i<tmpChar;++i){
            streamRead(ss,tmpT,version);
            bond->getStereoAtoms().push_back(static_cast<int>(tmpT));
          }
        } else {
          bond->setStereo(Bond::STEREONONE);
        }
      }
    }
    if(version>5000 && hasQuery) {
      Tags tag;
      if(bond){
        Bond *tbond=bond;
        bond = new QueryBond(*bond);
        delete tbond;
      } else {
        bond = new QueryBond();
      }

      // we have a query:
      streamRead(ss,tag,version);
      if(tag != BEGINQUERY){
        throw MolPicklerException("Bad pickle format: BEGINQUERY tag not found.");
      }
      static_cast<QueryBond *>(bond)->setQuery(unpickleQuery(ss,bond,version));
      streamRead(ss,tag,version);
      if(tag != ENDQUERY){
        throw MolPicklerException("Bad pickle format: ENDQUERY tag not found.");
      }
    }
    if(bond){
      bond->setBeginAtomIdx(begIdx);
      bond->setEndAtomIdx(endIdx);
      mol->addBond(bond,true);
    }
    return bond;
  }

  //--------------------------------------
  //
  //            Rings
  //
  //--------------------------------------
  template <typename T>
  void MolPickler::_pickleSSSR(std::ostream &ss,const RingInfo *ringInfo,
			       std::map<int,int> &atomIdxMap){
    PRECONDITION(ringInfo,"missing ring info");
    T tmpT;
    tmpT = ringInfo->numRings();
    streamWrite(ss,tmpT);
    for(unsigned int i=0;i<ringInfo->numRings();i++){
      INT_VECT ring;
      ring = ringInfo->atomRings()[i];
      tmpT = static_cast<T>(ring.size());
      streamWrite(ss,tmpT);
      for(unsigned int j=0;j<ring.size();j++){
	tmpT = static_cast<T>(atomIdxMap[ring[j]]);
	streamWrite(ss,tmpT);
      }
#if 0
      ring = ringInfo->bondRings()[i];
      tmpT = static_cast<T>(ring.size());
      streamWrite(ss,tmpT);
      for(unsigned int j=0;j<ring.size();j++){
	tmpT = static_cast<T>(ring[j]);
	streamWrite(ss,tmpT);
      }
#endif
    }
  }

  template <typename T>
  void MolPickler::_addRingInfoFromPickle(std::istream &ss,ROMol *mol,
					  int version,
					  bool directMap){
    PRECONDITION(mol,"empty molecule");
    RingInfo *ringInfo=mol->getRingInfo();
    if(!ringInfo->isInitialized()) ringInfo->initialize();

    T numRings;
    streamRead(ss,numRings,version);

    if(numRings>0){
      ringInfo->preallocate(mol->getNumAtoms(),mol->getNumBonds());
      for(unsigned int i=0;i<static_cast<unsigned int>(numRings);i++){
        T tmpT;
        T ringSize;
        streamRead(ss,ringSize,version);

        INT_VECT atoms(static_cast<int>(ringSize));
        INT_VECT bonds(static_cast<int>(ringSize));
        for(unsigned int j=0; j<static_cast<unsigned int>(ringSize); j++){
          streamRead(ss,tmpT,version);
          if(directMap){
            atoms[j] = static_cast<int>(tmpT);
          } else {
            atoms[j] = mol->getAtomWithBookmark(static_cast<int>(tmpT))->getIdx();
          }
        }
        if(version<7000){
          for(unsigned int j=0; j<static_cast<unsigned int>(ringSize); j++){
            streamRead(ss,tmpT,version);
            if(directMap){
              bonds[j] = static_cast<int>(tmpT);
            } else {
              bonds[j] = mol->getBondWithBookmark(static_cast<int>(tmpT))->getIdx();
            }
          }
        } else {
          for(unsigned int j=1;j<static_cast<unsigned int>(ringSize);++j){
            bonds[j-1]=mol->getBondBetweenAtoms(atoms[j-1],atoms[j])->getIdx();
          }
          bonds[ringSize-1]=mol->getBondBetweenAtoms(atoms[0],atoms[ringSize-1])->getIdx();
        }
        ringInfo->addRing(atoms,bonds);
      }
    }

  }


  //--------------------------------------
  //
  //            Version 1 Pickler:
  //
  //  NOTE: this is not 64bit clean, but it shouldn't be used anymore anyway
  //
  //--------------------------------------

  void MolPickler::_pickleV1(const ROMol *mol,std::ostream &ss){
    PRECONDITION(mol,"empty molecule");
    ROMol::ConstAtomIterator atIt;
    const Conformer *conf = 0;
    if (mol->getNumConformers() > 0) {
      conf = &(mol->getConformer());
    }
    for(atIt=mol->beginAtoms();atIt!=mol->endAtoms();atIt++){
      const Atom *atom = *atIt;

      streamWrite(ss,BEGINATOM);
      streamWrite(ss,ATOM_NUMBER,atom->getAtomicNum());

      streamWrite(ss,ATOM_INDEX,atom->getIdx());

      streamWrite(ss,ATOM_POS);
      RDGeom::Point3D p;
      if (conf) {
        p = conf->getAtomPos(atom->getIdx());
      } 
      streamWrite(ss,p.x);
      streamWrite(ss,p.y);
      streamWrite(ss,p.z);
    

      if(atom->getFormalCharge() != 0){
	streamWrite(ss,ATOM_CHARGE,atom->getFormalCharge());
      }
      if(atom->getNumExplicitHs() != 0){
	streamWrite(ss,ATOM_NEXPLICIT,atom->getNumExplicitHs());
      }
      if(atom->getChiralTag() != 0){
	streamWrite(ss,ATOM_CHIRALTAG,atom->getChiralTag());
      }
      if(atom->getMass() != 0.0){
	streamWrite(ss,ATOM_MASS,atom->getMass());
      }
      if(atom->getIsAromatic()){
	streamWrite(ss,ATOM_ISAROMATIC,static_cast<char>(atom->getIsAromatic()));
      }
      streamWrite(ss,ENDATOM);
    }
  
    ROMol::ConstBondIterator bondIt;
    for(bondIt=mol->beginBonds();bondIt!=mol->endBonds();bondIt++){
      const Bond *bond = *bondIt;
      streamWrite(ss,BEGINBOND);
      streamWrite(ss,BOND_INDEX,bond->getIdx());
      streamWrite(ss,BOND_BEGATOMIDX,bond->getBeginAtomIdx());
      streamWrite(ss,BOND_ENDATOMIDX,bond->getEndAtomIdx());
      streamWrite(ss,BOND_TYPE,bond->getBondType());
      if(bond->getBondDir()){
	streamWrite(ss,BOND_DIR,bond->getBondDir());
      }
      streamWrite(ss,ENDBOND);
    }
    streamWrite(ss,ENDMOL);
  }

  void MolPickler::_depickleV1(std::istream &ss,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    Tags tag;

    Conformer *conf = new Conformer();
    mol->addConformer(conf);
    streamRead(ss,tag,1);
    while(tag != ENDMOL){
      switch(tag){
      case BEGINATOM:
	_addAtomFromPickleV1(ss,mol);
	break;
      case BEGINBOND:
	_addBondFromPickleV1(ss,mol);
	break;
      default:
	UNDER_CONSTRUCTION("bad tag in pickle");
      }
      streamRead(ss,tag,1);
    }
    mol->clearAllAtomBookmarks();
    mol->clearAllBondBookmarks();
  }
  
  void MolPickler::_addAtomFromPickleV1(std::istream &ss,ROMol *mol) {
    PRECONDITION(mol,"empty molecule");
    Tags tag;
    int intVar;
    double dblVar;
    char charVar;
    int version=1;
    streamRead(ss,tag,version);
    Atom *atom = new Atom();
    Conformer &conf = mol->getConformer();
    RDGeom::Point3D pos;
    while(tag != ENDATOM){
      switch(tag){
      case ATOM_INDEX:
	streamRead(ss,intVar,version);
	mol->setAtomBookmark(atom,intVar);
	break;
      case ATOM_NUMBER:
	streamRead(ss,intVar,version);
	atom->setAtomicNum(intVar);
	break;
      case ATOM_POS:
	streamRead(ss,pos.x,version);
	streamRead(ss,pos.y,version);
	streamRead(ss,pos.z,version);
        break;
      case ATOM_CHARGE:
	streamRead(ss,intVar,version);
	atom->setFormalCharge(intVar);
	break;
      case ATOM_NEXPLICIT:
	streamRead(ss,intVar,version);
	atom->setNumExplicitHs(intVar);
	break;
      case ATOM_CHIRALTAG:
	streamRead(ss,intVar,version);
	atom->setChiralTag(static_cast<Atom::ChiralType>(intVar));
	break;
      case ATOM_MASS:
	streamRead(ss,dblVar,version);
	atom->setMass(dblVar);
	break;
      case ATOM_ISAROMATIC:
	streamRead(ss,charVar,version);
	atom->setIsAromatic(charVar);
	break;
      default:
	ASSERT_INVARIANT(0,"bad tag in atom block of pickle");
      }
      streamRead(ss,tag,version);
    }
    unsigned int id = mol->addAtom(atom,false,true);
    conf.setAtomPos(id, pos);
  }
  void MolPickler::_addBondFromPickleV1(std::istream &ss,ROMol *mol){
    PRECONDITION(mol,"empty molecule");
    Tags tag;
    int intVar,idx=-1;
    int version=1;
    Bond::BondType bt;
    Bond::BondDir bd;
    streamRead(ss,tag,version);
    Bond *bond = new Bond();
    while(tag != ENDBOND){
      switch(tag){
      case BOND_INDEX:
	streamRead(ss,idx,version);
	break;
      case BOND_BEGATOMIDX:
	streamRead(ss,intVar,version);
	bond->setBeginAtomIdx(mol->getAtomWithBookmark(intVar)->getIdx());
	break;
      case BOND_ENDATOMIDX:
	streamRead(ss,intVar,version);
	bond->setEndAtomIdx(mol->getAtomWithBookmark(intVar)->getIdx());
	break;
      case BOND_TYPE:
	streamRead(ss,bt,version);
	bond->setBondType(bt);
	break;
      case BOND_DIR:
	streamRead(ss,bd,version);
	bond->setBondDir(bd);
	break;
      default:
	ASSERT_INVARIANT(0,"bad tag in bond block of pickle");
      }
      streamRead(ss,tag,version);
    }
    mol->addBond(bond,true);
  }
};
