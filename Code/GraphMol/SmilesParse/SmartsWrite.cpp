// $Id$
//
//  Copyright (C) 2002-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SmartsWrite.h"
#include <sstream>
#include "SmilesWrite.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Canon.h>
#include <RDGeneral/RDLog.h>

namespace RDKit {
  using namespace Canon;

  // local utility namespace
  namespace {

    std::string _recurseGetSmarts(const QueryAtom::QUERYATOM_QUERY *node, bool negate);
    std::string _recurseBondSmarts(const QueryBond::QUERYBOND_QUERY *node, bool negate);

    bool _checkForOrAndLowAnd(std::string smarts) {
      int orLoc, andLoc;
      // if we're a pure recursive smarts, we don't need to worry about this
      if(smarts[0]=='$' && smarts[smarts.size()-1]==')') return false;
      orLoc = smarts.find(",");
      andLoc = smarts.find(";");
      if ( (orLoc > 0) && (andLoc > 0)) {
        return true;
      }
      else {
        return false;
      }
    }

    std::string _combineChildSmarts(std::string cs1, std::string cs2, 
                                    std::string descrip) {
      std::string res = "";
      if ((descrip.find("Or") > 0) && (descrip.find("Or") < descrip.length()) ) {
        // if either of child smarts already have a "," and ";" we can't have one
        // more OR here
        if ( _checkForOrAndLowAnd(cs1) || _checkForOrAndLowAnd(cs2) ) {
          throw "This is a non-smartable query - OR above and below AND in the binary tree";
        }
        res += cs1;
        res += ",";
        res += cs2;
      }
      else if ( (descrip.find("And") > 0)  && (descrip.find("And") < descrip.length()) ){
        int orLoc1, orLoc2;
        std::string symb;
        orLoc1 = cs1.find(',');
        orLoc2 = cs2.find(',');
        if ( (orLoc1 > 0) || (orLoc2 > 0) ) {
          symb = ";";
        }
        else {
          symb = "&";
        }
        res += cs1;
        res += symb;
        res += cs2;
      }
      else {
        std::stringstream err;
        err << "Don't know how to combine using " << descrip;
        throw err.str();
      }
      return res;
    }

    template <typename T>
    void describeQuery(const T *query,std::string leader="\t"){
      //BOOST_LOG(rdInfoLog) << leader << query->getDescription() << std::endl;
      typename T::CHILD_VECT_CI iter;
      for(iter=query->beginChildren();iter!=query->endChildren();++iter){
        describeQuery(iter->get(),leader+"\t");
      }
    }
    
    // called with the children of AND queries:
    std::string smartsOrganicAtom(const QueryAtom::QUERYATOM_QUERY *child1, 
                                  const QueryAtom::QUERYATOM_QUERY *child2) {
      PRECONDITION(child1&&child2,"bad query");

      std::string desc1 = child1->getDescription();
      std::string desc2 = child2->getDescription();
      const QueryAtom::QUERYATOM_QUERY *origA, *otherA;
      if (desc1 == "AtomAtomicNum") {
        origA = child1;
        otherA = child2;
      }
      else {
        origA = child2;
        otherA = child1;
      }
    
      std::string odsc = otherA->getDescription();
      CHECK_INVARIANT( ((odsc == "AtomIsAliphatic") || (odsc == "AtomIsAromatic")),
                       "Should be either an aromatic or Aliphatic atom");

      const ATOM_EQUALS_QUERY *torig = static_cast<const ATOM_EQUALS_QUERY *>(origA);
      int val = torig->getVal();
      std::string res = PeriodicTable::getTable()->getElementSymbol(val);
      if (odsc == "AtomIsAromatic") {
        // if aromatic convert the first leter to a small letter
        res[0] += ('a' - 'A');
        if(res.length()>1){
          res = "["+res+"]";
        }
      }
      if (torig->getNegation()) {
        res = "!" + res;
      }
      return res;
    }

    std::string getAtomSmartsSimple(const ATOM_EQUALS_QUERY *query, bool &needParen) {
      PRECONDITION(query,"bad query");
    
      std::string descrip = query->getDescription();
      std::stringstream res;
      if (descrip == "AtomImplicitHCount") {
        res << "h" << query->getVal();
        needParen = true;
      }else if (descrip == "AtomHasImplicitH") {
        res << "h";
        needParen = true;
      } else if (descrip == "AtomTotalValence") {
        res << "v" << query->getVal();
        needParen = true;
      } else if (descrip == "AtomAtomicNum") {
        res << "#" << query->getVal();
        needParen = true;
      } else if (descrip == "AtomExplicitDegree") {
        res << "D" << query->getVal();
        needParen = true;
      } else if (descrip == "AtomTotalDegree") {
        res << "X" << query->getVal();
        needParen = true;
      } else if (descrip == "AtomHasRingBond") {
        res << "x";
        needParen = true;
      } else if (descrip == "AtomHCount") {
        res << "H" << query->getVal();
        needParen = true;
      } else if (descrip == "AtomIsAliphatic") {
        res << "A";
        needParen = false;
      } else if (descrip == "AtomIsAromatic") {
        res << "a";
        needParen = false;
      } else if (descrip == "AtomNull") {
        res << "*";
        needParen = false;
      } else if (descrip == "AtomInRing") {
        res << "R";
        needParen = true;
      } else if (descrip == "AtomMinRingSize") {
        res << "r" << query->getVal();
        needParen = true;
      } else if (descrip == "AtomInNRings") {
        res << "R";
        if(query->getVal()>=0){
          res << query->getVal();
        }
        needParen = true;
      } else if (descrip == "AtomFormalCharge") {
        int val = query->getVal();
        if (val < 0) {
          res << "-";
        }
        else {
          res << "+";
        }
        if (abs(val) != 1) {
          res << abs(val);
        }
        needParen = true;
      }
      else if (descrip == "AtomHybridization") {
        res << "^";
        switch(query->getVal()){
        case Atom::S:   res << "0"; break;
        case Atom::SP:  res << "1"; break;
        case Atom::SP2: res << "2"; break;
        case Atom::SP3: res << "3"; break;
        }
      } else if (descrip == "AtomMass") {
        res << query->getVal()/massIntegerConversionFactor <<"*";
        needParen = true;
      } else if (descrip == "AtomIsotope") {
        res << query->getVal()<<"*";
        needParen = true;
      } else if (descrip == "AtomHasRingBond") {
        res << "x";
        needParen = true;
      } else if (descrip == "AtomRingBondCount") {
        res << "x" << query->getVal();
        needParen = true;
      } else if (descrip == "AtomUnsaturated") {
        res<<"$(*=,:,#*)";
        needParen = true;
      } else {
        BOOST_LOG(rdWarningLog)<<"Cannot write SMARTS for query type : " << descrip << ". Ignoring it." <<std::endl;
        res<<"*";
      }
      return res.str();
    }

    std::string getRecursiveStructureQuerySmarts(const QueryAtom::QUERYATOM_QUERY *query){
      PRECONDITION(query,"bad query");
      PRECONDITION(query->getDescription()=="RecursiveStructure","bad query");
      const RecursiveStructureQuery *rquery = static_cast<const RecursiveStructureQuery *>(query);
      ROMol *qmol = const_cast<ROMol *>(rquery->getQueryMol());
      std::string res = MolToSmarts(*qmol);
      res = "$(" + res + ")";
      if (rquery->getNegation()) {
        res = "!" + res;
      }
      return res;
    }
    
    std::string getBondSmartsSimple(const BOND_EQUALS_QUERY *bquery) {
      PRECONDITION(bquery,"bad query");
      std::string descrip = bquery->getDescription();
      std::string res = "";
      if (descrip == "BondNull") {
        res += "~";
      }
      else if (descrip == "BondInRing") {
        res += "@";
      }
      else if (descrip == "BondDir") {
        int val = bquery->getVal();
        if (val == static_cast<int>(Bond::ENDDOWNRIGHT) ) {
          res += "\\";
        }
        else if (val == static_cast<int>(Bond::ENDUPRIGHT) ) {
          res += "/";
        }
        else {
          throw "Can't write smarts for this bond dir type";
        }
      }
      else if (descrip == "BondOrder") {
        int val = bquery->getVal();
        if (val == static_cast<int>(Bond::SINGLE) ) {
          res += "-";
        }
        else if (val == static_cast<int>(Bond::DOUBLE) ) {
          res += "=";
        }
        else if (val == static_cast<int>(Bond::TRIPLE)) {
          res += "#";
        }
        else if (val == static_cast<int>(Bond::AROMATIC) ) {
          res += ":";
        }
      }
      else {
        std::stringstream msg;
        msg << "Canot write smarts for this query bond type : " << descrip;
        throw msg.str().c_str();
      }
      return res;
    }
    
    std::string _recurseGetSmarts(const QueryAtom::QUERYATOM_QUERY *node, bool negate) {
      PRECONDITION(node,"bad node");
      // the algorithm goes something like this
      // - recursively get the smarts for the child queries
      // - combine the child smarts using the following rules:
      //      - if we are currently at an OR query, combine the subqueries with a ",",
      //        but only if neither of child smarts do not contain "," and ";" 
      //        This situation leads to a no smartable situation and throw an error
      //      - if we are currently at an and query, combine the child smarts with "&" 
      //        if neither of the child smarts contain a "," - otherwise combine them
      //        the child smarts with a ";"
      // 
      // There is an additional complication with composite nodes that carry a negation - in this
      // case we will propogate the neagtion to the child nodes using the following rules
      //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
      //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))

      std::string descrip = node->getDescription();
      std::string res="";

      const QueryAtom::QUERYATOM_QUERY *child1;
      const QueryAtom::QUERYATOM_QUERY *child2;
      QueryAtom::QUERYATOM_QUERY::CHILD_VECT_CI chi;
      chi = node->beginChildren();
      child1 = chi->get();
      chi++;
      child2 = chi->get();
      chi++;
      // OK we should be at the end of vector by now - since we can have only two children, 
      // well - at least in this case 
      CHECK_INVARIANT(chi == node->endChildren(), "Too many children on the query");
    
      std::string dsc1, dsc2;
      dsc1 = child1->getDescription();
      dsc2 = child2->getDescription();
      std::string csmarts1, csmarts2;

      bool needParen;

      // deal with any special AND cases
      //  1. This "node" is an AtomAnd between a AliphaticAtom (or AromaticAtom) and 
      //      an organic atom e.g. "C"
      if (descrip == "AtomAnd") {
        bool specialCase = false;
        //case 1
        if ( ((dsc1 == "AtomAtomicNum") &&
              ((dsc2 == "AtomIsAliphatic") ||
               (dsc2 == "AtomIsAromatic")))
             ||
             ((dsc2 == "AtomAtomicNum") &&
              ((dsc1 == "AtomIsAliphatic") ||
               (dsc1 == "AtomIsAromatic"))) ) {
          // we trap this one because it's nicer to see
          //   "CC" in the output than "[#6&A][#6&A]"
          res = smartsOrganicAtom(child1, child2);
          specialCase = true;
        }
        if (specialCase) {
          if (negate) {
            res = "!" + res;
          }
          return res;
        }
      }
      
      // deal with the first child
      if (dsc1=="RecursiveStructure") {
        csmarts1 = getRecursiveStructureQuerySmarts(child1);
      } else if ((dsc1 != "AtomOr") && (dsc1 != "AtomAnd")) {
        // child 1 is a simple node
        const ATOM_EQUALS_QUERY *tchild = static_cast<const ATOM_EQUALS_QUERY *>(child1);
        csmarts1 = getAtomSmartsSimple(tchild, needParen);
        bool nneg = (negate)^(tchild->getNegation());
        if (nneg) {
          csmarts1 = "!" + csmarts1;
        }
      } else {
        // child 1 is composite node - recurse
        bool nneg=(negate)^(child1->getNegation());
        csmarts1 = _recurseGetSmarts(child1, nneg);
      }
    
      // deal with the second child
      if (dsc2=="RecursiveStructure") {
        csmarts2 = getRecursiveStructureQuerySmarts(child2);
      } else if ((dsc2 != "AtomOr") && (dsc2 != "AtomAnd")) {
        // child 2 is a simple node
        const ATOM_EQUALS_QUERY *tchild = static_cast<const ATOM_EQUALS_QUERY *>(child2);
        csmarts2 = getAtomSmartsSimple(tchild, needParen);
        bool nneg = (negate)^(tchild->getNegation());
        if (nneg) {
          csmarts2 = "!" + csmarts2;
        }
      }
      else {
        bool nneg=(negate)^(child2->getNegation());
        csmarts2 = _recurseGetSmarts(child2, nneg);
      }
    
      // ok if we have a negation and we have an OR , we have to change to
      // an AND since we propogated the negation
      // i.e NOT (A OR B) = (NOT (A)) AND (NOT(B))
      if (negate) {
        if (descrip == "AtomOr") {
          descrip = "AtomAnd";
        }
        else if (descrip == "AtomAnd") {
          descrip = "AtomOr";
        }
      }

      res += _combineChildSmarts(csmarts1, csmarts2, descrip);
      return res;
    }

    std::string _recurseBondSmarts(const QueryBond::QUERYBOND_QUERY *node, bool negate) {
      // the algorithm goes something like this
      // - recursively get the smarts for the child query bonds
      // - combine the child smarts using the following rules:
      //      - if we are currently at an OR query, combine the subqueries with a ",",
      //        but only if neither of child smarts do not contain "," and ";" 
      //        This situation leads to a no smartable situation and throw an error
      //      - if we are currently at an and query, combine the child smarts with "&" 
      //        if neither of the child smarts contain a "," - otherwise combine them
      //        the child smarts with a ";"
      // 
      // There is an additional complication with composite nodes that carry a negation - in this
      // case we will propogate the neagtion to the child nodes using the following rules
      //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
      //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))
      PRECONDITION(node,"bad node");
      std::string descrip = node->getDescription();
      std::string res="";

      const QueryBond::QUERYBOND_QUERY *child1;
      const QueryBond::QUERYBOND_QUERY *child2;
      QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI chi;

      chi = node->beginChildren();
      child1 = chi->get();
      chi++;
      child2 = chi->get();
      chi++;
      // OK we should be at the end of vector by now - since we can have only two children, 
      // well - atleat in this case 
      CHECK_INVARIANT(chi == node->endChildren(), "Too many children on the query");

      std::string dsc1, dsc2;
      dsc1 = child1->getDescription();
      dsc2 = child2->getDescription();
      std::string csmarts1, csmarts2;

      if ((dsc1 != "BondOr") && (dsc1 != "BondAnd") ) {
        // child1 is  simple node get the smarts directly
        const BOND_EQUALS_QUERY *tchild  = static_cast<const BOND_EQUALS_QUERY *>(child1);
        csmarts1 = getBondSmartsSimple(tchild);
        bool nneg = (negate)^(tchild->getNegation());
        if (nneg) {
          csmarts1 = "!" + csmarts1;
        }
      }
      else {
        // child1 is a composite node recurse further
        bool nneg=(negate)^(child1->getNegation());
        csmarts1 = _recurseBondSmarts(child1, nneg);
      }
    
      // now deal with the second child
      if ((dsc2 != "BondOr") && (dsc2 != "BondAnd")) {
        // child 2 is a simple node
        const BOND_EQUALS_QUERY *tchild  = static_cast<const BOND_EQUALS_QUERY *>(child2);
        csmarts2 = getBondSmartsSimple(tchild);
        bool nneg = (negate)^(tchild->getNegation());
        if (nneg) {
          csmarts2 = "!" + csmarts2;
        }
      }
      else {
        // child two is a composite node - recurse
        bool nneg=(negate)^(child2->getNegation());
        csmarts1 = _recurseBondSmarts(child2, nneg);
      }
    
      // ok if we have a negation and we have to change the underlying logic, since we propogated the negation
      // i.e NOT (A OR B) = (NOT (A)) AND (NOT(B))
      if (negate) {
        if (descrip == "BondOr") {
          descrip = "BondAnd";
        }
        else if (descrip == "BondAnd") {
          descrip = "BondOr";
        }
      }
      res += _combineChildSmarts(csmarts1, csmarts2, descrip);
      return res;
    }
    
    std::string FragmentSmartsConstruct(ROMol &mol,unsigned int atomIdx,
                                        std::vector<Canon::AtomColors> &colors,
                                        INT_VECT &ranks){
      Canon::MolStack molStack;
      molStack.reserve(mol.getNumAtoms() + mol.getNumBonds());
      std::stringstream res;
    
      // this is dirty trick get around the fact that canonicalizeFragment
      // thinks we already called findSSSR - to do some atom ranking
      // but for smarts we are going to ignore that part. We will artificially
      // set the "SSSR" property to an empty property
      VECT_INT_VECT rings;
      mol.getRingInfo()->reset();
      mol.getRingInfo()->initialize();
      Canon::canonicalizeFragment(mol,atomIdx,colors,ranks,
                                  molStack);

      // now clear the "SSSR" property
      mol.getRingInfo()->reset();
      Canon::MolStack::const_iterator msCI;
      for (msCI = molStack.begin(); msCI != molStack.end(); msCI++) {
        switch(msCI->type) {
        case Canon::MOL_STACK_ATOM:
          {
            QueryAtom *qatm = static_cast<QueryAtom *>(msCI->obj.atom);
            res << SmartsWrite::GetAtomSmarts(qatm);
            break;
          }
        case Canon::MOL_STACK_BOND:
          {
            QueryBond *qbnd = static_cast<QueryBond *>(msCI->obj.bond);
            res << SmartsWrite::GetBondSmarts(qbnd);
            break;
          }
        case Canon::MOL_STACK_RING :
          {
            if (msCI->number < 10) 
              res << msCI->number;
            else
              res << "%" << msCI->number;
            break;
          }
        case Canon::MOL_STACK_BRANCH_OPEN:
          {
            res << "(";
            break;
          }
        case Canon::MOL_STACK_BRANCH_CLOSE:
          {
            res << ")";
            break;
          }
        default:
          break;
        }
      }
      return res.str();
    }

    // this is the used when converting a SMILES or
    // non-query atom from a mol file into SMARTS.
    std::string getNonQueryAtomSmarts(const QueryAtom *qatom){
      PRECONDITION(qatom,"bad atom");
      PRECONDITION(!qatom->hasQuery(),"atom should not have query");
      std::stringstream res;
      res<<"[";
      if(SmilesWrite::inOrganicSubset(qatom->getAtomicNum())){
        res<<"#"<<qatom->getAtomicNum();
      } else {
        res<<qatom->getSymbol();
      }
      int hs=qatom->getNumExplicitHs();
      // FIX: probably should be smarter about Hs:
      if(hs){
        res<<"H";
        if(hs>1) res<<hs;
      }
      int chg=qatom->getFormalCharge();
      if(chg){
        if(chg==-1){
          res<<"-";
        } else if(chg==1){
          res<<"+";
        } else if(chg<0){
          res<<qatom->getFormalCharge();
        } else {
          res<<"+"<<qatom->getFormalCharge();
        }
      }
      res<<"]";
      return res.str();
    }
    
    // this is the used when converting a SMILES or
    // non-query bond from a mol file into SMARTS.
    std::string getNonQueryBondSmarts(const QueryBond *qbond){
      PRECONDITION(qbond,"bad bond");
      PRECONDITION(!qbond->hasQuery(),"bond should not have query");
      std::string res;

      if(qbond->getIsAromatic()){
        res=":";
      } else {
        switch(qbond->getBondType()){
        case Bond::SINGLE:
          res="-";
          break;
        case Bond::DOUBLE:
          res="=";
          break;
        case Bond::TRIPLE:
          res="#";
          break;
        case Bond::AROMATIC: 
          res=":"; 
          break;
        default:
          // do nothing (i.e. match anything)
          res="";
          break;
        }
      }
      return res;
    }

  } // end of local utility namespace

  namespace SmartsWrite {      
    std::string GetAtomSmarts(const QueryAtom *qatom) {
      PRECONDITION(qatom,"bad atom");
      std::string res;
      bool needParen=false;
    
      //BOOST_LOG(rdInfoLog)<<"Atom: " <<qatom->getIdx()<<std::endl;
      if(!qatom->hasQuery()){
        res =getNonQueryAtomSmarts(qatom);
        //BOOST_LOG(rdInfoLog)<<"\tno query:" <<res;
        return res;
      }  
      QueryAtom::QUERYATOM_QUERY *query = qatom->getQuery();
      //describeQuery(query);

      PRECONDITION(query,"atom has no query");
      std::string descrip = qatom->getQuery()->getDescription();
      if (descrip == ""){
        // we have simple atom - just generate the smiles and return
        res = SmilesWrite::GetAtomSmiles(qatom);
	if(res[0]=='['){
	  // chop the brackets off, we'll put them back on later:
	  needParen=true;
	  res = res.substr(1,res.size()-2);
	}
      }
      else if ((descrip == "AtomOr") || (descrip == "AtomAnd")) {
        // we have a composite query 
        needParen = true;
        res = _recurseGetSmarts(query, query->getNegation());
        if (res.length() == 1) { // single atom symbol we don't need parens
          needParen = false;
        }
      } else if(descrip == "RecursiveStructure"){
        // it's a bare recursive structure query:
        res = getRecursiveStructureQuerySmarts(query);
        needParen = true;
      }
      else {// we have a simple smarts 
        ATOM_EQUALS_QUERY *tquery = static_cast<ATOM_EQUALS_QUERY *>(qatom->getQuery());
        res = getAtomSmartsSimple(tquery, needParen);
        if (tquery->getNegation()) {
          res = "!" + res;
        }
      }
      if(qatom->hasProp("molAtomMapNumber")){
	needParen=true;
	std::string mapNum;
	qatom->getProp("molAtomMapNumber",mapNum);
	res += ":"+mapNum;
      }
      if (needParen ) {
        res = "[" + res + "]";
      }
      return res;
    }

    std::string GetBondSmarts(const QueryBond *bond) {
      PRECONDITION(bond,"bad bond");
      std::string res = "";
    
      //BOOST_LOG(rdInfoLog)<<"bond: " <<bond->getIdx()<<std::endl;;
      // it is possible that we are regular single bond and we don't need to write anything
      if(!bond->hasQuery()){
        res=getNonQueryBondSmarts(bond);
        //BOOST_LOG(rdInfoLog)<<"\tno query:" <<res;
        return res;
      }
      //describeQuery(bond->getQuery());
      if ((typeid(*bond) == typeid(Bond)) && 
          ( (bond->getBondType() == Bond::SINGLE) || (bond->getBondType() == Bond::AROMATIC)) ) {
        return res;
      }
      const QueryBond::QUERYBOND_QUERY *query = bond->getQuery();
      PRECONDITION(query,"bond has no query");
      std::string descrip = query->getDescription();

      if ((descrip == "BondAnd") || (descrip == "BondOr") ) {
        // composite query
        res = _recurseBondSmarts(query, query->getNegation());
      }
      else {
        //simple query
        if (query->getNegation()) {
          res = "!";
        }
        const BOND_EQUALS_QUERY *tquery  = static_cast<const BOND_EQUALS_QUERY *>(query);
        res += getBondSmartsSimple(tquery);
      }
      return res;
    }
  } // end of namespace SmartsWrite
  
  std::string MolToSmarts(ROMol &inmol, bool doIsomericSmiles) {
    std::string res;
    unsigned int nAtoms = inmol.getNumAtoms();
    if(!nAtoms) return "";

    ROMol mol(inmol);
    INT_VECT ranks;
    ranks.resize(nAtoms);
    // For smiles writing we would be canonicalizing but we will not do that here.
    // We will simple use the atom indices as the rank
    for(ROMol::AtomIterator atIt=mol.beginAtoms();atIt!=mol.endAtoms();atIt++) {
      ranks.push_back((*atIt)->getIdx());
    }
    
    std::vector<AtomColors> colors;
    colors.resize(nAtoms);
    std::vector<AtomColors>::iterator colorIt;
    for(colorIt=colors.begin();colorIt!=colors.end();colorIt++)
      *colorIt=Canon::WHITE_NODE;
    
    colorIt = std::find(colors.begin(),colors.end(),Canon::WHITE_NODE);
    while (colorIt != colors.end()) {
      unsigned int nextAtomIdx=0;
      int nextRank;
      std::string subSmi;
      nextRank = nAtoms + 1;
      for (unsigned int i = 0; i < nAtoms; i++) {
        if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
          nextRank = ranks[i];
          nextAtomIdx = i;
        }
      }
      
      subSmi = FragmentSmartsConstruct(mol, nextAtomIdx, colors, 
                                       ranks);
      res += subSmi;
      
      colorIt = std::find(colors.begin(),colors.end(),Canon::WHITE_NODE);
      if (colorIt != colors.end()) {
        res += ".";
      }
    }
    return res;
  }

}
