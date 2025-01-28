%{

  //
  //  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
  //
  //   @@ All Rights Reserved  @@
  //
#include <cstring>
#include <iostream>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <RDGeneral/RDLog.h>

#define YYDEBUG 1
#include "smarts.tab.hpp"

extern int yysmarts_lex(YYSTYPE *,void *, int &);

using namespace RDKit;
namespace {
 void yyErrorCleanup(std::vector<RDKit::RWMol *> *molList){
  for(std::vector<RDKit::RWMol *>::iterator iter=molList->begin();
      iter != molList->end(); ++iter){
     SmilesParseOps::CleanupAfterParseError(*iter);
     delete *iter;
  }
  molList->clear();
  molList->resize(0);
 }
}
void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                RDKit::Atom* &,
                RDKit::Bond* &,
                unsigned int &,unsigned int &,
                std::list<unsigned int> *,
		void *,int , const char *msg  )
{
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error: " << msg << " while parsing: " << input << std::endl;
}

void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::list<unsigned int> *,
		void *,int, const char * msg )
{
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error: " << msg << " while parsing: " << input << std::endl;
}


%}

%define api.pure full
%lex-param   {yyscan_t *scanner}
%lex-param   {int& start_token}
%parse-param {const char *input}
%parse-param {std::vector<RDKit::RWMol *> *molList}
%parse-param {RDKit::Atom* &lastAtom}
%parse-param {RDKit::Bond* &lastBond}
%parse-param {unsigned &numAtomsParsed}
%parse-param {unsigned &numBondsParsed}
%parse-param {std::list<unsigned int> *branchPoints}
%parse-param {void *scanner}
%parse-param {int& start_token}

%code provides {
#ifndef YY_DECL
#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)
#endif
}

%union {
  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  RDKit::Atom::ChiralType chiraltype;
  int                      ival;
}

%token START_MOL START_ATOM START_BOND;
%token <ival> AROMATIC_ATOM_TOKEN ORGANIC_ATOM_TOKEN
%token <atom> ATOM_TOKEN
%token <atom> SIMPLE_ATOM_QUERY_TOKEN COMPLEX_ATOM_QUERY_TOKEN
%token <atom> RINGSIZE_ATOM_QUERY_TOKEN RINGBOND_ATOM_QUERY_TOKEN IMPLICIT_H_ATOM_QUERY_TOKEN
%token <atom> HYB_TOKEN HETERONEIGHBOR_ATOM_QUERY_TOKEN ALIPHATIC ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN
%token <ival> ZERO_TOKEN NONZERO_DIGIT_TOKEN
%token GROUP_OPEN_TOKEN GROUP_CLOSE_TOKEN SEPARATOR_TOKEN
%token RANGE_OPEN_TOKEN RANGE_CLOSE_TOKEN
%token HASH_TOKEN MINUS_TOKEN PLUS_TOKEN
%token H_TOKEN AT_TOKEN PERCENT_TOKEN
%token ATOM_OPEN_TOKEN ATOM_CLOSE_TOKEN
%token NOT_TOKEN AND_TOKEN OR_TOKEN SEMI_TOKEN BEGIN_RECURSE END_RECURSE
%token COLON_TOKEN UNDERSCORE_TOKEN
%token <bond> BOND_TOKEN
%token <chiraltype> CHI_CLASS_TOKEN 
%type <moli>  mol
%type <atom> atomd simple_atom hydrogen_atom
%type <atom> atom_expr point_query atom_query recursive_query possible_range_query
%type <ival> ring_number nonzero_number number charge_spec digit
%type <bond> bondd bond_expr bond_query
%token EOS_TOKEN

%left SEMI_TOKEN
%left OR_TOKEN
%left AND_TOKEN
%right NOT_TOKEN

%destructor { delete $$; } <atom>
%destructor { delete $$; } <bond>

%start meta_start

%%

/* --------------------------------------------------------------- */
meta_start:
START_MOL mol {
// the molList has already been updated, no need to do anything
}
| START_ATOM atomd EOS_TOKEN {
  lastAtom = $2;
  YYACCEPT;
}
| START_ATOM bad_atom_def {
  YYABORT;
}
| START_ATOM {
  YYABORT;
}
| START_BOND bond_expr EOS_TOKEN {
  lastBond = $2;
  YYACCEPT;
}
| START_BOND bond_expr {
  delete $2;
  YYABORT;
}
| START_BOND {
  YYABORT;
}
| meta_start error EOS_TOKEN{
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
| meta_start EOS_TOKEN {
  YYACCEPT;
}
| error EOS_TOKEN {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
;

bad_atom_def:
ATOM_OPEN_TOKEN bad_atom_def
| ATOM_CLOSE_TOKEN bad_atom_def
| COLON_TOKEN bad_atom_def
| atom_expr {
  delete $1;
  YYABORT;
}
;

/* --------------------------------------------------------------- */
// FIX: mol MINUS DIGIT
mol: atomd {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  $1->setProp(RDKit::common_properties::_SmilesStart,1);
  (*molList)[ sz ]->addAtom($1,true,true);
  //delete $1;
  $$ = sz;
}
| mol atomd       {
  RWMol *mp = (*molList)[$$];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom($2,true,true);

  QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(a1,mp->getAtomWithIdx(atomIdx2));
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond(newB,true);
}

| mol bond_expr atomd  {
  RWMol *mp = (*molList)[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($3,true,true);
  if( $2->getBondType() == Bond::DATIVER ){
    $2->setBeginAtomIdx(atomIdx1);
    $2->setEndAtomIdx(atomIdx2);
    $2->setBondType(Bond::DATIVE);
  }else if ( $2->getBondType() == Bond::DATIVEL ){
    $2->setBeginAtomIdx(atomIdx2);
    $2->setEndAtomIdx(atomIdx1);
    $2->setBondType(Bond::DATIVE);
  } else {
    $2->setBeginAtomIdx(atomIdx1);
    $2->setEndAtomIdx(atomIdx2);
  }
  $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond($2);
  delete $2;
}

| mol SEPARATOR_TOKEN atomd {
  RWMol *mp = (*molList)[$$];
  $3->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom($3,true,true);
}

| mol ring_number {
  RWMol * mp = (*molList)[$$];
  Atom *atom=mp->getActiveAtom();

  QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(atom, nullptr);
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atom->getIdx());
  mp->setBondBookmark(newB,$2);
  if(!(mp->getAllBondsWithBookmark($2).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }
  mp->setAtomBookmark(atom,$2);

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  if(atom->hasProp(RDKit::common_properties::_RingClosures)){
    atom->getProp(RDKit::common_properties::_RingClosures,tmp);
  }
  tmp.push_back(-($2+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);

}

| mol bond_expr ring_number {
  RWMol * mp = (*molList)[$$];
  Atom *atom=mp->getActiveAtom();

  mp->setBondBookmark($2,$3);
  $2->setOwningMol(mp);
  $2->setBeginAtomIdx(atom->getIdx());
  $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->setAtomBookmark(atom,$3);

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  if(atom->hasProp(RDKit::common_properties::_RingClosures)){
    atom->getProp(RDKit::common_properties::_RingClosures,tmp);
  }
  tmp.push_back(-($3+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);

}

| mol GROUP_OPEN_TOKEN atomd {
  RWMol *mp = (*molList)[$$];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom($3,true,true);

  QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(a1,mp->getAtomWithIdx(atomIdx2));
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond(newB);
  delete newB;

  branchPoints->push_back(atomIdx1);
}

| mol GROUP_OPEN_TOKEN bond_expr atomd  {
  RWMol *mp = (*molList)[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($4,true,true);
  if( $3->getBondType() == Bond::DATIVER ){
    $3->setBeginAtomIdx(atomIdx1);
    $3->setEndAtomIdx(atomIdx2);
    $3->setBondType(Bond::DATIVE);
  }else if ( $3->getBondType() == Bond::DATIVEL ){
    $3->setBeginAtomIdx(atomIdx2);
    $3->setEndAtomIdx(atomIdx1);
    $3->setBondType(Bond::DATIVE);
  } else {
    $3->setBeginAtomIdx(atomIdx1);
    $3->setEndAtomIdx(atomIdx2);
  }
  $3->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond($3,true);
  branchPoints->push_back(atomIdx1);

}


| mol GROUP_CLOSE_TOKEN {
  if(branchPoints->empty()){
     yyerror(input,molList,branchPoints,scanner,start_token,"extra close parentheses");
     yyErrorCleanup(molList);
     YYABORT;
  }
  RWMol *mp = (*molList)[$$];
  mp->setActiveAtom(branchPoints->back());
  branchPoints->pop_back();
}

;

/* --------------------------------------------------------------- */
atomd:	simple_atom
| hydrogen_atom
| ATOM_OPEN_TOKEN atom_expr ATOM_CLOSE_TOKEN
{
  $$ = $2;
}
| ATOM_OPEN_TOKEN atom_expr COLON_TOKEN number ATOM_CLOSE_TOKEN
{
  $$ = $2;
  $$->setProp(RDKit::common_properties::molAtomMapNumber,$4);
}
;


/* --------------------------------------------------------------- */
/*
Some ugliness here around how Hs are handled. This is due to this paragraph
in the SMIRKS docs, of all places:
http://www.daylight.com/dayhtml/doc/theory/theory.smirks.html

Hence, a single change to SMARTS interpretation, for expressions of the form:
[<weight>]H<charge><map>]. In SMARTS, these expressions now are interpreted as
a hydrogen atom, rather than as any atom with one hydrogen attached.
All other SMARTS hydrogen expressions retain their pre-4.51 meanings.

Thanks to John Mayfield for pointing this out.
*/

/* --------------------------------------------------------------- */
hydrogen_atom:	ATOM_OPEN_TOKEN H_TOKEN ATOM_CLOSE_TOKEN
{
  $$ = new QueryAtom(1);
}
| ATOM_OPEN_TOKEN H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN
{
  $$ = new QueryAtom(1);
  $$->setProp(RDKit::common_properties::molAtomMapNumber,$4);
}
| ATOM_OPEN_TOKEN number H_TOKEN ATOM_CLOSE_TOKEN {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope($2);
  newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
  $$=newQ;
}
| ATOM_OPEN_TOKEN number H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope($2);
  newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,$5);

  $$=newQ;
}

| ATOM_OPEN_TOKEN H_TOKEN charge_spec ATOM_CLOSE_TOKEN {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge($3);
  newQ->expandQuery(makeAtomFormalChargeQuery($3),Queries::COMPOSITE_AND,true);
  $$=newQ;
}
| ATOM_OPEN_TOKEN H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge($3);
  newQ->expandQuery(makeAtomFormalChargeQuery($3),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,$5);

  $$=newQ;
}
| ATOM_OPEN_TOKEN number H_TOKEN charge_spec ATOM_CLOSE_TOKEN {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope($2);
  newQ->setFormalCharge($4);
  newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery($4),Queries::COMPOSITE_AND,true);
  $$=newQ;
}
| ATOM_OPEN_TOKEN number H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope($2);
  newQ->setFormalCharge($4);
  newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery($4),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,$6);

  $$=newQ;
}
;

/* --------------------------------------------------------------- */
atom_expr: atom_expr AND_TOKEN atom_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if($1->getChiralTag()==Atom::CHI_UNSPECIFIED) $1->setChiralTag($3->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps($1);
  delete $3;
}
| atom_expr OR_TOKEN atom_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if($1->getChiralTag()==Atom::CHI_UNSPECIFIED) $1->setChiralTag($3->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps($1);
  $1->setAtomicNum(0);
  delete $3;
}
| atom_expr SEMI_TOKEN atom_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if($1->getChiralTag()==Atom::CHI_UNSPECIFIED) $1->setChiralTag($3->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps($1);
  delete $3;
}
| atom_expr point_query {
  $1->expandQuery($2->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if($1->getChiralTag()==Atom::CHI_UNSPECIFIED) $1->setChiralTag($2->getChiralTag());
  if($2->getNumExplicitHs()){
    if(!$1->getNumExplicitHs()){
      $1->setNumExplicitHs($2->getNumExplicitHs());
      $1->setNoImplicit(true);
    } else if($1->getNumExplicitHs()!=$2->getNumExplicitHs()){
      // conflicting queries...
      $1->setNumExplicitHs(0);
      $1->setNoImplicit(false);
    }
  }
  if($2->getFormalCharge()){
    if(!$1->getFormalCharge()){
      $1->setFormalCharge($2->getFormalCharge());
    } else if($1->getFormalCharge()!=$2->getFormalCharge()){
      // conflicting queries...
      $1->setFormalCharge(0);
    }
  }
  delete $2;
}
| point_query
;

point_query: NOT_TOKEN point_query {
  $2->getQuery()->setNegation(!($2->getQuery()->getNegation()));
  $2->setAtomicNum(0);
  SmilesParseOps::ClearAtomChemicalProps($2);
  $$ = $2;
}
| recursive_query
| atom_query
;

/* --------------------------------------------------------------- */
recursive_query: BEGIN_RECURSE mol END_RECURSE {
  // this is a recursive SMARTS expression
  QueryAtom *qA = new QueryAtom();
  //  FIX: there's maybe a leak here
  RWMol *molP = (*molList)[$2];
  // close any rings in the molecule:
  SmilesParseOps::CloseMolRings(molP,0);

  //molP->debugMol(std::cout);
  qA->setQuery(new RecursiveStructureQuery(molP));
  //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
  int sz = molList->size();
  if ( sz==$2+1) {
    molList->resize( sz-1 );
  }
  $$ = qA;
}
| BEGIN_RECURSE mol END_RECURSE UNDERSCORE_TOKEN  nonzero_number{
  // UNDOCUMENTED EXTENSION:
  // this is a recursive SMARTS expression with a serial number
  // please don't write your own SMARTS that include this extension:
  // the RDKit smarts parsing code will automatically insert serial
  // numbers for recursive smarts patterns.
  QueryAtom *qA = new QueryAtom();
  //  FIX: there's maybe a leak here
  RWMol *molP = (*molList)[$2];
  // close any rings in the molecule:
  SmilesParseOps::CloseMolRings(molP,0);

  //molP->debugMol(std::cout);
  qA->setQuery(new RecursiveStructureQuery(molP,$5));
  //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
  int sz = molList->size();
  if ( sz==$2+1) {
    molList->resize( sz-1 );
  }
  $$ = qA;
}
;

/* --------------------------------------------------------------- */
atom_query:	simple_atom
| number simple_atom {
  $2->setIsotope($1);
  $2->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
  $$=$2;
}
| ATOM_TOKEN
| number ATOM_TOKEN {
  $2->setIsotope($1);
  $2->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
  $$=$2;
}
| HASH_TOKEN number { $$ = new QueryAtom($2); }
| number HASH_TOKEN number {
  $$ = new QueryAtom($3);
  $$->setIsotope($1);
  $$->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
}
| COMPLEX_ATOM_QUERY_TOKEN
| HETERONEIGHBOR_ATOM_QUERY_TOKEN
| ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN
| RINGSIZE_ATOM_QUERY_TOKEN
| RINGBOND_ATOM_QUERY_TOKEN
| IMPLICIT_H_ATOM_QUERY_TOKEN
| COMPLEX_ATOM_QUERY_TOKEN number {
  static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
}
| HETERONEIGHBOR_ATOM_QUERY_TOKEN number {
  $1->setQuery(makeAtomNumHeteroatomNbrsQuery($2));
}
| ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN number {
  $1->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery($2));
}
| RINGSIZE_ATOM_QUERY_TOKEN number {
  $1->setQuery(makeAtomMinRingSizeQuery($2));
}
| RINGBOND_ATOM_QUERY_TOKEN number {
  $1->setQuery(makeAtomRingBondCountQuery($2));
}
| IMPLICIT_H_ATOM_QUERY_TOKEN number {
  $1->setQuery(makeAtomImplicitHCountQuery($2));
}
| possible_range_query RANGE_OPEN_TOKEN MINUS_TOKEN number RANGE_CLOSE_TOKEN {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>($1->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>($4,oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  $1->setQuery(nq);
}
| possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN RANGE_CLOSE_TOKEN {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>($1->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>($3,oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  $1->setQuery(nq);
}
| possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN number RANGE_CLOSE_TOKEN {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>($1->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery($3,$5,false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  $1->setQuery(nq);
}
| number H_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery($1));
  newQ->setIsotope($1);
  newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs(1);
  $$=newQ;
}
| number H_TOKEN number {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery($1));
  newQ->setIsotope($1);
  newQ->expandQuery(makeAtomHCountQuery($3),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs($3);
  $$=newQ;
}
| H_TOKEN number {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery($2));
  newQ->setNumExplicitHs($2);
  $$=newQ;
}
| H_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  $$=newQ;
}
| charge_spec {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery($1));
  newQ->setFormalCharge($1);
  $$=newQ;
}
| AT_TOKEN AT_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  $$=newQ;
}
| AT_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  $$=newQ;
}
| CHI_CLASS_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag($1);
  newQ->setProp(common_properties::_chiralPermutation,0);
  $$=newQ;
}
| CHI_CLASS_TOKEN number {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag($1);
  newQ->setProp(common_properties::_chiralPermutation,$2);
  $$=newQ;
}
| HYB_TOKEN
| number {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery($1));
  $$=newQ;
}
;

possible_range_query : COMPLEX_ATOM_QUERY_TOKEN
| HETERONEIGHBOR_ATOM_QUERY_TOKEN {
  $1->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
}
| ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN {
  $1->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
}
| RINGSIZE_ATOM_QUERY_TOKEN {
  $1->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
}
| RINGBOND_ATOM_QUERY_TOKEN {
  $1->setQuery(makeAtomRingBondCountQuery(0));
}
| IMPLICIT_H_ATOM_QUERY_TOKEN {
  $1->setQuery(makeAtomImplicitHCountQuery(0));
}
| PLUS_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery(0));
  $$ = newQ;
}
| MINUS_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNegativeFormalChargeQuery(0));
  $$ = newQ;
}
;

/* --------------------------------------------------------------- */
simple_atom: 	ORGANIC_ATOM_TOKEN {
  //
  // This construction (and some others) may seem odd, but the
  // SMARTS definition requires that an atom which is aliphatic on
  // input (i.e. something in the "organic subset" that is given with
  // a capital letter) only match aliphatic atoms.
  //
  // The following rule applies a similar logic to aromatic atoms.
  //
  $$ = new QueryAtom($1);
  $$->setQuery(makeAtomTypeQuery($1,false));
}
| AROMATIC_ATOM_TOKEN {
  $$ = new QueryAtom($1);
  $$->setIsAromatic(true);
  $$->setQuery(makeAtomTypeQuery($1,true));
}
| SIMPLE_ATOM_QUERY_TOKEN
;


/* --------------------------------------------------------------- */
bond_expr:bond_expr AND_TOKEN bond_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete $3;
}
| bond_expr OR_TOKEN bond_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete $3;
}
| bond_expr SEMI_TOKEN bond_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete $3;
}
| bond_query
;

bond_query: bondd
| bond_query bondd {
  $1->expandQuery($2->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete $2;
}
;

/* --------------------------------------------------------------- */
bondd: BOND_TOKEN
| MINUS_TOKEN {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  $$ = newB;
}
| HASH_TOKEN {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  $$ = newB;
}
| COLON_TOKEN {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  $$ = newB;
}
| AT_TOKEN {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  $$ = newB;
}
| NOT_TOKEN bondd {
  $2->getQuery()->setNegation(!($2->getQuery()->getNegation()));
  $$ = $2;
}
;

/* --------------------------------------------------------------- */
charge_spec: PLUS_TOKEN PLUS_TOKEN { $$=2; }
| PLUS_TOKEN number { $$=$2; }
| PLUS_TOKEN { $$=1; }
| MINUS_TOKEN MINUS_TOKEN { $$=-2; }
| MINUS_TOKEN number { $$=-$2; }
| MINUS_TOKEN { $$=-1; }
;

/* --------------------------------------------------------------- */
ring_number:  digit
| PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit { $$ = $2*10+$3; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN { $$ = $3; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN { $$ = $3*10+$4; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*100+$4*10+$5; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*1000+$4*100+$5*10+$6; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*10000+$4*1000+$5*100+$6*10+$7; }
;


/* --------------------------------------------------------------- */
number:  ZERO_TOKEN
| nonzero_number
;

/* --------------------------------------------------------------- */
nonzero_number:  NONZERO_DIGIT_TOKEN
| nonzero_number digit { 
    if($1 >= std::numeric_limits<std::int32_t>::max()/10 || 
     $1*10 >= std::numeric_limits<std::int32_t>::max()-$2 ){
     yysmarts_error(input,molList,lastAtom,lastBond,numAtomsParsed,numBondsParsed,branchPoints,scanner,start_token,"number too large");
     YYABORT;
  }
  $$ = $1*10 + $2; }
;

digit: NONZERO_DIGIT_TOKEN
| ZERO_TOKEN
;

%%
