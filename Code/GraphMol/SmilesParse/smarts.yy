%{

  // $Id$
  //
  //  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
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

extern int yysmarts_lex();

#define YYDEBUG 1

void
yysmarts_error( const char * msg )
{

}

using namespace RDKit;

namespace RDKit {
  namespace SmilesParse {  
    extern std::vector<RWMol *> molList_g;
  }
}

 void yyErrorCleanup(){
  for(std::vector<RDKit::RWMol *>::iterator iter=SmilesParse::molList_g.begin();
	         iter != SmilesParse::molList_g.end(); ++iter){
     delete *iter;
  }
  
  SmilesParse::molList_g.clear();
  SmilesParse::molList_g.resize(0);

 }

%}
 
 
%union {
  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;
}

%token <atom> AROMATIC_ATOM_TOKEN ATOM_TOKEN ORGANIC_ATOM_TOKEN
%token <atom> SIMPLE_ATOM_QUERY_TOKEN COMPLEX_ATOM_QUERY_TOKEN RINGSIZE_ATOM_QUERY_TOKEN
%token <ival> DIGIT_TOKEN
%token GROUP_OPEN_TOKEN GROUP_CLOSE_TOKEN SEPARATOR_TOKEN
%token HASH_TOKEN MINUS_TOKEN PLUS_TOKEN 
%token CHIRAL_MARKER_TOKEN CHI_CLASS_TOKEN CHI_CLASS_OH_TOKEN
%token H_TOKEN AT_TOKEN PERCENT_TOKEN
%token ATOM_OPEN_TOKEN ATOM_CLOSE_TOKEN 
%token NOT_TOKEN AND_TOKEN OR_TOKEN SEMI_TOKEN BEGIN_RECURSE END_RECURSE
%token HYB_TOKEN COLON_TOKEN
%token <bond> BOND_TOKEN
%type <moli> cmpd mol branch
%type <atom> atomd element simple_atom
%type <atom> atom_expr point_query atom_query recursive_query
%type <ival> ring_number number charge_spec
%type <bond> bondd bond_expr bond_query
%token EOS_TOKEN

%left SEMI_TOKEN
%left OR_TOKEN
%left AND_TOKEN
%right NOT_TOKEN

%%

/* --------------------------------------------------------------- */
cmpd: mol
| cmpd SEPARATOR_TOKEN mol {
  RWMol *m1_p = SmilesParse::molList_g[$1],*m2_p=SmilesParse::molList_g[$3];
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::IONIC,Bond::NONE,true);
  delete m2_p;
  int sz = SmilesParse::molList_g.size();
  if ( sz==$3+1) {
    SmilesParse::molList_g.resize( sz-1 );
  }
}
| cmpd error EOS_TOKEN{
  yyclearin;
  yyerrok;
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error" << std::endl;
  yyErrorCleanup();
  YYABORT;
}
| cmpd EOS_TOKEN {
  YYACCEPT;
}
| error EOS_TOKEN{
  yyclearin;
  yyerrok;
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error" << std::endl;
  
  yyErrorCleanup();
  YYABORT;
}
;

/* --------------------------------------------------------------- */
// FIX: mol MINUS DIGIT
mol: atomd {
  int sz     = SmilesParse::molList_g.size();
  SmilesParse::molList_g.resize( sz + 1);
  SmilesParse::molList_g[ sz ] = new RWMol();
  SmilesParse::molList_g[ sz ]->addAtom($1,true,true);
  //delete $1;
  $$ = sz;
}
| mol atomd       {
  RWMol *mp = SmilesParse::molList_g[$$];
  RWMol::GRAPH_NODE_TYPE a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom($2,true,true);

  QueryBond *newB;
  // this is a bit of a hack to try and get nicer "SMILES" from 
  // a SMARTS molecule:
  if(!(a1->getIsAromatic() && $2->getIsAromatic())){
    newB = new QueryBond(Bond::SINGLE);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
  			    Queries::COMPOSITE_OR,
  			    true);
  } else {
    newB = new QueryBond(Bond::AROMATIC);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::SINGLE),
  			    Queries::COMPOSITE_OR,
  			    true);
  }
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  mp->addBond(newB);
  delete newB;
  //delete $2;
}

| mol bond_expr atomd  {
  RWMol *mp = SmilesParse::molList_g[$$];
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
  mp->addBond($2);
  delete $2;

}

| mol ring_number {
  RWMol * mp = SmilesParse::molList_g[$$];
  Atom *atom=mp->getActiveAtom();

  // this is a bit of a hack to try and get nicer "SMILES" from 
  // a SMARTS molecule:
  QueryBond * newB;
  if(!atom->getIsAromatic()){
    newB = new QueryBond(Bond::SINGLE);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
  			    Queries::COMPOSITE_OR,
  			    true);
  } else {
    newB = new QueryBond(Bond::AROMATIC);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::SINGLE),
  			    Queries::COMPOSITE_OR,
  			    true);
  }
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atom->getIdx());
  mp->setBondBookmark(newB,$2);

  mp->setAtomBookmark(atom,$2);
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-($2+1));
  atom->setProp("_RingClosures",tmp);

}

| mol bond_expr ring_number {
  RWMol * mp = SmilesParse::molList_g[$$];
  Atom *atom=mp->getActiveAtom();

  mp->setBondBookmark($2,$3);
  $2->setOwningMol(mp);
  $2->setBeginAtomIdx(atom->getIdx());

  mp->setAtomBookmark(atom,$3);
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-($3+1));
  atom->setProp("_RingClosures",tmp);

}

| mol branch {
  RWMol *m1_p = SmilesParse::molList_g[$$],*m2_p=SmilesParse::molList_g[$2];
  // FIX: handle generic bonds here
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED,Bond::NONE,false,true);
  delete m2_p;
  int sz = SmilesParse::molList_g.size();
  if ( sz==$2+1) {
    SmilesParse::molList_g.resize( sz-1 );
  }
}
; 

/* --------------------------------------------------------------- */
branch:	GROUP_OPEN_TOKEN mol GROUP_CLOSE_TOKEN { $$ = $2; }
| GROUP_OPEN_TOKEN bond_expr mol GROUP_CLOSE_TOKEN {
  // FIX: this needs to handle arbitrary bond_exprs
  $$ = $3;
  int sz     = SmilesParse::molList_g.size();
  $2->setOwningMol(SmilesParse::molList_g[ sz-1 ]);
  $2->setBeginAtomIdx(0);
  SmilesParse::molList_g[ sz-1 ]->setBondBookmark($2,ci_LEADING_BOND);
}
;

/* --------------------------------------------------------------- */
atomd:	simple_atom
| ATOM_OPEN_TOKEN H_TOKEN ATOM_CLOSE_TOKEN
{
  $$ = new QueryAtom(1);
}
| ATOM_OPEN_TOKEN atom_expr ATOM_CLOSE_TOKEN
{
  $$ = $2;
}
;

/* --------------------------------------------------------------- */
atom_expr: atom_expr AND_TOKEN atom_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete $3;
}
| atom_expr OR_TOKEN atom_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete $3;
}
| atom_expr SEMI_TOKEN atom_expr {
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete $3;
}
| NOT_TOKEN atom_expr {
  $2->getQuery()->setNegation(!($2->getQuery()->getNegation()));
  $$ = $2;
}
| atom_expr NOT_TOKEN point_query {
  // FIX: this stuff (formerly element NOT_TOKEN atom_expr_piece) is making us shift-reduce-a-riffic
  $3->getQuery()->setNegation(!($3->getQuery()->getNegation()));
  $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete $3;
  }
| atom_expr point_query {
  // FIX: this stuff (formerly element atom_expr_piece) is making us shift-reduce-a-riffic
  $1->expandQuery($2->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete $2;
}
| atom_expr COLON_TOKEN number {
  if($1->hasProp("molAtomMapNumber")){
    BOOST_LOG(rdWarningLog) << "Warning: Atom-map index (:%d) specified multiple times for one atom, additional specifications ignored" << std::endl;
  } else {
    $1->setProp("molAtomMapNumber",$3);
  }
}
| point_query
| element AT_TOKEN AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); $$=$1;}
| element AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); $$=$1;}
| element
;

/* --------------------------------------------------------------- */
point_query: atom_query
| recursive_query
;

/* --------------------------------------------------------------- */
recursive_query: BEGIN_RECURSE mol END_RECURSE {
  // this is a recursive SMARTS expression
  QueryAtom *qA = new QueryAtom();
  //  FIX: there's maybe a leak here
  RWMol *molP = SmilesParse::molList_g[$2];
  // close any rings in the molecule:
  SmilesParseOps::CloseMolRings(molP,0);

  //molP->debugMol(std::cout);
  qA->setQuery(new RecursiveStructureQuery(molP));
  //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
  int sz = SmilesParse::molList_g.size();
  if ( sz==$2+1) {
    SmilesParse::molList_g.resize( sz-1 );
  }
  $$ = qA;
}
;

/* --------------------------------------------------------------- */
atom_query: COMPLEX_ATOM_QUERY_TOKEN
| COMPLEX_ATOM_QUERY_TOKEN number {
  static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
}
| RINGSIZE_ATOM_QUERY_TOKEN
| RINGSIZE_ATOM_QUERY_TOKEN number {
  delete $1->getQuery();
  $1->setQuery(makeAtomMinRingSizeQuery($2));
}
| charge_spec {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery($1));
  $$=newQ;
}
| H_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  $$=newQ;
}
| H_TOKEN number {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery($2));
  $$=newQ;
}
| HYB_TOKEN DIGIT_TOKEN {
  QueryAtom *newQ = new QueryAtom();
  Atom::HybridizationType hyb=Atom::UNSPECIFIED;
  // we're following the oelib "standard", where:
  // ^n implies that the atom is sp^n hybridized
  switch($2){
  case 0: hyb=Atom::S;break;
  case 1: hyb=Atom::SP;break;
  case 2: hyb=Atom::SP2;break;
  case 3: hyb=Atom::SP3;break;
  }
  newQ->setQuery(makeAtomHybridizationQuery(hyb));
  $$=newQ;
}
;

/* --------------------------------------------------------------- */
element: simple_atom
| ATOM_TOKEN
| HASH_TOKEN number { $$ = new QueryAtom($2); }
| number simple_atom {
  $2->expandQuery(makeAtomMassQuery($1),Queries::COMPOSITE_AND);
  $$ = $2;
}
| number ATOM_TOKEN {
  $2->expandQuery(makeAtomMassQuery($1),Queries::COMPOSITE_AND);
  $$ = $2;
}
| number HASH_TOKEN number {
  $$ = new QueryAtom($3);
  $$->expandQuery(makeAtomMassQuery($1),Queries::COMPOSITE_AND);
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
  $1->expandQuery(makeAtomAliphaticQuery(),Queries::COMPOSITE_AND);
}
| AROMATIC_ATOM_TOKEN{
  $1->expandQuery(makeAtomAromaticQuery(),Queries::COMPOSITE_AND);
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
/*
chival:	CHI_CLASS_TOKEN number
	| AT_TOKEN
        ;
*/
/* --------------------------------------------------------------- */
ring_number:  DIGIT_TOKEN
| PERCENT_TOKEN DIGIT_TOKEN DIGIT_TOKEN { $$ = $2*10 + $3; }
;

/* --------------------------------------------------------------- */
number:  DIGIT_TOKEN
| number DIGIT_TOKEN { $$ = $1*10 + $2; }
;

%%


