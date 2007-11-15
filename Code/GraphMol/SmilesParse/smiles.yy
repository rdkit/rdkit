%{

  // $Id: smiles.y 4974 2006-02-18 00:49:21Z glandrum $
  //
  //  Copyright (C) 2001-2006 Randal Henne, Greg Landrum and Rational Discovery LLC
  //
  //   @@ All Rights Reserved  @@
  //

#include <iostream>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include "SmilesParseOps.h"  
#include <RDGeneral/RDLog.h>

extern int yysmiles_lex();

#define YYDEBUG 1

void
yysmiles_error( const char * msg )
{

}

using namespace RDKit;

namespace RDKit {
 namespace SmilesParse{
  extern std::vector<RDKit::RWMol *> molList_g;
 }
}
static RWMol * curMol_gps = 0;


%}
 
 
%union {
  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  int                      ival;
}

%token <atom> AROMATIC_ATOM_TOKEN ATOM_TOKEN ORGANIC_ATOM_TOKEN
%token <ival> DIGIT_TOKEN
%token GROUP_OPEN_TOKEN GROUP_CLOSE_TOKEN SEPARATOR_TOKEN LOOP_CONNECTOR_TOKEN
%token MINUS_TOKEN PLUS_TOKEN CHIRAL_MARKER_TOKEN CHI_CLASS_TOKEN CHI_CLASS_OH_TOKEN
%token H_TOKEN AT_TOKEN PERCENT_TOKEN
%token <bond> BOND_TOKEN
%type <moli> cmpd mol branch
%type <atom> atomd element chiral_element h_element charge_element simple_atom
%type <ival>  number ring_number
%token ATOM_OPEN_TOKEN ATOM_CLOSE_TOKEN
%token EOS_TOKEN

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
  BOOST_LOG(rdErrorLog) << "SMILES Parse Error" << std::endl;
  SmilesParse::molList_g.clear();
  SmilesParse::molList_g.resize(0);
  YYABORT;
}
| cmpd EOS_TOKEN {
  YYACCEPT;
}
| error EOS_TOKEN {
  yyclearin;
  yyerrok;
  BOOST_LOG(rdErrorLog) << "SMILES Parse Error" << std::endl;
  SmilesParse::molList_g.clear();
  SmilesParse::molList_g.resize(0);
  YYABORT;
}
;

/* --------------------------------------------------------------- */
// FIX: mol MINUS DIGIT
mol: atomd {
  int sz     = SmilesParse::molList_g.size();
  SmilesParse::molList_g.resize( sz + 1);
  SmilesParse::molList_g[ sz ] = new RWMol();
  curMol_gps = SmilesParse::molList_g[ sz ];
  curMol_gps->addAtom($1);
  delete $1;
  $$ = sz;
}

| mol atomd       {
  RWMol *mp = SmilesParse::molList_g[$$];
  RWMol::GRAPH_NODE_TYPE a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom($2);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  delete $2;
}

| mol BOND_TOKEN atomd  {
  RWMol *mp = SmilesParse::molList_g[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($3);
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
  mp->addBond($2,true);
  delete $3;
}

| mol MINUS_TOKEN atomd {
  RWMol *mp = SmilesParse::molList_g[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($3);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  delete $3;
}

| mol ring_number {
  RWMol * mp = SmilesParse::molList_g[$$];
  Atom *atom=mp->getActiveAtom();
  mp->setAtomBookmark(atom,$2);

  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::UNSPECIFIED);
  mp->setBondBookmark(newB,$2);
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-($2+1));
  atom->setProp("_RingClosures",tmp);
}

| mol BOND_TOKEN ring_number {
  RWMol * mp = SmilesParse::molList_g[$$];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     $2->getBondType());
  newB->setBondDir($2->getBondDir());
  mp->setAtomBookmark(atom,$3);
  mp->setBondBookmark(newB,$3);
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-($3+1));
  atom->setProp("_RingClosures",tmp);
  delete $2;
}

| mol MINUS_TOKEN ring_number {
  RWMol * mp = SmilesParse::molList_g[$$];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::SINGLE);
  mp->setAtomBookmark(atom,$3);
  mp->setBondBookmark(newB,$3);
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-($3+1));
  atom->setProp("_RingClosures",tmp);
}

| mol branch {
  RWMol *m1_p = SmilesParse::molList_g[$$],*m2_p=SmilesParse::molList_g[$2];
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED,Bond::NONE,false);
  delete m2_p;
  int sz = SmilesParse::molList_g.size();
  if ( sz==$2+1) {
    SmilesParse::molList_g.resize( sz-1 );
  }
}
;
/* --------------------------------------------------------------- */
branch:	GROUP_OPEN_TOKEN mol GROUP_CLOSE_TOKEN { $$ = $2; }
| GROUP_OPEN_TOKEN BOND_TOKEN mol GROUP_CLOSE_TOKEN {
  $$ = $3;
  int sz     = SmilesParse::molList_g.size();
  curMol_gps = SmilesParse::molList_g[ sz-1 ];

  Bond *partialBond = curMol_gps->createPartialBond(0,$2->getBondType());
  partialBond->setBondDir($2->getBondDir());
  curMol_gps->setBondBookmark(partialBond,
			      ci_LEADING_BOND);
  delete $2;
}
| GROUP_OPEN_TOKEN MINUS_TOKEN mol GROUP_CLOSE_TOKEN {
  $$ = $3;
  int sz     = SmilesParse::molList_g.size();
  curMol_gps = SmilesParse::molList_g[ sz-1 ];

  Bond *partialBond = curMol_gps->createPartialBond(0,Bond::SINGLE);
  curMol_gps->setBondBookmark(partialBond,
			      ci_LEADING_BOND);
}
;

/* --------------------------------------------------------------- */
atomd:	simple_atom
| ATOM_OPEN_TOKEN charge_element ATOM_CLOSE_TOKEN
{
  $$ = $2;
  $2->setNoImplicit(true);
}
;

/* --------------------------------------------------------------- */
charge_element:	h_element
| h_element PLUS_TOKEN { $1->setFormalCharge(1); }
| h_element PLUS_TOKEN PLUS_TOKEN { $1->setFormalCharge(2); }
| h_element PLUS_TOKEN number { $1->setFormalCharge($3); }
| h_element MINUS_TOKEN { $1->setFormalCharge(-1); }
| h_element MINUS_TOKEN MINUS_TOKEN { $1->setFormalCharge(-2); }
| h_element MINUS_TOKEN number { $1->setFormalCharge(-$3); }
		;

/* --------------------------------------------------------------- */
h_element:      H_TOKEN { $$ = new Atom(1); }
                | number H_TOKEN { $$ = new Atom(1); $$->setMass($1); }
                | chiral_element
		| chiral_element H_TOKEN		{ $$ = $1; $1->setNumExplicitHs(1);}
		| chiral_element H_TOKEN number	{ $$ = $1; $1->setNumExplicitHs($3);}
		;

/* --------------------------------------------------------------- */
chiral_element:	 element
| element AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
| element AT_TOKEN AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
;

/* --------------------------------------------------------------- */
element:	simple_atom
		|	number simple_atom { $2->setMass( $1 ); $$ = $2; }
		|	ATOM_TOKEN
		|	number ATOM_TOKEN	   { $2->setMass( $1 ); $$ = $2; }
		;

/* --------------------------------------------------------------- */
simple_atom:      ORGANIC_ATOM_TOKEN
                | AROMATIC_ATOM_TOKEN
                ;


/* --------------------------------------------------------------- */
ring_number:  DIGIT_TOKEN
| PERCENT_TOKEN DIGIT_TOKEN DIGIT_TOKEN { $$ = $2*10 + $3; }
;

/* --------------------------------------------------------------- */
number:  DIGIT_TOKEN
| number DIGIT_TOKEN { $$ = $1*10 + $2; }
;


/*
  chival:	CHI_CLASS_TOKEN DIGIT_TOKEN
	| AT_TOKEN
        ;
*/


%%


