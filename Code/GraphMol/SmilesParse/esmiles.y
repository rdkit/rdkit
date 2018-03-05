%{

//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#pragma warning(disable:4786)

#include <iostream>
#include <string>


#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <RDGeneral/Invariant.h>
#include "SmilesParseOps.h"  

#include <vector>

using namespace std;

extern int yylex();

#define YYDEBUG 1

void
yyerror( const char * msg )
{
  std::cerr << msg << "\n";
}

using namespace RDKit;

extern vector<Mol *> molList_g;
static Mol * curMol_gps = 0;

int label_gs = -1;


%}
 
 
%union {
  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  int                      ival;
}

%token <atom> AROMATIC_ATOM ATOM ORGANIC_ATOM
%token <ival> DIGIT
%token GROUP_OPEN GROUP_CLOSE SEPARATOR LOOP_CONNECTOR
%token MINUS PLUS CHIRAL_MARKER CHI_CLASS CHI_CLASS_OH
%token H AT PERCENT
%token <bond> BOND

#if 0
%type <moli> mol branch hapto_branch dbl_hapto_branch
%type <atom> atomd general_atom element chiral_element h_element charge_element
%type <ival>  number
#else
%type <moli> mol branch hapto_branch dbl_hapto_branch
%type <atom> atomd general_atom element chiral_element h_element charge_element
%type <ival>  number
#endif


%token ATOM_OPEN ATOM_CLOSE
%token LT GT DATIVE_MARK

//%left SEPARATOR GROUP_CLOSE DIGIT
//%right GROUP_OPEN
 
%%

// FIX: mol MINUS DIGIT, mol LT DIGIT, mol GT DIGIT

mol: general_atom {
  int sz     = molList_g.size();
  molList_g.resize( sz + 1);
  molList_g[ sz ] = new Mol();
  curMol_gps = molList_g[ sz ];
  curMol_gps->addAtom($1);
  delete $1;
  $$ = sz;
}

| mol general_atom       {
  Mol *mp = molList_g[$$];
  Mol::GRAPH_NODE_TYPE a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom($2);
  if(a1->getIsAromatic() && $2->getIsAromatic())
    mp->addBond(atomIdx1,atomIdx2,Bond::AROMATIC);
  else
    mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  delete $2;
}

| mol BOND general_atom  {
  Mol *mp = molList_g[$$];
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
  mp->addBond($2);

  delete $2;
  delete $3;
}

| mol MINUS general_atom {
  Mol *mp = molList_g[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($3);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  delete $3;
}

| mol LT general_atom {
  Mol *mp = molList_g[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($3);
  mp->addBond(atomIdx2,atomIdx1,Bond::DATIVE);
  delete $3;
}

| mol GT general_atom {
  Mol *mp = molList_g[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($3);
  mp->addBond(atomIdx1,atomIdx2,Bond::DATIVE);
  delete $3;
}

| mol number {
  Mol * mp = molList_g[$$];
  mp->setAtomBookmark(mp->getActiveAtom(),$2);
}

| mol BOND number {
  Mol * mp = molList_g[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  Bond *newB_sp = mp->createPartialBond(atomIdx1,
					$2->getBondType());
  newB_sp->setBondDir($2->getBondDir());
  mp->setAtomBookmark(mp->getActiveAtom(),$3);
  mp->setBondBookmark(newB_sp,$3);

  delete $2;
}

| mol branch {
  Mol *m1_p = molList_g[$$],*m2_p=molList_g[$2];
  AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED);
  delete m2_p;
  int sz = molList_g.size();
  if ( sz==$2+1) {
    molList_g.resize( sz-1 );
  }
}

| mol BOND branch {
  // FIX: should check that there's no partial bond already formed
  // in the branch!
  Mol *m1_p = molList_g[$$],*m2_p=molList_g[$3];
  AddFragToMol(m1_p,m2_p,$2->getBondType(),$2->getBondDir());

  delete m2_p;
  delete $2;

  int sz = molList_g.size();
  if ( sz==$3+1) {
    molList_g.resize( sz-1 );
  }
}

| mol MINUS branch {
  // FIX: should check that there's no partial bond already formed
  // in the branch!
  Mol *m1_p = molList_g[$$],*m2_p=molList_g[$3];
  AddFragToMol(m1_p,m2_p,Bond::SINGLE);

  delete m2_p;
  int sz = molList_g.size();
  if ( sz==$3+1) {
    molList_g.resize( sz-1 );
  }
}

| mol LT branch {
  // FIX: should check that there's no partial bond already formed
  // in the branch!
  Mol *m1_p = molList_g[$$],*m2_p=molList_g[$3];
  AddFragToMol(m1_p,m2_p,Bond::DATIVEL);

  delete m2_p;
  int sz = molList_g.size();
  if ( sz==$3+1) {
    molList_g.resize( sz-1 );
  }
}

| mol GT branch {
  // FIX: should check that there's no partial bond already formed
  // in the branch!
  Mol *m1_p = molList_g[$$],*m2_p=molList_g[$3];
  AddFragToMol(m1_p,m2_p,Bond::DATIVER);

  delete m2_p;
  int sz = molList_g.size();
  if ( sz==$3+1) {
    molList_g.resize( sz-1 );
  }
}

| mol SEPARATOR mol {
  Mol *m1_p = molList_g[$1],*m2_p=molList_g[$3];
  m1_p->insertMol(m2_p);
  delete m2_p;
  int sz = molList_g.size();
  if ( sz==$3+1) {
    molList_g.resize( sz-1 );
  }
}

| mol LT dbl_hapto_branch GT mol {
  Mol *m1_p = molList_g[$1],*m2_p=molList_g[$3],*m3_p=molList_g[$5];

  AddDblHaptoBranch(m1_p,m2_p,m3_p,true,true,1,1);
  
  delete m2_p;
  molList_g.erase(molList_g.begin()+$3);
}

| mol LT hapto_branch {
  Mol *m1_p = molList_g[$1],*m2_p=molList_g[$3];

  AddHaptoBranch(m1_p,m2_p,true,true,1);
  
  delete m2_p;
  int sz = molList_g.size();
  if ( sz==$3+1) {
    molList_g.resize( sz-1 );
  }
}

| hapto_branch GT mol {
  Mol *m1_p = molList_g[$3],*m2_p=molList_g[$1];

  AddHaptoBranch(m1_p,m2_p,false,true,1);
  
  delete m2_p;
  int sz = molList_g.size();
  molList_g.erase(molList_g.begin()+$1);
}

| mol GROUP_OPEN LT hapto_branch GROUP_CLOSE {
  Mol *m1_p = molList_g[$1],*m2_p=molList_g[$4];

  AddHaptoBranch(m1_p,m2_p,false,true,1);
  
  delete m2_p;
  int sz = molList_g.size();
  if ( sz==$4+1) {
    molList_g.resize( sz-1 );
  }
}

| error {
  yyclearin;
  std::cerr << "SMILES Parse Error" << std::endl;
  molList_g.clear();
  molList_g.resize(0);
  YYABORT;
}

; 


number:  DIGIT
| PERCENT DIGIT DIGIT { $$ = $2*10 + $3; }
;


element:	ORGANIC_ATOM
                | AROMATIC_ATOM
		|	DIGIT ORGANIC_ATOM { $2->setMass( $1 ); $$ = $2; }
		|	ATOM
		|	DIGIT ATOM	   { $2->setMass( $1 ); $$ = $2; }
                |	H { $$ = new Atom(1); }
		|	DIGIT H	   { $$ = new Atom(1); $$->setMass( $1 ); }


		;

chiral_element:	element			
				| element chival
				;

h_element:	chiral_element      
		| chiral_element H		{ $$ = $1; $1->setNumExplicitHs(1);}
		| chiral_element H DIGIT	{ $$ = $1; $1->setNumExplicitHs($3);}
		;

/*sign:		PLUS
		| MINUS
		;
*/

charge_element:	h_element
| h_element PLUS { $1->setFormalCharge(1); }
| h_element PLUS PLUS { $1->setFormalCharge(2); }
| h_element MINUS { $1->setFormalCharge(-1); }
| h_element MINUS MINUS { $1->setFormalCharge(-2); }
| h_element PLUS DIGIT { $1->setFormalCharge($3); }
| h_element MINUS DIGIT { $1->setFormalCharge(-$3); }
/*| h_element sign DIGIT*/
		;


general_atom: atomd
{
  $$ = $1;
}
| atomd DATIVE_MARK
{
  $$ = $1;
  $1->setDativeFlag(1);
}
;



atomd:	ORGANIC_ATOM
{
  $$ = $1;
  label_gs = -1;
}
| AROMATIC_ATOM
{
  $$ = $1;
  label_gs = -1;
}
| ATOM_OPEN charge_element ATOM_CLOSE
{
  $$ = $2;
  label_gs = -1;
  $2->setNoImplicit(true);
}
	;

chival:	CHI_CLASS DIGIT
	| AT
        ;

dbl_hapto_branch: hapto_branch DATIVE_MARK {
  $$ = $1;
}
| dbl_hapto_branch GROUP_OPEN GT mol GROUP_CLOSE{
  Mol *branch = molList_g[$1],*mol=molList_g[$4];
  Atom *origActive = mol->getActiveAtom();
  mol->setActiveAtom(0);
  AddHaptoBranch(mol,branch,true,false,1);
  
  $$ = $4-1;
  molList_g.erase(molList_g.begin()+$1);
}
;

hapto_branch: DATIVE_MARK GROUP_OPEN mol GROUP_CLOSE {
  $$ = $3;
}
| GROUP_OPEN mol GROUP_CLOSE DATIVE_MARK {
  $$ = $2;
}
;

branch:	GROUP_OPEN mol GROUP_CLOSE { $$ = $2; }
| GROUP_OPEN BOND mol GROUP_CLOSE {
  $$ = $3;
  int sz     = molList_g.size();
  curMol_gps = molList_g[ sz-1 ];

  Bond *partialBond = curMol_gps->createPartialBond(0,$2->getBondType());
  partialBond->setBondDir($2->getBondDir());
  curMol_gps->setBondBookmark(partialBond,
			      ci_LEADING_BOND);
  delete $2;
}
| GROUP_OPEN MINUS mol GROUP_CLOSE {
  $$ = $3;
  int sz     = molList_g.size();
  curMol_gps = molList_g[ sz-1 ];

  Bond *partialBond = curMol_gps->createPartialBond(0,Bond::SINGLE);
  curMol_gps->setBondBookmark(partialBond,
			      ci_LEADING_BOND);
}
| GROUP_OPEN GT mol GROUP_CLOSE {
  $$ = $3;
  int sz     = molList_g.size();
  curMol_gps = molList_g[ sz-1 ];

  Bond *partialBond = curMol_gps->createPartialBond(0,Bond::DATIVER);
  curMol_gps->setBondBookmark(partialBond,
			      ci_LEADING_BOND);
}
| GROUP_OPEN LT mol GROUP_CLOSE {
  $$ = $3;
  int sz     = molList_g.size();
  curMol_gps = molList_g[ sz-1 ];

  Bond *partialBond = curMol_gps->createPartialBond(0,Bond::DATIVEL);
  curMol_gps->setBondBookmark(partialBond,
			      ci_LEADING_BOND);
}
	;

%%


