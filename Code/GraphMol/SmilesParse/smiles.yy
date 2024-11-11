%{

  // $Id$
  //
  //  Copyright (C) 2001-2016 Randal Henne, Greg Landrum and Rational Discovery LLC
  //
  //   @@ All Rights Reserved  @@
  //

#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <RDGeneral/RDLog.h>

#define YYDEBUG 1
#include "smiles.tab.hpp"

extern int yysmiles_lex(YYSTYPE *,void *,int &, unsigned int&);

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
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                RDKit::Atom* &,
                RDKit::Bond* &,
                unsigned int &,unsigned int &,
                std::vector<std::pair<unsigned int, unsigned int>>&,
		void *,int, unsigned int bad_token_position, const char * msg )
{
  yyErrorCleanup(ms);
  SmilesParseOps::detail::printSyntaxErrorMessage(input, msg, bad_token_position);
}

void
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::vector<std::pair<unsigned int, unsigned int>>&,
		void *,int, unsigned int bad_token_position, const char * msg )
{
  yyErrorCleanup(ms);
  SmilesParseOps::detail::printSyntaxErrorMessage(input, msg, bad_token_position);
}


%}

%define api.pure full
%lex-param   {yyscan_t *scanner}
%lex-param   {int& start_token}
%lex-param   {unsigned int& current_token_position}
%parse-param {const char *input}
%parse-param {std::vector<RDKit::RWMol *> *molList}
%parse-param {RDKit::Atom* &lastAtom}
%parse-param {RDKit::Bond* &lastBond}
%parse-param {unsigned &numAtomsParsed}
%parse-param {unsigned &numBondsParsed}
%parse-param {std::vector<std::pair<unsigned int, unsigned int>>& branchPoints}
%parse-param {void *scanner}
%parse-param {int& start_token}
%parse-param {unsigned int& current_token_position}

%code provides {
#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token, unsigned int& current_token_position)
}

%union {
  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  RDKit::Atom::ChiralType chiraltype;
  int                      ival;
}

%token START_MOL START_ATOM START_BOND;
%token <atom> AROMATIC_ATOM_TOKEN ATOM_TOKEN ORGANIC_ATOM_TOKEN
%token <ival> NONZERO_DIGIT_TOKEN ZERO_TOKEN
%token GROUP_OPEN_TOKEN GROUP_CLOSE_TOKEN SEPARATOR_TOKEN LOOP_CONNECTOR_TOKEN
%token MINUS_TOKEN PLUS_TOKEN
%token H_TOKEN AT_TOKEN PERCENT_TOKEN COLON_TOKEN HASH_TOKEN
%token <bond> BOND_TOKEN
%token <chiraltype> CHI_CLASS_TOKEN
%type <moli> mol
%type <atom> atomd element chiral_element h_element charge_element simple_atom
%type <bond> bondd
%type <ival> nonzero_number number ring_number digit branch_open_token
%token ATOM_OPEN_TOKEN ATOM_CLOSE_TOKEN
%token EOS_TOKEN

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
| START_BOND bondd EOS_TOKEN {
  lastBond = $2;
  YYACCEPT;
}
| START_BOND bondd {
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
| charge_element {
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
  RDKit::RWMol *curMol = (*molList)[ sz ];
  $1->setProp(RDKit::common_properties::_SmilesStart,1);
  curMol->addAtom($1, true, true);
  //delete $1;
  $$ = sz;
}

| mol atomd       {
  RWMol *mp = (*molList)[$$];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom($2,true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $2;
}

| mol BOND_TOKEN atomd  {
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
  mp->addBond($2,true);
  //delete $3;
}

| mol MINUS_TOKEN atomd {
  RWMol *mp = (*molList)[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($3,true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $3;
}

| mol SEPARATOR_TOKEN atomd {
  RWMol *mp = (*molList)[$$];
  $3->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom($3,true,true);
}

| mol ring_number {
  RWMol * mp = (*molList)[$$];
  Atom *atom=mp->getActiveAtom();
  mp->setAtomBookmark(atom,$2);

  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::UNSPECIFIED);
  mp->setBondBookmark(newB,$2);
  newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  if(!(mp->getAllBondsWithBookmark($2).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-($2+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
}

| mol BOND_TOKEN ring_number {
  RWMol * mp = (*molList)[$$];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     $2->getBondType());
  if($2->hasProp(RDKit::common_properties::_unspecifiedOrder)){
    newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  }
  newB->setBondDir($2->getBondDir());
  mp->setAtomBookmark(atom,$3);
  mp->setBondBookmark(newB,$3);
  if(!(mp->getAllBondsWithBookmark($3).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-($3+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
  delete $2;
}

| mol MINUS_TOKEN ring_number {
  RWMol * mp = (*molList)[$$];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::SINGLE);
  mp->setAtomBookmark(atom,$3);
  mp->setBondBookmark(newB,$3);
  if(!(mp->getAllBondsWithBookmark($3).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-($3+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
}

| mol branch_open_token atomd {
  RWMol *mp = (*molList)[$$];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom($3,true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  branchPoints.push_back({atomIdx1, $2});
}
| mol branch_open_token BOND_TOKEN atomd  {
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

  branchPoints.push_back({atomIdx1, $2});
}
| mol branch_open_token MINUS_TOKEN atomd {
  RWMol *mp = (*molList)[$$];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom($4,true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  branchPoints.push_back({atomIdx1, $2});
}
| mol GROUP_CLOSE_TOKEN {
  if(branchPoints.empty()){
     yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"extra close parentheses");
     yyErrorCleanup(molList);
     YYABORT;
  }
  RWMol *mp = (*molList)[$$];

  mp->setActiveAtom(branchPoints.back().first);
  branchPoints.pop_back();
}
;

// We'll use the token position for unclosed branch syntax error messages
branch_open_token: GROUP_OPEN_TOKEN { $$ = current_token_position; };

/* --------------------------------------------------------------- */
bondd:      BOND_TOKEN
          | MINUS_TOKEN {
          $$ = new Bond(Bond::SINGLE);
          }
;

/* --------------------------------------------------------------- */
atomd:	simple_atom
| ATOM_OPEN_TOKEN charge_element COLON_TOKEN number ATOM_CLOSE_TOKEN
{
  $$ = $2;
  $$->setNoImplicit(true);
  $$->setProp(RDKit::common_properties::molAtomMapNumber,$4);
}
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
                | number H_TOKEN { $$ = new Atom(1); $$->setIsotope($1); }
                | H_TOKEN H_TOKEN { $$ = new Atom(1); $$->setNumExplicitHs(1); }
                | number H_TOKEN H_TOKEN { $$ = new Atom(1); $$->setIsotope($1); $$->setNumExplicitHs(1);}
                | H_TOKEN H_TOKEN number { $$ = new Atom(1); $$->setNumExplicitHs($3); }
                | number H_TOKEN H_TOKEN number { $$ = new Atom(1); $$->setIsotope($1); $$->setNumExplicitHs($4);}
                | chiral_element
		| chiral_element H_TOKEN		{ $$ = $1; $1->setNumExplicitHs(1);}
		| chiral_element H_TOKEN number	{ $$ = $1; $1->setNumExplicitHs($3);}
		;

/* --------------------------------------------------------------- */
chiral_element:	 element
| element AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
| element AT_TOKEN AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
| element CHI_CLASS_TOKEN { $1->setChiralTag($2); $1->setProp(common_properties::_chiralPermutation,0); }
| element CHI_CLASS_TOKEN number { $1->setChiralTag($2); $1->setProp(common_properties::_chiralPermutation,$3); }
;

/* --------------------------------------------------------------- */
element:	simple_atom
		|	number simple_atom { $2->setIsotope( $1 ); $$ = $2; }
		|	ATOM_TOKEN
		|	number ATOM_TOKEN	   { $2->setIsotope( $1 ); $$ = $2; }
		|	HASH_TOKEN	number   { $$ = new Atom($2); }
		|	number HASH_TOKEN	number   { $$ = new Atom($3); $$->setIsotope($1); }
		;

/* --------------------------------------------------------------- */
simple_atom:       ORGANIC_ATOM_TOKEN
                | AROMATIC_ATOM_TOKEN
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
     yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"number too large");
     yyErrorCleanup(molList);
     YYABORT;
  }
  $$ = $1*10 + $2;
  }
;

digit: NONZERO_DIGIT_TOKEN
| ZERO_TOKEN
;

/*
  chival:	CHI_CLASS_TOKEN DIGIT_TOKEN
	| AT_TOKEN
        ;
*/


%%
