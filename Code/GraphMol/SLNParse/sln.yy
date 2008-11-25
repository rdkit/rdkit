%{

  // $Id: sln.yy 1043 2008-07-25 12:00:41Z landrgr1 $
  //
  //  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
  //  All rights reserved.
  // 
  // Redistribution and use in source and binary forms, with or without
  // modification, are permitted provided that the following conditions are
  // met: 
  //
  //     * Redistributions of source code must retain the above copyright 
  //       notice, this list of conditions and the following disclaimer.
  //     * Redistributions in binary form must reproduce the above
  //       copyright notice, this list of conditions and the following 
  //       disclaimer in the documentation and/or other materials provided 
  //       with the distribution.
  //     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
  //       nor the names of its contributors may be used to endorse or promote 
  //       products derived from this software without specific prior
  //       written permission.
  //
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  //
  // Created by Greg Landrum, September 2006
  //

#include <cstring>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SLNParse/SLNParseOps.h>
#include <GraphMol/SLNParse/SLNAttribs.h>
#include <GraphMol/SLNParse/SLNParse.h>
#include <RDGeneral/RDLog.h>

extern int yysln_lex();

#define YYDEBUG 1

void
yysln_error( const char * msg )
{

}

 extern bool slnParserDoQueries;
 namespace RDKit{
   namespace SLNParse{
     extern std::vector<RDKit::RWMol *> molList_g;
   }
 }
 namespace SLNParse = RDKit::SLNParse;

%}
 
 
%union {
  int                      mol_T;
  RDKit::Atom *            atom_T;
  RDKit::Bond *            bond_T;
  int                      ival_T;
  std::string*             text_T;
  char                     char_T;
  RDKit::SLNParse::AttribType       *attrib_T;
  RDKit::SLNParse::AttribListType   *attriblist_T;
}



%type <atom_T> atom primatom
%type <bond_T> bond primbond onebond
%type <attriblist_T> attriblist
%type <attrib_T> attrib recursivequery
%type <mol_T> mol cmpd primmol
%token <text_T> TEXT_BLOCK
%token<char_T> CHAR_TOKEN
%token <ival_T> NUMBER_TOKEN
%token H_TOKEN AT_TOKEN
%token <atom_T> ATOM_TOKEN
%token <text_T> COMPARE_TOKEN
%token OPEN_PAREN_TOKEN CLOSE_PAREN_TOKEN
%token OPEN_BRACKET_TOKEN CLOSE_BRACKET_TOKEN
%token OPEN_ANGLE_TOKEN CLOSE_ANGLE_TOKEN
%token SEPARATOR_TOKEN
%token ASTERIX_TOKEN
%token EOS_TOKEN
%token PLUS_TOKEN MINUS_TOKEN
%token COLON_TOKEN EQUALS_TOKEN TILDE_TOKEN HASH_TOKEN COMMA_TOKEN
%token NOT_TOKEN AND_TOKEN OR_TOKEN SEMI_TOKEN
%token RECURSE_TOKEN NEG_RECURSE_TOKEN
%token ERROR_TOKEN

%left SEMI_TOKEN
%left OR_TOKEN
%left AND_TOKEN
%right NOT_TOKEN

%%

/* --------------------------------------------------------------- */
cmpd: mol
| cmpd SEPARATOR_TOKEN mol {
  $$=SLNParse::addFragToMol(SLNParse::molList_g,$1,$3);
}
| cmpd error EOS_TOKEN {
  yyclearin;
  yyerrok;
  YYABORT;
}
| cmpd EOS_TOKEN {
  YYACCEPT;
}
| error EOS_TOKEN {
  yyclearin;
  yyerrok;
  YYABORT;
}
;

mol: primmol
| primmol OPEN_ANGLE_TOKEN attriblist CLOSE_ANGLE_TOKEN
;


primmol: atom {
  $$=SLNParse::startMol(SLNParse::molList_g,$1);
}
| primmol atom {
  SLNParse::addAtomToMol(SLNParse::molList_g,$$,$2);
  $$=$1;
}
| primmol bond atom{
  SLNParse::addAtomToMol(SLNParse::molList_g,$$,$3,$2);
  $$=$1;
}
| primmol AT_TOKEN NUMBER_TOKEN {
  SLNParse::closeRingBond(SLNParse::molList_g,$$,$3);
  $$=$1;
}
| primmol bond AT_TOKEN NUMBER_TOKEN {
  // closeRingBond() takes ownership of the bond
  SLNParse::closeRingBond(SLNParse::molList_g,$$,$4,$2);
  $$=$1;
}
| primmol OPEN_PAREN_TOKEN primmol CLOSE_PAREN_TOKEN {
  SLNParse::addBranchToMol(SLNParse::molList_g,$$,$3);
  $$=$1;
}
| primmol OPEN_PAREN_TOKEN bond primmol CLOSE_PAREN_TOKEN{
  // addBranchToMol() takes ownership of the bond and deletes the
  // branch, so no leaks here'
  SLNParse::addBranchToMol(SLNParse::molList_g,$$,$4,$3);
  $$=$1;
}
| primmol OPEN_PAREN_TOKEN bond AT_TOKEN NUMBER_TOKEN CLOSE_PAREN_TOKEN{
  SLNParse::closeRingBond(SLNParse::molList_g,$$,$5,$3);
  $$=$1;
}
;

atom: primatom
| H_TOKEN {
  $$ = new RDKit::Atom(1);
}
| primatom H_TOKEN {
  $1->setNumExplicitHs(1);
  $$=$1;
}
| primatom H_TOKEN NUMBER_TOKEN {
  $1->setNumExplicitHs($3);
  $$=$1;
}
;

primatom: ATOM_TOKEN 
| primatom ASTERIX_TOKEN{
  $$->setProp("_starred",1,true);
}
| primatom OPEN_BRACKET_TOKEN NUMBER_TOKEN CLOSE_BRACKET_TOKEN {
  $1->setProp("_AtomID",static_cast<unsigned int>($3));
  $$=$1;
}
| primatom OPEN_BRACKET_TOKEN NUMBER_TOKEN COLON_TOKEN attriblist CLOSE_BRACKET_TOKEN {
  $1->setProp("_AtomID",static_cast<unsigned int>($3));
  SLNParse::parseAtomAttribs($1,*$5,slnParserDoQueries);
  delete $5;
  $$=$1;
}
| primatom OPEN_BRACKET_TOKEN attriblist CLOSE_BRACKET_TOKEN {
  SLNParse::parseAtomAttribs($1,*$3,slnParserDoQueries);
  delete $3;
  $$=$1;
}
;


bond: primbond
| primbond OPEN_BRACKET_TOKEN attriblist CLOSE_BRACKET_TOKEN {
  SLNParse::parseBondAttribs($1,*$3,slnParserDoQueries);
  delete $3;
  $$ = $1;
}

/* tildes can't be mixed with regular bonds in expressions: */
| TILDE_TOKEN { 
  RDKit::Bond *bond=new RDKit::QueryBond();
  bond->setQuery(RDKit::makeBondNullQuery());   
  $$ = bond;
}
| TILDE_TOKEN OPEN_BRACKET_TOKEN attriblist CLOSE_BRACKET_TOKEN {
  RDKit::Bond *bond=new RDKit::QueryBond();
  bond->setQuery(RDKit::makeBondNullQuery());   
  SLNParse::parseBondAttribs(bond,*$3,slnParserDoQueries);
  delete $3;
  $$ = bond;
}
;

primbond: onebond
| primbond onebond {
	if(!slnParserDoQueries){
    BOOST_LOG(rdErrorLog) << "SLN Parse Error: sequential bonds not allowed in non-queries" << std::endl;
    SLNParse::molList_g.clear();
    SLNParse::molList_g.resize(0);
    YYABORT;
	} else {
	  RDKit::QueryBond *b1=static_cast<RDKit::QueryBond *>($1);
	  RDKit::QueryBond *b2=static_cast<RDKit::QueryBond *>($2);
	  b1->expandQuery(b2->getQuery()->copy(),Queries::COMPOSITE_OR,true);
		delete b2;
	}
}
;

onebond: MINUS_TOKEN {
  RDKit::Bond *bond;
  if(slnParserDoQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::SINGLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::SINGLE);
  }
  $$ = bond;
}
| EQUALS_TOKEN {
  RDKit::Bond *bond;
  if(slnParserDoQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::DOUBLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::DOUBLE);
  }
  $$ = bond;
}
| HASH_TOKEN {
  RDKit::Bond *bond;
  if(slnParserDoQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::TRIPLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::TRIPLE);
  }
  $$ = bond;

}
| COLON_TOKEN {
  RDKit::Bond *bond;
  if(slnParserDoQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::AROMATIC);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::AROMATIC);
  }
  $$ = bond;
}
;


attriblist: attrib {
  $$ = new SLNParse::AttribListType();
  $$->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>($1)));
}
| attriblist SEMI_TOKEN attrib{
  $$->push_back(std::make_pair(SLNParse::AttribLowPriAnd,
                               boost::shared_ptr<SLNParse::AttribType>($3)));
}
| attriblist AND_TOKEN attrib{
  $$->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>($3)));
}
| attriblist OR_TOKEN attrib{
  $$->push_back(std::make_pair(SLNParse::AttribOr,
                               boost::shared_ptr<SLNParse::AttribType>($3)));
}
;

attrib: TEXT_BLOCK {
  $$ = new SLNParse::AttribType();
  $$->first = *$1;
  boost::to_lower($$->first);
  $$->op = "";
  $$->second = "";
  delete $1;
}
| NOT_TOKEN attrib {
  $2->negated=true;
  $$=$2;
}
| TEXT_BLOCK COMPARE_TOKEN TEXT_BLOCK {
  $$ = new SLNParse::AttribType();
  $$->first = *$1;
  $$->op = *$2;
  $$->second = *$3;
  boost::to_lower($$->first);
  boost::to_lower($$->second);
  delete $1;
  delete $2;
  delete $3;
}
| PLUS_TOKEN {
  $$ = new SLNParse::AttribType();
  $$->first = "charge";
  $$->op = "=";
  $$->second = "+1";
}
| PLUS_TOKEN NUMBER_TOKEN {
  $$ = new SLNParse::AttribType();
  $$->first = "charge";
  $$->op = "=";
  $$->second = SLNParse::convertToString($2);
}
| MINUS_TOKEN {
  $$ = new SLNParse::AttribType();
  $$->first = "charge";
  $$->op = "=";
  $$->second = "-1";
}
| MINUS_TOKEN NUMBER_TOKEN {
  $$ = new SLNParse::AttribType();
  $$->first = "charge";
  $$->op = "=";
  $$->second = SLNParse::convertToString(-$2);
}
| recursivequery {
  $$ = $1;
}
;

recursivequery: RECURSE_TOKEN cmpd {
   int sz = SLNParse::molList_g.size();
   RDKit::ROMol *mol=SLNParse::molList_g[$2];
   SLNParse::molList_g.resize( sz-1 );
   RDKit::RWMol *tmp=SLNParse::finalizeQueryMol(mol,true);
   delete mol;
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(tmp);
   RDKit::ATOM_OR_QUERY *orq=new RDKit::ATOM_OR_QUERY();
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   $$ = new SLNParse::AttribType();
   $$->first="is";
   $$->op = "=";
   $$->second = "";
   $$->structQuery=static_cast<void *>(orq);
}
| NEG_RECURSE_TOKEN cmpd {
   int sz = SLNParse::molList_g.size();
   RDKit::ROMol *mol=SLNParse::molList_g[$2];
   SLNParse::molList_g.resize( sz-1 );
   RDKit::RWMol *tmp=SLNParse::finalizeQueryMol(mol,true);
   delete mol;
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(tmp);
   RDKit::ATOM_OR_QUERY *orq=new RDKit::ATOM_OR_QUERY();
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   orq->setNegation(true);

   $$ = new SLNParse::AttribType();
   $$->first="is";
   $$->op = "=";
   $$->second = "";
   $$->structQuery=static_cast<void *>(orq);
} 
| recursivequery COMMA_TOKEN cmpd {
   int sz = SLNParse::molList_g.size();
   RDKit::ROMol *mol=SLNParse::molList_g[$3];
   SLNParse::molList_g.resize( sz-1 );
   RDKit::RWMol *tmp=SLNParse::finalizeQueryMol(mol,true);
   delete mol;
   
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(tmp);

   RDKit::ATOM_OR_QUERY *orq=static_cast<RDKit::ATOM_OR_QUERY *>($1->structQuery);
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   $$=$1;
}
;
%%


