%option reentrant
%option bison-bridge
%option noyywrap

%{

// $Id$
//
//  Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <cstdio>
#ifdef WIN32
#include <io.h> 	 
#endif

#include <RDBoost/Exceptions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>

#include <string>
#include <cstring>
#include "smarts.tab.hpp"

using namespace RDKit;

//static PeriodicTable * gl_ptab = PeriodicTable::getTable();

void setup_smarts_string(const std::string &text,yyscan_t yyscanner){
  YY_BUFFER_STATE buff=yysmarts__scan_string(text.c_str(),yyscanner);
  POSTCONDITION(buff,"invalid buffer");
}

#define YY_FATAL_ERROR(msg) smarts_lexer_error(msg)

void smarts_lexer_error(const char *msg) {
     BOOST_LOG(rdErrorLog) << msg<<std::endl;
     throw ValueErrorException(msg);
}
%}
%option stack
%s IN_ATOM_STATE
%s IN_BRANCH_STATE
%s IN_RECURSION_STATE
%%

@[' ']*TH |
@[' ']*AL |
@[' ']*SQ |
@[' ']*BP |
@[' ']*OH 	{ return CHI_CLASS_TOKEN; }

@		{ return AT_TOKEN; }


<IN_ATOM_STATE>He |
<IN_ATOM_STATE>Li |
<IN_ATOM_STATE>Be |
<IN_ATOM_STATE>Ne |
<IN_ATOM_STATE>Na |
<IN_ATOM_STATE>Mg |
<IN_ATOM_STATE>Al |
<IN_ATOM_STATE>Si |
<IN_ATOM_STATE>Ar |
<IN_ATOM_STATE>K |
<IN_ATOM_STATE>Ca |
<IN_ATOM_STATE>Sc |
<IN_ATOM_STATE>Ti |
<IN_ATOM_STATE>V |
<IN_ATOM_STATE>Cr |
<IN_ATOM_STATE>Mn |
<IN_ATOM_STATE>Co |
<IN_ATOM_STATE>Fe |
<IN_ATOM_STATE>Ni |
<IN_ATOM_STATE>Cu |
<IN_ATOM_STATE>Zn |
<IN_ATOM_STATE>Ga |
<IN_ATOM_STATE>Ge |
<IN_ATOM_STATE>As |
<IN_ATOM_STATE>Se |
<IN_ATOM_STATE>Kr |
<IN_ATOM_STATE>Rb |
<IN_ATOM_STATE>Sr |
<IN_ATOM_STATE>Y |
<IN_ATOM_STATE>Zr |
<IN_ATOM_STATE>Nb |
<IN_ATOM_STATE>Mo |
<IN_ATOM_STATE>Tc |
<IN_ATOM_STATE>Ru |
<IN_ATOM_STATE>Rh |
<IN_ATOM_STATE>Pd |
<IN_ATOM_STATE>Ag |
<IN_ATOM_STATE>Cd |
<IN_ATOM_STATE>In |
<IN_ATOM_STATE>Sn |
<IN_ATOM_STATE>Sb |
<IN_ATOM_STATE>Te |
<IN_ATOM_STATE>Xe |
<IN_ATOM_STATE>Cs |
<IN_ATOM_STATE>Ba |
<IN_ATOM_STATE>La |
<IN_ATOM_STATE>Ce |
<IN_ATOM_STATE>Pr |
<IN_ATOM_STATE>Nd |
<IN_ATOM_STATE>Pm |
<IN_ATOM_STATE>Sm |
<IN_ATOM_STATE>Eu |
<IN_ATOM_STATE>Gd |
<IN_ATOM_STATE>Tb |
<IN_ATOM_STATE>Dy |
<IN_ATOM_STATE>Ho |
<IN_ATOM_STATE>Er |
<IN_ATOM_STATE>Tm |
<IN_ATOM_STATE>Yb |
<IN_ATOM_STATE>Lu |
<IN_ATOM_STATE>Hf |
<IN_ATOM_STATE>Ta |
<IN_ATOM_STATE>W |
<IN_ATOM_STATE>Re |
<IN_ATOM_STATE>Os |
<IN_ATOM_STATE>Ir |
<IN_ATOM_STATE>Pt |
<IN_ATOM_STATE>Au |
<IN_ATOM_STATE>Hg |
<IN_ATOM_STATE>Tl |
<IN_ATOM_STATE>Pb |
<IN_ATOM_STATE>Bi |
<IN_ATOM_STATE>Po |
<IN_ATOM_STATE>At |
<IN_ATOM_STATE>Rn |
<IN_ATOM_STATE>Fr |
<IN_ATOM_STATE>Ra |
<IN_ATOM_STATE>Ac |
<IN_ATOM_STATE>Th |
<IN_ATOM_STATE>Pa |
<IN_ATOM_STATE>U |
<IN_ATOM_STATE>Np |
<IN_ATOM_STATE>Pu |
<IN_ATOM_STATE>Am |
<IN_ATOM_STATE>Cm |
<IN_ATOM_STATE>Bk |
<IN_ATOM_STATE>Cf |
<IN_ATOM_STATE>Es |
<IN_ATOM_STATE>Fm |
<IN_ATOM_STATE>Md |
<IN_ATOM_STATE>No |
<IN_ATOM_STATE>Lr 	{   yylval->atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>D {
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomExplicitDegreeQuery(1));
	return COMPLEX_ATOM_QUERY_TOKEN;
}

<IN_ATOM_STATE>X {
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomTotalDegreeQuery(1));
	return COMPLEX_ATOM_QUERY_TOKEN;
}

<IN_ATOM_STATE>x {
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomHasRingBondQuery());
	return RINGBOND_ATOM_QUERY_TOKEN;
}

<IN_ATOM_STATE>v {
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomTotalValenceQuery(1));
	return COMPLEX_ATOM_QUERY_TOKEN;
}

<IN_ATOM_STATE>h {
	yylval->atom = new QueryAtom();
        yylval->atom->setQuery(makeAtomHasImplicitHQuery());
	return IMPLICIT_H_ATOM_QUERY_TOKEN;
}

<IN_ATOM_STATE>R {
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(new AtomRingQuery(-1));
	return COMPLEX_ATOM_QUERY_TOKEN;
}

<IN_ATOM_STATE>r {
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomInRingQuery());
	return RINGSIZE_ATOM_QUERY_TOKEN;
}

H			{  return H_TOKEN;  }


B			{  yylval->ival = 5;  return ORGANIC_ATOM_TOKEN;  }

C			{  yylval->ival = 6;  return ORGANIC_ATOM_TOKEN;  }

N			{  yylval->ival = 7;  return ORGANIC_ATOM_TOKEN;  }

O			{  yylval->ival = 8;  return ORGANIC_ATOM_TOKEN;  }

F			{  yylval->ival = 9;  return ORGANIC_ATOM_TOKEN;  }

P			{  yylval->ival = 15;  return ORGANIC_ATOM_TOKEN;  }

S			{  yylval->ival = 16;  return ORGANIC_ATOM_TOKEN;  }

Cl			{  yylval->ival = 17;  return ORGANIC_ATOM_TOKEN;  }

Br			{  yylval->ival = 35;  return ORGANIC_ATOM_TOKEN;  }

I			{  yylval->ival = 53;  return ORGANIC_ATOM_TOKEN;  }


b			{  yylval->ival = 5;  return AROMATIC_ATOM_TOKEN;  }

c			{  yylval->ival = 6;  return AROMATIC_ATOM_TOKEN;  }

n			{  yylval->ival = 7;  return AROMATIC_ATOM_TOKEN;  }

o			{  yylval->ival = 8;  return AROMATIC_ATOM_TOKEN;  }

p			{  yylval->ival = 15;  return AROMATIC_ATOM_TOKEN;  }

s			{  yylval->ival = 16;  return AROMATIC_ATOM_TOKEN;  }

<IN_ATOM_STATE>se	{  yylval->ival = 34;  return AROMATIC_ATOM_TOKEN;  }

<IN_ATOM_STATE>te	{  yylval->ival = 52;  return AROMATIC_ATOM_TOKEN;  }



\*			{
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomNullQuery());
	return SIMPLE_ATOM_QUERY_TOKEN;
}

a			{
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomAromaticQuery());
	yylval->atom->setIsAromatic(true);
	return SIMPLE_ATOM_QUERY_TOKEN;
}

A			{
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomAliphaticQuery());
	return SIMPLE_ATOM_QUERY_TOKEN;
}


\: 			{ return COLON_TOKEN; }

\_ 			{ return UNDERSCORE_TOKEN; }

\-			{ return MINUS_TOKEN; }

\+			{ return PLUS_TOKEN; }

\#			{ return HASH_TOKEN; }

\=	{ yylval->bond = new QueryBond(Bond::DOUBLE);
	yylval->bond->setQuery(makeBondOrderEqualsQuery(Bond::DOUBLE));
	return BOND_TOKEN;  }

\~	{ yylval->bond = new QueryBond();
	yylval->bond->setQuery(makeBondNullQuery());
	return BOND_TOKEN;  }

[\\]{1,2}    { yylval->bond = new QueryBond(Bond::SINGLE);
	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
	return BOND_TOKEN;  }
	
[\/]    { yylval->bond = new QueryBond(Bond::SINGLE);
	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
	return BOND_TOKEN;  }


<IN_ATOM_STATE>\$\(              { yy_push_state(IN_RECURSION_STATE,yyscanner); return BEGIN_RECURSE; }

\(       	{ yy_push_state(IN_BRANCH_STATE,yyscanner); return GROUP_OPEN_TOKEN; }
<IN_BRANCH_STATE>\)       	{ yy_pop_state(yyscanner); return GROUP_CLOSE_TOKEN; }
<IN_RECURSION_STATE>\)       	{ yy_pop_state(yyscanner); return END_RECURSE; }


\[			{ yy_push_state(IN_ATOM_STATE,yyscanner); return ATOM_OPEN_TOKEN; }
<IN_ATOM_STATE>\]	{ yy_pop_state(yyscanner); return ATOM_CLOSE_TOKEN; }
\]			{ /* FIX: ???
                           This rule is here because otherwise recursive SMARTS queries like:
	                   [$(C(=O)[O,N])] lex improperly (no ATOM_CLOSE token is returned).
 			   I am not 100% sure that the approach we're using here will work
                           all the time, but I'm hoping that any problems caused here in
                           the lexer will get caught in the parser. 
			  */	
                          return ATOM_CLOSE_TOKEN; }

\.       	{ return SEPARATOR_TOKEN; }

\%              { return PERCENT_TOKEN; }

[0]		{ yylval->ival = 0;  return ZERO_TOKEN; }
[1-9]		{ yylval->ival = yytext[0]-'0';  return NONZERO_DIGIT_TOKEN; }

\!			{ return NOT_TOKEN; }

\;			{ return SEMI_TOKEN; }

\&			{ return AND_TOKEN; }

\,			{ return OR_TOKEN; }

\^0		{
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::S));
	return HYB_TOKEN;
}

\^1		{
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP));
	return HYB_TOKEN;
}

\^2		{
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP2));
	return HYB_TOKEN;
}

\^3		{
	yylval->atom = new QueryAtom();
	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3));
	return HYB_TOKEN;
}

\n		return EOS_TOKEN;

<<EOF>>		{ return EOS_TOKEN; }
.		return yytext[0];

%%

#undef yysmarts_wrap
int yysmarts_wrap( void ) { return 1; }




