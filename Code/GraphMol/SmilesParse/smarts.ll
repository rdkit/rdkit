%{

// $Id$
//
//  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <cstdio>

#undef YY_INPUT

#include <GraphMol/SmilesParse/InputFiller.h>

extern INPUT_FUNC_TYPE gp_myInput;

#define YY_INPUT( b, r, ms) (r = gp_myInput( b, ms ))

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>

#include <string>
#include <cstring>
#include "smarts.tab.hpp"

using namespace RDKit;

//static PeriodicTable * gl_ptab = PeriodicTable::getTable();

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
<IN_ATOM_STATE>Lr 	{   yysmarts_lval.atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>D {
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomExplicitDegreeQuery(0));
	return COMPLEX_ATOM_QUERY_TOKEN;
}
<IN_ATOM_STATE>X {
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomTotalDegreeQuery(0));
	return COMPLEX_ATOM_QUERY_TOKEN;
}
<IN_ATOM_STATE>v {
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomTotalValenceQuery(0));
	return COMPLEX_ATOM_QUERY_TOKEN;
}
<IN_ATOM_STATE>h {
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomImplicitValenceQuery(1));
	return COMPLEX_ATOM_QUERY_TOKEN;
}
<IN_ATOM_STATE>R {
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(new AtomRingQuery(-1));
	return COMPLEX_ATOM_QUERY_TOKEN;
}
<IN_ATOM_STATE>r {
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomInRingQuery());
	return RINGSIZE_ATOM_QUERY_TOKEN;
}

H			{
				return H_TOKEN; 
			}


B  |
C  |
N  |
O  |
P  |
S  |
F  |
Cl |
Br | 
I			{	yysmarts_lval.atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
				return ORGANIC_ATOM_TOKEN;
			}

c		    {	yysmarts_lval.atom = new QueryAtom ( 6 );
			yysmarts_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
n		    {	yysmarts_lval.atom = new QueryAtom( 7 );
			yysmarts_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
o		    {	yysmarts_lval.atom = new QueryAtom( 8 );
			yysmarts_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
p		    {	yysmarts_lval.atom = new QueryAtom( 15 );
			yysmarts_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
s		    {	yysmarts_lval.atom = new QueryAtom( 16 );
			yysmarts_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
<IN_ATOM_STATE>se   {	yysmarts_lval.atom = new QueryAtom( 34 );
			yysmarts_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
<IN_ATOM_STATE>te   {	yysmarts_lval.atom = new QueryAtom( 52 );
			yysmarts_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}



\*			{
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomNullQuery());
	return SIMPLE_ATOM_QUERY_TOKEN;
}

a			{
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomAromaticQuery());
	yysmarts_lval.atom->setIsAromatic(true);
	return SIMPLE_ATOM_QUERY_TOKEN;
}
A			{
	yysmarts_lval.atom = new QueryAtom();
	yysmarts_lval.atom->setQuery(makeAtomAliphaticQuery());
	return SIMPLE_ATOM_QUERY_TOKEN;
}


\: 			{ return COLON_TOKEN; }

\-			{ return MINUS_TOKEN; }

\+			{ return PLUS_TOKEN; }

\#			{ return HASH_TOKEN; }

[\=\~]    { yysmarts_lval.bond = new QueryBond();
              switch(yytext[0]){
	      case '=':
		yysmarts_lval.bond->setBondType(Bond::DOUBLE);
                yysmarts_lval.bond->setQuery(makeBondOrderEqualsQuery(Bond::DOUBLE));
		break;
	      case '~':
		yysmarts_lval.bond->setQuery(makeBondNullQuery());
		break;
	      }
	return BOND_TOKEN; }

[\\]    { yysmarts_lval.bond = new QueryBond(Bond::SINGLE);
	yysmarts_lval.bond->setBondDir(Bond::ENDDOWNRIGHT);
        yysmarts_lval.bond->setQuery(makeBondDirEqualsQuery(Bond::ENDDOWNRIGHT));
	return BOND_TOKEN;  }
	
[\/]    { yysmarts_lval.bond = new QueryBond(Bond::SINGLE);
	yysmarts_lval.bond->setBondDir(Bond::ENDUPRIGHT);
        yysmarts_lval.bond->setQuery(makeBondDirEqualsQuery(Bond::ENDUPRIGHT));
	return BOND_TOKEN;  }


<IN_ATOM_STATE>\$\(              { yy_push_state(IN_RECURSION_STATE); return BEGIN_RECURSE; }

\(       	{ yy_push_state(IN_BRANCH_STATE); return GROUP_OPEN_TOKEN; }
<IN_BRANCH_STATE>\)       	{ yy_pop_state(); return GROUP_CLOSE_TOKEN; }
<IN_RECURSION_STATE>\)       	{ yy_pop_state(); return END_RECURSE; }


\[			{ yy_push_state(IN_ATOM_STATE); return ATOM_OPEN_TOKEN; }
<IN_ATOM_STATE>\]	{ yy_pop_state(); return ATOM_CLOSE_TOKEN; }
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

[0-9]		{ yysmarts_lval.ival = atoi( yytext ); return DIGIT_TOKEN; }

\!			{ return NOT_TOKEN; }

\;			{ return SEMI_TOKEN; }

\&			{ return AND_TOKEN; }

\,			{ return OR_TOKEN; }

\^			{ return HYB_TOKEN; }

\n		return EOS_TOKEN;

<<EOF>>		{ return EOS_TOKEN; }
.		return yytext[0];

%%

#undef yysmarts_wrap
int yysmarts_wrap( void ) { return 1; }




