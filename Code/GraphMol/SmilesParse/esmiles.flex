%{
// $Id: esmiles.flex 4974 2006-02-18 00:49:21Z glandrum $
//
//  Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#pragma warning(disable:4786)


#include <stdio.h>

#undef yywrap 
#undef YY_INPUT

#include "InputFiller.h"

extern INPUT_FUNC_TYPE gp_myInput;

#define YY_INPUT( b, r, ms) (r = gp_myInput( b, ms ))

#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/PeriodicTable.h>

#ifdef WIN32
#include <io.h>
#endif
#include <string>
#include "smiles.tab.h"

#ifdef WIN32
int
isatty( int fd )
{ 
	return _isatty( fd );
}
#endif
using namespace RDKit;

static PeriodicTable * gl_ptab = PeriodicTable::getTable();

%}
%s IN_ATOM_STATE
%%

@[' ']*TH |
@[' ']*AL |
@[' ']*SQ |
@[' ']*BP |
@[' ']*OH 	{ return CHI_CLASS; }

@+		{ return AT; }


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
<IN_ATOM_STATE>Lr 	{   yylval.atom = new Atom( gl_ptab->getAtomicNumber( yytext ) );
				return ATOM; 
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
I			{	yylval.atom = new Atom( gl_ptab->getAtomicNumber( yytext ) );
				return ORGANIC_ATOM;
			}

H			{
				return H; 
			}

c		    {	yylval.atom = new Atom ( 6 );
			yylval.atom->setIsAromatic(true);
				return AROMATIC_ATOM; 
			}
n		    {	yylval.atom = new Atom( 7 );
			yylval.atom->setIsAromatic(true);
				return AROMATIC_ATOM; 
			}
o		    {	yylval.atom = new Atom( 8 );
			yylval.atom->setIsAromatic(true);
				return AROMATIC_ATOM; 
			}
p		    {	yylval.atom = new Atom( 15 );
			yylval.atom->setIsAromatic(true);
				return AROMATIC_ATOM; 
			}
s		    {	yylval.atom = new Atom( 16 );
			yylval.atom->setIsAromatic(true);
				return AROMATIC_ATOM; 
			}

\-			{ return MINUS; }

\+			{ return PLUS; }

\<			{ return LT; }

\>			{ return GT; }

[\=\#\:]    { yylval.bond = new Bond();
              Bond::BondType bt;
              switch(yytext[0]){
	      case '=':
		bt = Bond::DOUBLE;
		break;
	      case '#':
		bt = Bond::TRIPLE;
		break;
	      case ':':
		bt = Bond::AROMATIC;
		break;
	      }
	      yylval.bond->setBondType(bt);
	return BOND; }

[\\]    { yylval.bond = new Bond(Bond::SINGLE);
	yylval.bond->setBondDir(Bond::ENDDOWNRIGHT);
	return BOND;  }
	
[\/]    { yylval.bond = new Bond(Bond::SINGLE);
	yylval.bond->setBondDir(Bond::ENDUPRIGHT);
	return BOND;  }

\(       	{ return GROUP_OPEN; }
\)       	{ return GROUP_CLOSE; }


\[			{ BEGIN IN_ATOM_STATE; return ATOM_OPEN; }
<IN_ATOM_STATE>\]	{ BEGIN INITIAL; return ATOM_CLOSE; }

\.       	{ return SEPARATOR; }

\_	        { return DATIVE_MARK; }

\%              { return PERCENT; }

[0-9]		{ yylval.ival = atoi( yytext ); return DIGIT; }



\n		return 0;

.		return yytext[0];

%%

int yywrap( void ) { return 1; }


