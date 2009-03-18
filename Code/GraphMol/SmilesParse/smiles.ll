%{

// $Id$
//
//  Copyright (C) 2001-2008 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <cstdio>
#ifdef WIN32
#include <io.h> 	 
#endif

#undef YY_INPUT

#include <GraphMol/SmilesParse/InputFiller.h>

extern INPUT_FUNC_TYPE gp_myInput;

#define YY_INPUT( b, r, ms) (r = gp_myInput( b, ms ))

#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/PeriodicTable.h>

#include <string>
#include <cstring>
#include "smiles.tab.hpp"

using namespace RDKit;

//static PeriodicTable * gl_ptab = PeriodicTable::getTable();

%}
%s IN_ATOM_STATE
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
<IN_ATOM_STATE>Lr 	{   yysmiles_lval.atom = new Atom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
				return ATOM_TOKEN; 
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
I			{	yysmiles_lval.atom = new Atom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
				return ORGANIC_ATOM_TOKEN;
			}

H			{
				return H_TOKEN; 
			}

c		    {	yysmiles_lval.atom = new Atom ( 6 );
			yysmiles_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
n		    {	yysmiles_lval.atom = new Atom( 7 );
			yysmiles_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
o		    {	yysmiles_lval.atom = new Atom( 8 );
			yysmiles_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
p		    {	yysmiles_lval.atom = new Atom( 15 );
			yysmiles_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
s		    {	yysmiles_lval.atom = new Atom( 16 );
			yysmiles_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}

<IN_ATOM_STATE>se   {	yysmiles_lval.atom = new Atom( 34 );
			yysmiles_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
<IN_ATOM_STATE>te   {	yysmiles_lval.atom = new Atom( 52 );
			yysmiles_lval.atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}

\* 	            {   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
                                                        std::string("*"));
                                // must be ORGANIC_ATOM_TOKEN because
                                // we aren't in square brackets:
				return ORGANIC_ATOM_TOKEN; 
			}
<IN_ATOM_STATE>X 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
                                                        std::string("X"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xa 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
std::string("Xa"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xb 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
std::string("Xb"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xc 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
std::string("Xc"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xd 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
std::string("Xd"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xf 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
std::string("Xf"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xg 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
std::string("Xg"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xh 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp("dummyLabel",
std::string("Xh"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xi 	{   yysmiles_lval.atom = new Atom( 0 );
		            yysmiles_lval.atom->setProp((const char *)"dummyLabel",
std::string("Xi"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}

\-			{ return MINUS_TOKEN; }

\+			{ return PLUS_TOKEN; }

[\=\#\:]    { yysmiles_lval.bond = new Bond();
              Bond::BondType bt=Bond::UNSPECIFIED;
              switch(yytext[0]){
	      case '=':
		bt = Bond::DOUBLE;
		break;
	      case '#':
		bt = Bond::TRIPLE;
		break;
	      case ':':
		bt = Bond::AROMATIC;
                yysmiles_lval.bond->setIsAromatic(true);
		break;
              default:
                CHECK_INVARIANT(0,"cannot get here");
	      }
	      yysmiles_lval.bond->setBondType(bt);
	return BOND_TOKEN; }	

[\\]    { yysmiles_lval.bond = new Bond(Bond::SINGLE);
	yysmiles_lval.bond->setBondDir(Bond::ENDDOWNRIGHT);
	return BOND_TOKEN;  }
	
[\/]    { yysmiles_lval.bond = new Bond(Bond::SINGLE);
	yysmiles_lval.bond->setBondDir(Bond::ENDUPRIGHT);
	return BOND_TOKEN;  }

\(       	{ return GROUP_OPEN_TOKEN; }
\)       	{ return GROUP_CLOSE_TOKEN; }


\[			{ BEGIN IN_ATOM_STATE; return ATOM_OPEN_TOKEN; }
<IN_ATOM_STATE>\]	{ BEGIN INITIAL; return ATOM_CLOSE_TOKEN; }

\.       	{ return SEPARATOR_TOKEN; }

\%              { return PERCENT_TOKEN; }

[0]		{ yysmiles_lval.ival = 0; return ZERO_TOKEN; }
[1-9]		{ yysmiles_lval.ival = atoi( yytext ); return NONZERO_DIGIT_TOKEN; }



\n		return 0;

<<EOF>>		{ return EOS_TOKEN; }
.		return yytext[0];

%%

#undef yysmiles_wrap
int yysmiles_wrap( void ) { return 1; }




