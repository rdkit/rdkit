%option reentrant
%option bison-bridge
%option noyywrap

%{

// $Id$
//
//  Copyright (C) 2001-2010 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <cstdio>
#ifdef WIN32
#include <io.h> 	 
#endif

#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/PeriodicTable.h>

#include <string>
#include <cstring>
#include "smiles.tab.hpp"

using namespace RDKit;

void setup_smiles_string(const std::string &text,yyscan_t yyscanner){
  YY_BUFFER_STATE buff=yysmiles__scan_string(text.c_str(),yyscanner);
  POSTCONDITION(buff,"invalid buffer");
}
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
<IN_ATOM_STATE>Lr 	{   yylval->atom = new Atom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
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
I			{	yylval->atom = new Atom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
				return ORGANIC_ATOM_TOKEN;
			}

H			{
				return H_TOKEN; 
			}

c		    {	yylval->atom = new Atom ( 6 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
n		    {	yylval->atom = new Atom( 7 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
o		    {	yylval->atom = new Atom( 8 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
p		    {	yylval->atom = new Atom( 15 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
s		    {	yylval->atom = new Atom( 16 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}

<IN_ATOM_STATE>se   {	yylval->atom = new Atom( 34 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}
<IN_ATOM_STATE>te   {	yylval->atom = new Atom( 52 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN; 
			}

\* 	            {   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
                                                        std::string("*"));
                                // must be ORGANIC_ATOM_TOKEN because
                                // we aren't in square brackets:
				return ORGANIC_ATOM_TOKEN; 
			}
<IN_ATOM_STATE>X 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
                                                        std::string("X"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xa 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
std::string("Xa"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xb 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
std::string("Xb"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xc 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
std::string("Xc"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xd 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
std::string("Xd"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xf 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
std::string("Xf"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xg 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
std::string("Xg"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xh 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp("dummyLabel",
std::string("Xh"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}
<IN_ATOM_STATE>Xi 	{   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp((const char *)"dummyLabel",
std::string("Xi"));
  BOOST_LOG(rdWarningLog)<<"Deprecation Warning: using X as a dummy-atom symbol is deprecated."<<std::endl;
				return ATOM_TOKEN; 
			}

\-			{ return MINUS_TOKEN; }

\+			{ return PLUS_TOKEN; }

[\=\#\:]    { yylval->bond = new Bond();
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
                yylval->bond->setIsAromatic(true);
		break;
              default:
                CHECK_INVARIANT(0,"cannot get here");
	      }
	      yylval->bond->setBondType(bt);
	return BOND_TOKEN; }	

[\\]    { yylval->bond = new Bond(Bond::SINGLE);
	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
	return BOND_TOKEN;  }
	
[\/]    { yylval->bond = new Bond(Bond::SINGLE);
	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
	return BOND_TOKEN;  }

\(       	{ return GROUP_OPEN_TOKEN; }
\)       	{ return GROUP_CLOSE_TOKEN; }


\[			{ BEGIN IN_ATOM_STATE; return ATOM_OPEN_TOKEN; }
<IN_ATOM_STATE>\]	{ BEGIN INITIAL; return ATOM_CLOSE_TOKEN; }

\.       	{ return SEPARATOR_TOKEN; }

\%              { return PERCENT_TOKEN; }

[0]		{ yylval->ival = 0; return ZERO_TOKEN; }
[1-9]		{ yylval->ival = atoi( yytext ); return NONZERO_DIGIT_TOKEN; }



\n		return 0;

<<EOF>>		{ return EOS_TOKEN; }
.		return yytext[0];

%%

#undef yysmiles_wrap
int yysmiles_wrap( void ) { return 1; }




