%{

// $Id$
//
//  Copyright (C) 2001-2010 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include <cstdio>

#include <RDGeneral/Exceptions.h>
#include <RDGeneral/types.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RDKitQueries.h>

#include <string>
#include <cstring>
#include "smiles.tab.hpp"
#include "SmilesScanner.h"

using namespace RDKit;

#undef YY_FATAL_ERROR
#define YY_FATAL_ERROR(msg) smiles_lexer_error(msg)

void smiles_lexer_error(const char *msg) {
     BOOST_LOG(rdErrorLog) << msg<<std::endl;
     throw ValueErrorException(msg);
}

#undef YY_DECL
#define YY_DECL int SmilesScanner::lex(SmilesParser::semantic_type* const lval, SmilesParser::location_type* location, int& start_token)

// we're doing this to track the position of the current token.
#undef YY_USER_ACTION
#define YY_USER_ACTION location->begin += yyleng; location->columns(yyleng);

// using an alias since the C++ parser uses lval
#undef yylval
#define yylval lval

using token = SmilesParser::token_kind_type;

%}

%option c++
%option full
%option never-interactive
%option noyywrap
%option yylineno
%option yyclass="RDKit::SmilesScanner"
%option nodefault
%option nostdinit


%s IN_ATOM_STATE
%%

%{
  if (start_token)
    {
      int t = start_token;
      start_token = 0;
      return t;
    }
%}

@[' ']*TH { yylval->chiraltype = Atom::ChiralType::CHI_TETRAHEDRAL; return token::CHI_CLASS_TOKEN; }
@[' ']*AL { yylval->chiraltype = Atom::ChiralType::CHI_ALLENE; return token::CHI_CLASS_TOKEN; }
@[' ']*SP { yylval->chiraltype = Atom::ChiralType::CHI_SQUAREPLANAR; return token::CHI_CLASS_TOKEN; }
@[' ']*TB { yylval->chiraltype = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL; return token::CHI_CLASS_TOKEN; }
@[' ']*OH { yylval->chiraltype = Atom::ChiralType::CHI_OCTAHEDRAL; return token::CHI_CLASS_TOKEN; }

@		{ return token::AT_TOKEN; }


<IN_ATOM_STATE>He	{ yylval->atom = new Atom(2); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Li	{ yylval->atom = new Atom(3); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Be	{ yylval->atom = new Atom(4); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ne	{ yylval->atom = new Atom(10); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Na	{ yylval->atom = new Atom(11); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Mg	{ yylval->atom = new Atom(12); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Al	{ yylval->atom = new Atom(13); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Si	{ yylval->atom = new Atom(14); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ar	{ yylval->atom = new Atom(18); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>K	{ yylval->atom = new Atom(19); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ca	{ yylval->atom = new Atom(20); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Sc	{ yylval->atom = new Atom(21); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ti	{ yylval->atom = new Atom(22); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>V	{ yylval->atom = new Atom(23); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Cr	{ yylval->atom = new Atom(24); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Mn	{ yylval->atom = new Atom(25); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Fe	{ yylval->atom = new Atom(26); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Co	{ yylval->atom = new Atom(27); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ni	{ yylval->atom = new Atom(28); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Cu	{ yylval->atom = new Atom(29); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Zn	{ yylval->atom = new Atom(30); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ga	{ yylval->atom = new Atom(31); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ge	{ yylval->atom = new Atom(32); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>As	{ yylval->atom = new Atom(33); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Se	{ yylval->atom = new Atom(34); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Kr	{ yylval->atom = new Atom(36); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Rb	{ yylval->atom = new Atom(37); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Sr	{ yylval->atom = new Atom(38); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Y	{ yylval->atom = new Atom(39); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Zr	{ yylval->atom = new Atom(40); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Nb	{ yylval->atom = new Atom(41); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Mo	{ yylval->atom = new Atom(42); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Tc	{ yylval->atom = new Atom(43); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ru	{ yylval->atom = new Atom(44); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Rh	{ yylval->atom = new Atom(45); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Pd	{ yylval->atom = new Atom(46); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ag	{ yylval->atom = new Atom(47); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Cd	{ yylval->atom = new Atom(48); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>In	{ yylval->atom = new Atom(49); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Sn	{ yylval->atom = new Atom(50); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Sb	{ yylval->atom = new Atom(51); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Te	{ yylval->atom = new Atom(52); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Xe	{ yylval->atom = new Atom(54); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Cs	{ yylval->atom = new Atom(55); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ba	{ yylval->atom = new Atom(56); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>La	{ yylval->atom = new Atom(57); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ce	{ yylval->atom = new Atom(58); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Pr	{ yylval->atom = new Atom(59); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Nd	{ yylval->atom = new Atom(60); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Pm	{ yylval->atom = new Atom(61); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Sm	{ yylval->atom = new Atom(62); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Eu	{ yylval->atom = new Atom(63); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Gd	{ yylval->atom = new Atom(64); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Tb	{ yylval->atom = new Atom(65); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Dy	{ yylval->atom = new Atom(66); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ho	{ yylval->atom = new Atom(67); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Er	{ yylval->atom = new Atom(68); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Tm	{ yylval->atom = new Atom(69); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Yb	{ yylval->atom = new Atom(70); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Lu	{ yylval->atom = new Atom(71); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Hf	{ yylval->atom = new Atom(72); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ta	{ yylval->atom = new Atom(73); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>W	{ yylval->atom = new Atom(74); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Re	{ yylval->atom = new Atom(75); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Os	{ yylval->atom = new Atom(76); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ir	{ yylval->atom = new Atom(77); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Pt	{ yylval->atom = new Atom(78); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Au	{ yylval->atom = new Atom(79); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Hg	{ yylval->atom = new Atom(80); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Tl	{ yylval->atom = new Atom(81); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Pb	{ yylval->atom = new Atom(82); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Bi	{ yylval->atom = new Atom(83); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Po	{ yylval->atom = new Atom(84); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>At	{ yylval->atom = new Atom(85); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Rn	{ yylval->atom = new Atom(86); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Fr	{ yylval->atom = new Atom(87); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ra	{ yylval->atom = new Atom(88); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ac	{ yylval->atom = new Atom(89); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Th	{ yylval->atom = new Atom(90); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Pa	{ yylval->atom = new Atom(91); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>U	{ yylval->atom = new Atom(92); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Np	{ yylval->atom = new Atom(93); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Pu	{ yylval->atom = new Atom(94); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Am	{ yylval->atom = new Atom(95); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Cm	{ yylval->atom = new Atom(96); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Bk	{ yylval->atom = new Atom(97); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Cf	{ yylval->atom = new Atom(98); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Es	{ yylval->atom = new Atom(99); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Fm	{ yylval->atom = new Atom(100); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Md	{ yylval->atom = new Atom(101); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>No	{ yylval->atom = new Atom(102); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Lr	{ yylval->atom = new Atom(103); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Rf	{ yylval->atom = new Atom(104); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Db	{ yylval->atom = new Atom(105); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Sg	{ yylval->atom = new Atom(106); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Bh	{ yylval->atom = new Atom(107); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Hs	{ yylval->atom = new Atom(108); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Mt	{ yylval->atom = new Atom(109); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ds	{ yylval->atom = new Atom(110); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Rg	{ yylval->atom = new Atom(111); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Cn	{ yylval->atom = new Atom(112); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Nh	{ yylval->atom = new Atom(113); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Fl	{ yylval->atom = new Atom(114); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Mc	{ yylval->atom = new Atom(115); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Lv	{ yylval->atom = new Atom(116); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Ts	{ yylval->atom = new Atom(117); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Og	{ yylval->atom = new Atom(118); return token::ATOM_TOKEN; }

<IN_ATOM_STATE>Uun	{ yylval->atom = new Atom(110); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uuu	{ yylval->atom = new Atom(111); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uub	{ yylval->atom = new Atom(112); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uut	{ yylval->atom = new Atom(113); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uuq	{ yylval->atom = new Atom(114); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uup	{ yylval->atom = new Atom(115); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uuh	{ yylval->atom = new Atom(116); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uus	{ yylval->atom = new Atom(117); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>Uuo	{ yylval->atom = new Atom(118); return token::ATOM_TOKEN; }

B  { yylval->atom = new Atom(5);return token::ORGANIC_ATOM_TOKEN; }
C  { yylval->atom = new Atom(6);return token::ORGANIC_ATOM_TOKEN; }
N  { yylval->atom = new Atom(7);return token::ORGANIC_ATOM_TOKEN; }
O  { yylval->atom = new Atom(8);return token::ORGANIC_ATOM_TOKEN; }
P  { yylval->atom = new Atom(15);return token::ORGANIC_ATOM_TOKEN; }
S  { yylval->atom = new Atom(16);return token::ORGANIC_ATOM_TOKEN; }
F  { yylval->atom = new Atom(9);return token::ORGANIC_ATOM_TOKEN; }
Cl { yylval->atom = new Atom(17);return token::ORGANIC_ATOM_TOKEN; }
Br { yylval->atom = new Atom(35);return token::ORGANIC_ATOM_TOKEN; }
I  { yylval->atom = new Atom(53);return token::ORGANIC_ATOM_TOKEN; }

H			{
				return token::H_TOKEN;
			}

b		    {	yylval->atom = new Atom ( 5 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
c		    {	yylval->atom = new Atom ( 6 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
n		    {	yylval->atom = new Atom( 7 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
o		    {	yylval->atom = new Atom( 8 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
p		    {	yylval->atom = new Atom( 15 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
s		    {	yylval->atom = new Atom( 16 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}

<IN_ATOM_STATE>si   {	yylval->atom = new Atom( 14 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
<IN_ATOM_STATE>as   {	yylval->atom = new Atom( 33 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
<IN_ATOM_STATE>se   {	yylval->atom = new Atom( 34 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}
<IN_ATOM_STATE>te   {	yylval->atom = new Atom( 52 );
			yylval->atom->setIsAromatic(true);
				return token::AROMATIC_ATOM_TOKEN;
			}

\* 	            {   yylval->atom = new Atom( 0 );
		            yylval->atom->setProp(common_properties::dummyLabel,
                                                        std::string("*"));
                                // must be ORGANIC_ATOM_TOKEN because
                                // we aren't in square brackets:
				return token::ORGANIC_ATOM_TOKEN;
			}

<IN_ATOM_STATE>\: 	{ return token::COLON_TOKEN; }

<IN_ATOM_STATE>\# 	{ return token::HASH_TOKEN; }


%{
  // The next block is a workaround for a pathlogy in the SMILES produced
  // by some Biovia tools
%}
<IN_ATOM_STATE>\'Rf\'	{ yylval->atom = new Atom(104); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Db\'	{ yylval->atom = new Atom(105); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Sg\'	{ yylval->atom = new Atom(106); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Bh\'	{ yylval->atom = new Atom(107); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Hs\'	{ yylval->atom = new Atom(108); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Mt\'	{ yylval->atom = new Atom(109); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Ds\'	{ yylval->atom = new Atom(110); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Rg\'	{ yylval->atom = new Atom(111); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Cn\'	{ yylval->atom = new Atom(112); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Nh\'	{ yylval->atom = new Atom(113); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Fl\'	{ yylval->atom = new Atom(114); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Mc\'	{ yylval->atom = new Atom(115); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Lv\'	{ yylval->atom = new Atom(116); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Ts\'	{ yylval->atom = new Atom(117); return token::ATOM_TOKEN; }
<IN_ATOM_STATE>\'Og\'	{ yylval->atom = new Atom(118); return token::ATOM_TOKEN; }

\=	{ yylval->bond = new Bond(Bond::DOUBLE);
	  return token::BOND_TOKEN; }
\#	{ yylval->bond = new Bond(Bond::TRIPLE);
	  return token::BOND_TOKEN; }
\:	{ yylval->bond = new Bond(Bond::AROMATIC);
	  yylval->bond->setIsAromatic(true);
	  return token::BOND_TOKEN; }
\$	{ yylval->bond = new Bond(Bond::QUADRUPLE);
	  return token::BOND_TOKEN; }
\-\>	{ yylval->bond = new Bond(Bond::DATIVER);
	  return token::BOND_TOKEN; }
\<\-	{ yylval->bond = new Bond(Bond::DATIVEL);
	  return token::BOND_TOKEN; }
\~	{ yylval->bond = new QueryBond();
	  yylval->bond->setQuery(makeBondNullQuery());
	  return token::BOND_TOKEN;  }

[\\]{1,2}    { yylval->bond = new Bond(Bond::UNSPECIFIED);
	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
	return token::BOND_TOKEN;  }

[\/]    { yylval->bond = new Bond(Bond::UNSPECIFIED);
	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
	return token::BOND_TOKEN;  }

\-			{ return token::MINUS_TOKEN; }

\+			{ return token::PLUS_TOKEN; }

\(       	{ return token::GROUP_OPEN_TOKEN; }
\)       	{ return token::GROUP_CLOSE_TOKEN; }


\[			{ BEGIN IN_ATOM_STATE; return token::ATOM_OPEN_TOKEN; }
<IN_ATOM_STATE>\]	{ BEGIN INITIAL; return token::ATOM_CLOSE_TOKEN; }

\.       	{ return token::SEPARATOR_TOKEN; }

\%              { return token::PERCENT_TOKEN; }

[0]		{ yylval->ival = 0; return token::ZERO_TOKEN; }
[1-9]		{ yylval->ival = yytext[0] - '0'; return token::NONZERO_DIGIT_TOKEN; }



\n		return 0;

<<EOF>>		{ return token::EOS_TOKEN; }
.		return yytext[0];

%%
