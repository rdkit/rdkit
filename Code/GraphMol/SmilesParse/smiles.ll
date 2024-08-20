%option reentrant
%option bison-bridge
%option noyywrap
%option never-interactive
%option nodefault
%option nostdinit

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

using namespace RDKit;

#define YY_FATAL_ERROR(msg) smiles_lexer_error(msg)

void smiles_lexer_error(const char *msg) {
     BOOST_LOG(rdErrorLog) << msg<<std::endl;
     throw ValueErrorException(msg);
}

size_t setup_smiles_string(const std::string &text,yyscan_t yyscanner){
//  YY_BUFFER_STATE buff=yysmiles__scan_string(text.c_str()+pos,yyscanner);
  // Faster implementation of yysmiles__scan_string that handles trimming
  YY_BUFFER_STATE b;
  char *buf;
  yyconst char * yybytes = text.c_str();
  yy_size_t _yybytes_len=text.size(), n, start, end;
  /* Get memory for full buffer, including space for trailing EOB's. */
  n = _yybytes_len + 2;
  buf = (char *) yysmiles_alloc(n ,yyscanner );
  if ( ! buf )
    smiles_lexer_error( "out of dynamic memory in yysmiles__scan_bytes()" );

  // ltrim

  for(start = 0 ; start < _yybytes_len; ++start) {
    if (yybytes[start] > 32) break;
  }
  for(end = _yybytes_len ; end > start; --end) {
    if (yybytes[end] > 32) break;
  }

  _yybytes_len = end-start+1;
  n = _yybytes_len + 2;
  memcpy(buf, yybytes+start, _yybytes_len);


  buf[_yybytes_len] = buf[_yybytes_len+1] = YY_END_OF_BUFFER_CHAR;

  b = yysmiles__scan_buffer(buf,n ,yyscanner);
  if ( ! b )
    smiles_lexer_error( "bad buffer in yysmiles__scan_bytes()" );

  /* It's okay to grow etc. this buffer, and we should throw it
   * away when we're done.
   */
  b->yy_is_our_buffer = 1;


  POSTCONDITION(b,"invalid buffer");
  return start;

}
%}

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

@[' ']*TH { yylval->chiraltype = Atom::ChiralType::CHI_TETRAHEDRAL; return CHI_CLASS_TOKEN; }
@[' ']*AL { yylval->chiraltype = Atom::ChiralType::CHI_ALLENE; return CHI_CLASS_TOKEN; }
@[' ']*SP { yylval->chiraltype = Atom::ChiralType::CHI_SQUAREPLANAR; return CHI_CLASS_TOKEN; }
@[' ']*TB { yylval->chiraltype = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL; return CHI_CLASS_TOKEN; }
@[' ']*OH { yylval->chiraltype = Atom::ChiralType::CHI_OCTAHEDRAL; return CHI_CLASS_TOKEN; }

@		{ return AT_TOKEN; }


<IN_ATOM_STATE>He	{ yylval->atom = new Atom(2); return ATOM_TOKEN; }
<IN_ATOM_STATE>Li	{ yylval->atom = new Atom(3); return ATOM_TOKEN; }
<IN_ATOM_STATE>Be	{ yylval->atom = new Atom(4); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ne	{ yylval->atom = new Atom(10); return ATOM_TOKEN; }
<IN_ATOM_STATE>Na	{ yylval->atom = new Atom(11); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mg	{ yylval->atom = new Atom(12); return ATOM_TOKEN; }
<IN_ATOM_STATE>Al	{ yylval->atom = new Atom(13); return ATOM_TOKEN; }
<IN_ATOM_STATE>Si	{ yylval->atom = new Atom(14); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ar	{ yylval->atom = new Atom(18); return ATOM_TOKEN; }
<IN_ATOM_STATE>K	{ yylval->atom = new Atom(19); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ca	{ yylval->atom = new Atom(20); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sc	{ yylval->atom = new Atom(21); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ti	{ yylval->atom = new Atom(22); return ATOM_TOKEN; }
<IN_ATOM_STATE>V	{ yylval->atom = new Atom(23); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cr	{ yylval->atom = new Atom(24); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mn	{ yylval->atom = new Atom(25); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fe	{ yylval->atom = new Atom(26); return ATOM_TOKEN; }
<IN_ATOM_STATE>Co	{ yylval->atom = new Atom(27); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ni	{ yylval->atom = new Atom(28); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cu	{ yylval->atom = new Atom(29); return ATOM_TOKEN; }
<IN_ATOM_STATE>Zn	{ yylval->atom = new Atom(30); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ga	{ yylval->atom = new Atom(31); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ge	{ yylval->atom = new Atom(32); return ATOM_TOKEN; }
<IN_ATOM_STATE>As	{ yylval->atom = new Atom(33); return ATOM_TOKEN; }
<IN_ATOM_STATE>Se	{ yylval->atom = new Atom(34); return ATOM_TOKEN; }
<IN_ATOM_STATE>Kr	{ yylval->atom = new Atom(36); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rb	{ yylval->atom = new Atom(37); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sr	{ yylval->atom = new Atom(38); return ATOM_TOKEN; }
<IN_ATOM_STATE>Y	{ yylval->atom = new Atom(39); return ATOM_TOKEN; }
<IN_ATOM_STATE>Zr	{ yylval->atom = new Atom(40); return ATOM_TOKEN; }
<IN_ATOM_STATE>Nb	{ yylval->atom = new Atom(41); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mo	{ yylval->atom = new Atom(42); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tc	{ yylval->atom = new Atom(43); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ru	{ yylval->atom = new Atom(44); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rh	{ yylval->atom = new Atom(45); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pd	{ yylval->atom = new Atom(46); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ag	{ yylval->atom = new Atom(47); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cd	{ yylval->atom = new Atom(48); return ATOM_TOKEN; }
<IN_ATOM_STATE>In	{ yylval->atom = new Atom(49); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sn	{ yylval->atom = new Atom(50); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sb	{ yylval->atom = new Atom(51); return ATOM_TOKEN; }
<IN_ATOM_STATE>Te	{ yylval->atom = new Atom(52); return ATOM_TOKEN; }
<IN_ATOM_STATE>Xe	{ yylval->atom = new Atom(54); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cs	{ yylval->atom = new Atom(55); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ba	{ yylval->atom = new Atom(56); return ATOM_TOKEN; }
<IN_ATOM_STATE>La	{ yylval->atom = new Atom(57); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ce	{ yylval->atom = new Atom(58); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pr	{ yylval->atom = new Atom(59); return ATOM_TOKEN; }
<IN_ATOM_STATE>Nd	{ yylval->atom = new Atom(60); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pm	{ yylval->atom = new Atom(61); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sm	{ yylval->atom = new Atom(62); return ATOM_TOKEN; }
<IN_ATOM_STATE>Eu	{ yylval->atom = new Atom(63); return ATOM_TOKEN; }
<IN_ATOM_STATE>Gd	{ yylval->atom = new Atom(64); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tb	{ yylval->atom = new Atom(65); return ATOM_TOKEN; }
<IN_ATOM_STATE>Dy	{ yylval->atom = new Atom(66); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ho	{ yylval->atom = new Atom(67); return ATOM_TOKEN; }
<IN_ATOM_STATE>Er	{ yylval->atom = new Atom(68); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tm	{ yylval->atom = new Atom(69); return ATOM_TOKEN; }
<IN_ATOM_STATE>Yb	{ yylval->atom = new Atom(70); return ATOM_TOKEN; }
<IN_ATOM_STATE>Lu	{ yylval->atom = new Atom(71); return ATOM_TOKEN; }
<IN_ATOM_STATE>Hf	{ yylval->atom = new Atom(72); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ta	{ yylval->atom = new Atom(73); return ATOM_TOKEN; }
<IN_ATOM_STATE>W	{ yylval->atom = new Atom(74); return ATOM_TOKEN; }
<IN_ATOM_STATE>Re	{ yylval->atom = new Atom(75); return ATOM_TOKEN; }
<IN_ATOM_STATE>Os	{ yylval->atom = new Atom(76); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ir	{ yylval->atom = new Atom(77); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pt	{ yylval->atom = new Atom(78); return ATOM_TOKEN; }
<IN_ATOM_STATE>Au	{ yylval->atom = new Atom(79); return ATOM_TOKEN; }
<IN_ATOM_STATE>Hg	{ yylval->atom = new Atom(80); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tl	{ yylval->atom = new Atom(81); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pb	{ yylval->atom = new Atom(82); return ATOM_TOKEN; }
<IN_ATOM_STATE>Bi	{ yylval->atom = new Atom(83); return ATOM_TOKEN; }
<IN_ATOM_STATE>Po	{ yylval->atom = new Atom(84); return ATOM_TOKEN; }
<IN_ATOM_STATE>At	{ yylval->atom = new Atom(85); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rn	{ yylval->atom = new Atom(86); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fr	{ yylval->atom = new Atom(87); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ra	{ yylval->atom = new Atom(88); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ac	{ yylval->atom = new Atom(89); return ATOM_TOKEN; }
<IN_ATOM_STATE>Th	{ yylval->atom = new Atom(90); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pa	{ yylval->atom = new Atom(91); return ATOM_TOKEN; }
<IN_ATOM_STATE>U	{ yylval->atom = new Atom(92); return ATOM_TOKEN; }
<IN_ATOM_STATE>Np	{ yylval->atom = new Atom(93); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pu	{ yylval->atom = new Atom(94); return ATOM_TOKEN; }
<IN_ATOM_STATE>Am	{ yylval->atom = new Atom(95); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cm	{ yylval->atom = new Atom(96); return ATOM_TOKEN; }
<IN_ATOM_STATE>Bk	{ yylval->atom = new Atom(97); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cf	{ yylval->atom = new Atom(98); return ATOM_TOKEN; }
<IN_ATOM_STATE>Es	{ yylval->atom = new Atom(99); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fm	{ yylval->atom = new Atom(100); return ATOM_TOKEN; }
<IN_ATOM_STATE>Md	{ yylval->atom = new Atom(101); return ATOM_TOKEN; }
<IN_ATOM_STATE>No	{ yylval->atom = new Atom(102); return ATOM_TOKEN; }
<IN_ATOM_STATE>Lr	{ yylval->atom = new Atom(103); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rf	{ yylval->atom = new Atom(104); return ATOM_TOKEN; }
<IN_ATOM_STATE>Db	{ yylval->atom = new Atom(105); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sg	{ yylval->atom = new Atom(106); return ATOM_TOKEN; }
<IN_ATOM_STATE>Bh	{ yylval->atom = new Atom(107); return ATOM_TOKEN; }
<IN_ATOM_STATE>Hs	{ yylval->atom = new Atom(108); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mt	{ yylval->atom = new Atom(109); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ds	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rg	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cn	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
<IN_ATOM_STATE>Nh	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fl	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mc	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
<IN_ATOM_STATE>Lv	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ts	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
<IN_ATOM_STATE>Og	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }

<IN_ATOM_STATE>Uun	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuu	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uub	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uut	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuq	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uup	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuh	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uus	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuo	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }

B  { yylval->atom = new Atom(5);return ORGANIC_ATOM_TOKEN; }
C  { yylval->atom = new Atom(6);return ORGANIC_ATOM_TOKEN; }
N  { yylval->atom = new Atom(7);return ORGANIC_ATOM_TOKEN; }
O  { yylval->atom = new Atom(8);return ORGANIC_ATOM_TOKEN; }
P  { yylval->atom = new Atom(15);return ORGANIC_ATOM_TOKEN; }
S  { yylval->atom = new Atom(16);return ORGANIC_ATOM_TOKEN; }
F  { yylval->atom = new Atom(9);return ORGANIC_ATOM_TOKEN; }
Cl { yylval->atom = new Atom(17);return ORGANIC_ATOM_TOKEN; }
Br { yylval->atom = new Atom(35);return ORGANIC_ATOM_TOKEN; }
I  { yylval->atom = new Atom(53);return ORGANIC_ATOM_TOKEN; }

H			{
				return H_TOKEN;
			}

b		    {	yylval->atom = new Atom ( 5 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN;
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

<IN_ATOM_STATE>si   {	yylval->atom = new Atom( 14 );
			yylval->atom->setIsAromatic(true);
				return AROMATIC_ATOM_TOKEN;
			}
<IN_ATOM_STATE>as   {	yylval->atom = new Atom( 33 );
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
		            yylval->atom->setProp(common_properties::dummyLabel,
                                                        std::string("*"));
                                // must be ORGANIC_ATOM_TOKEN because
                                // we aren't in square brackets:
				return ORGANIC_ATOM_TOKEN;
			}

<IN_ATOM_STATE>\: 	{ return COLON_TOKEN; }

<IN_ATOM_STATE>\# 	{ return HASH_TOKEN; }


%{
  // The next block is a workaround for a pathlogy in the SMILES produced
  // by some Biovia tools
%}
<IN_ATOM_STATE>\'Rf\'	{ yylval->atom = new Atom(104); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Db\'	{ yylval->atom = new Atom(105); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Sg\'	{ yylval->atom = new Atom(106); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Bh\'	{ yylval->atom = new Atom(107); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Hs\'	{ yylval->atom = new Atom(108); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Mt\'	{ yylval->atom = new Atom(109); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Ds\'	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Rg\'	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Cn\'	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Nh\'	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Fl\'	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Mc\'	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Lv\'	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Ts\'	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Og\'	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }

\=	{ yylval->bond = new Bond(Bond::DOUBLE);
	  return BOND_TOKEN; }
\#	{ yylval->bond = new Bond(Bond::TRIPLE);
	  return BOND_TOKEN; }
\:	{ yylval->bond = new Bond(Bond::AROMATIC);
	  yylval->bond->setIsAromatic(true);
	  return BOND_TOKEN; }
\$	{ yylval->bond = new Bond(Bond::QUADRUPLE);
	  return BOND_TOKEN; }
\-\>	{ yylval->bond = new Bond(Bond::DATIVER);
	  return BOND_TOKEN; }
\<\-	{ yylval->bond = new Bond(Bond::DATIVEL);
	  return BOND_TOKEN; }
\~	{ yylval->bond = new QueryBond();
	  yylval->bond->setQuery(makeBondNullQuery());
	  return BOND_TOKEN;  }

[\\]{1,2}    { yylval->bond = new Bond(Bond::UNSPECIFIED);
	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
	return BOND_TOKEN;  }

[\/]    { yylval->bond = new Bond(Bond::UNSPECIFIED);
	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
	return BOND_TOKEN;  }

\-			{ return MINUS_TOKEN; }

\+			{ return PLUS_TOKEN; }

\(       	{ return GROUP_OPEN_TOKEN; }
\)       	{ return GROUP_CLOSE_TOKEN; }


\[			{ BEGIN IN_ATOM_STATE; return ATOM_OPEN_TOKEN; }
<IN_ATOM_STATE>\]	{ BEGIN INITIAL; return ATOM_CLOSE_TOKEN; }

\.       	{ return SEPARATOR_TOKEN; }

\%              { return PERCENT_TOKEN; }

[0]		{ yylval->ival = 0; return ZERO_TOKEN; }
[1-9]		{ yylval->ival = yytext[0] - '0'; return NONZERO_DIGIT_TOKEN; }



\n		return 0;

<<EOF>>		{ return EOS_TOKEN; }
.		return yytext[0];

%%

#undef yysmiles_wrap
int yysmiles_wrap( void ) { return 1; }
