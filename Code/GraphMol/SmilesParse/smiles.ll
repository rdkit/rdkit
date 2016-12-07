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
<IN_ATOM_STATE>Lr |
<IN_ATOM_STATE>Rf |
<IN_ATOM_STATE>Db |
<IN_ATOM_STATE>Sg |
<IN_ATOM_STATE>Bh |
<IN_ATOM_STATE>Hs |
<IN_ATOM_STATE>Mt |
<IN_ATOM_STATE>Ds |
<IN_ATOM_STATE>Rg |
<IN_ATOM_STATE>Cn |
<IN_ATOM_STATE>Uut |
<IN_ATOM_STATE>Fl |
<IN_ATOM_STATE>Uup |
<IN_ATOM_STATE>Lv {   yylval->atom = new Atom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
				return ATOM_TOKEN;
			}
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

[\=\#\:] {
              yylval->bond = new Bond();
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
\-\> {
    yylval->bond = new Bond(Bond::DATIVER);
    return BOND_TOKEN;
}
\<\- {
    yylval->bond = new Bond(Bond::DATIVEL);
    return BOND_TOKEN;
}

\~	{ yylval->bond = new QueryBond();
	yylval->bond->setQuery(makeBondNullQuery());
	return BOND_TOKEN;  }

[\\]{1,2}    { yylval->bond = new Bond(Bond::SINGLE);
	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
	return BOND_TOKEN;  }

[\/]    { yylval->bond = new Bond(Bond::SINGLE);
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
[1-9]		{ yylval->ival = atoi( yytext ); return NONZERO_DIGIT_TOKEN; }



\n		return 0;

<<EOF>>		{ return EOS_TOKEN; }
.		return yytext[0];

%%

#undef yysmiles_wrap
int yysmiles_wrap( void ) { return 1; }
