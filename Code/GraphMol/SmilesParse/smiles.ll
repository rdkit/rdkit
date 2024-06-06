%{

// $Id$
//
//  Copyright (C) 2001-2010 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

#include "SmilesScanner.h"

using smiles_parser::ast_parser::Parser;
using smiles_parser::ast_parser::Scanner;

#undef YY_DECL
#define YY_DECL int Scanner::lex(Parser::semantic_type* const lval, Parser::location_type* location)

// constructing the token this way allows us maintain the relationship between
// the input smiles and the current token
#define YY_USER_ACTION \
    lval->build<std::string_view>(input().substr(location->begin.column, yyleng)); \
    location->begin += yyleng; location->columns(yyleng);

using token = Parser::token_kind_type;
%}

%option c++
%option full
%option never-interactive
%option noyywrap
%option yylineno
%option yyclass="smiles_parser::ast_parser::Scanner"

%s IN_ATOM_STATE

%%

<*>[0-9]+ { return token::NUMBER; }

<*>H { return token::HYDROGEN; }
<*>Cl|Br|B|C|N|O|P|S|F|I|b|c|n|o|s|p { return token::ORGANIC_ATOM; }

<*>[a-z]  { return token::AROMATIC_SYMBOL; }
<IN_ATOM_STATE>si|as|se|te { return token::AROMATIC_SYMBOL; }

<IN_ATOM_STATE>[A-Z][a-z]*? {  return token::ELEMENT_SYMBOL ; }
@[' ']*[A-Z][A-Z] {  return token::CHIRAL_TAG;}

<INITIAL>\[  { BEGIN IN_ATOM_STATE;  return yytext[0]; }
<IN_ATOM_STATE>\]	{ BEGIN INITIAL; return yytext[0]; }

<*>.		{  return yytext[0]; }

%%
