%{

  // $Id$
  //
  //  Copyright (C) 2001-2016 Randal Henne, Greg Landrum and Rational Discovery LLC
  //
  //   @@ All Rights Reserved  @@
  //
  //
%}


%require  "4.0"

%debug
%defines
%define api.namespace {smiles_parser::ast_parser}
%define api.parser.class {Parser}
%define api.value.type variant
%defines "smiles.tab.hpp"
%output "smiles.tab.cpp"
%define api.location.file "smiles_location.tab.hpp"

%language "c++"
%locations


%code requires {
namespace smiles_parser {
namespace ast_parser {

class Scanner;
class Builder;

} // namespace ast_parser
} // namespace smiles_parser
} // code requires

%parse-param { Scanner& token_scanner }
%parse-param { Builder& ast }
%initial-action { @$.begin.column = 0; };

%code {
#include <charconv>
#include <utility>

#include "SmilesScanner.h"
#include "SmilesASTParser.h"

#undef yylex
#define yylex token_scanner.lex

} // code

%token <std::string_view> HYDROGEN ORGANIC_ATOM AROMATIC_SYMBOL ELEMENT_SYMBOL CHIRAL_TAG NUMBER;

%type <int> minus_signs plus_signs;


%start mol

%%
mol: /* empty */ | chain;

chain: atom
   | chain connector atom
   | chain ring_bond
   | chain branch_open connector chain branch_close;

connector: /* empty */ | dot | bond;

dot: '.' { ast.add_dot(); }
branch_open: '('  { ast.open_branch(); }
branch_close: ')'  { ast.close_branch(); }

bond: '-' { ast.add_bond("-"); }
    | '=' { ast.add_bond("="); }
    | '#' { ast.add_bond("#"); }
    | ':' { ast.add_bond(":"); }
    | '$' { ast.add_bond("$"); }
    | '~' { ast.add_bond("~"); }
    | '/' { ast.add_bond("/"); }
    | '\\'    { ast.add_bond("\\"); }
    | '-' '>' { ast.add_bond("->"); }
    | '<' '-' { ast.add_bond("<-"); }
    | '\\' '\\' { ast.add_bond("\\\\"); }

atom: organic | dummy_atom | bracket_atom { ast.set_no_implicit_hs(); };

bracket_atom: '[' isotope chirality hcount charge ']'
            | '[' isotope chirality hcount charge ':' map_number ']';

isotope: symbol | NUMBER symbol { ast.add_isotope_num($1); }

symbol: organic | element_symbol | dummy_atom;

organic: ORGANIC_ATOM { ast.add_atom($1); }

element_symbol: ELEMENT_SYMBOL { ast.add_atom($1); }
              | HYDROGEN { ast.add_atom($1); }
              /* some Biovia tools produce SMILES of this form */
              | '\'' ELEMENT_SYMBOL '\'' { ast.add_atom($2); }
              | AROMATIC_SYMBOL  { ast.add_atom($1); }
              | '#' NUMBER { ast.add_atom($2); }
              ;

dummy_atom: '*' { ast.add_atom("*"); }

chirality: /* empty */
         | '@' { ast.add_chirality_tag("@"); }
         | '@' '@' { ast.add_chirality_tag("@@"); }
         | CHIRAL_TAG { ast.add_chirality_class($1); }
         | CHIRAL_TAG NUMBER { ast.add_chirality_class($1, $2); }


hcount: /* empty */
      | HYDROGEN { ast.add_explicit_h(1); }
      | HYDROGEN NUMBER { ast.add_explicit_h($2); };

charge: /* empty */
      | plus_signs { ast.add_atom_charge($1); }
      | minus_signs { ast.add_atom_charge($1); }
      | '+' NUMBER { ast.add_atom_charge($2); }
      | '-' NUMBER { ast.add_atom_charge($2); };

plus_signs: '+' { $$ = 1; } | plus_signs '+' { $$ = $1 + 1; };
minus_signs: '-' { $$ = -1; } | minus_signs '-' { $$ = $1 - 1; };

map_number: NUMBER { ast.add_atom_map_number($1); };

ring_bond: ring_number | bond ring_number;

ring_number: NUMBER { ast.add_ring($1); }
           | '%' NUMBER { ast.add_ring($2); }
           | '%' '(' NUMBER ')' { ast.add_ring($3); };
%%

void smiles_parser::ast_parser::Parser::error(const location& loc, const std::string& msg) {
    auto bad_token_position = loc.begin.column - 1;
    ast.save_error_message(bad_token_position, msg);
}
