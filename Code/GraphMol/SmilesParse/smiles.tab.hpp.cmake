/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YYSMILES_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED
# define YY_YYSMILES_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yysmiles_debug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    START_MOL = 258,
    START_ATOM = 259,
    START_BOND = 260,
    AROMATIC_ATOM_TOKEN = 261,
    ATOM_TOKEN = 262,
    ORGANIC_ATOM_TOKEN = 263,
    NONZERO_DIGIT_TOKEN = 264,
    ZERO_TOKEN = 265,
    GROUP_OPEN_TOKEN = 266,
    GROUP_CLOSE_TOKEN = 267,
    SEPARATOR_TOKEN = 268,
    LOOP_CONNECTOR_TOKEN = 269,
    MINUS_TOKEN = 270,
    PLUS_TOKEN = 271,
    H_TOKEN = 272,
    AT_TOKEN = 273,
    PERCENT_TOKEN = 274,
    COLON_TOKEN = 275,
    HASH_TOKEN = 276,
    BOND_TOKEN = 277,
    CHI_CLASS_TOKEN = 278,
    ATOM_OPEN_TOKEN = 279,
    ATOM_CLOSE_TOKEN = 280,
    EOS_TOKEN = 281
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 82 "smiles.yy" /* yacc.c:1909  */

  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  RDKit::Atom::ChiralType chiraltype;
  int                      ival;

#line 89 "C:/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.hpp" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmiles_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token);
/* "%code provides" blocks.  */
#line 77 "smiles.yy" /* yacc.c:1909  */

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)

#line 106 "C:/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.hpp" /* yacc.c:1909  */

#endif /* !YY_YYSMILES_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED  */
