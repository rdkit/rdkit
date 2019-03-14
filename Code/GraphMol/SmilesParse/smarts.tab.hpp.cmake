/* A Bison parser, made by GNU Bison 3.1.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018 Free Software Foundation, Inc.

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

#ifndef YY_YYSMARTS_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
# define YY_YYSMARTS_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yysmarts_debug;
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
    ORGANIC_ATOM_TOKEN = 262,
    ATOM_TOKEN = 263,
    SIMPLE_ATOM_QUERY_TOKEN = 264,
    COMPLEX_ATOM_QUERY_TOKEN = 265,
    RINGSIZE_ATOM_QUERY_TOKEN = 266,
    RINGBOND_ATOM_QUERY_TOKEN = 267,
    IMPLICIT_H_ATOM_QUERY_TOKEN = 268,
    HYB_TOKEN = 269,
    HETERONEIGHBOR_ATOM_QUERY_TOKEN = 270,
    ALIPHATIC = 271,
    ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN = 272,
    ZERO_TOKEN = 273,
    NONZERO_DIGIT_TOKEN = 274,
    GROUP_OPEN_TOKEN = 275,
    GROUP_CLOSE_TOKEN = 276,
    SEPARATOR_TOKEN = 277,
    RANGE_OPEN_TOKEN = 278,
    RANGE_CLOSE_TOKEN = 279,
    HASH_TOKEN = 280,
    MINUS_TOKEN = 281,
    PLUS_TOKEN = 282,
    CHIRAL_MARKER_TOKEN = 283,
    CHI_CLASS_TOKEN = 284,
    CHI_CLASS_OH_TOKEN = 285,
    H_TOKEN = 286,
    AT_TOKEN = 287,
    PERCENT_TOKEN = 288,
    ATOM_OPEN_TOKEN = 289,
    ATOM_CLOSE_TOKEN = 290,
    NOT_TOKEN = 291,
    AND_TOKEN = 292,
    OR_TOKEN = 293,
    SEMI_TOKEN = 294,
    BEGIN_RECURSE = 295,
    END_RECURSE = 296,
    COLON_TOKEN = 297,
    UNDERSCORE_TOKEN = 298,
    BOND_TOKEN = 299,
    EOS_TOKEN = 300
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 62 "smarts.yy" /* yacc.c:1913  */

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;

#line 107 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.hpp" /* yacc.c:1913  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token);
/* "%code provides" blocks.  */
#line 57 "smarts.yy" /* yacc.c:1913  */

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)

#line 124 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.hpp" /* yacc.c:1913  */

#endif /* !YY_YYSMARTS_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED  */
