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

#ifndef YY_YYSMARTS_MNT_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
# define YY_YYSMARTS_MNT_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
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
    AROMATIC_ATOM_TOKEN = 258,
    ORGANIC_ATOM_TOKEN = 259,
    ATOM_TOKEN = 260,
    SIMPLE_ATOM_QUERY_TOKEN = 261,
    COMPLEX_ATOM_QUERY_TOKEN = 262,
    RANGE_ATOM_QUERY_TOKEN = 263,
    RINGSIZE_ATOM_QUERY_TOKEN = 264,
    RINGBOND_ATOM_QUERY_TOKEN = 265,
    IMPLICIT_H_ATOM_QUERY_TOKEN = 266,
    HYB_TOKEN = 267,
    ZERO_TOKEN = 268,
    NONZERO_DIGIT_TOKEN = 269,
    GROUP_OPEN_TOKEN = 270,
    GROUP_CLOSE_TOKEN = 271,
    SEPARATOR_TOKEN = 272,
    RANGE_OPEN_TOKEN = 273,
    RANGE_CLOSE_TOKEN = 274,
    HASH_TOKEN = 275,
    MINUS_TOKEN = 276,
    PLUS_TOKEN = 277,
    CHIRAL_MARKER_TOKEN = 278,
    CHI_CLASS_TOKEN = 279,
    CHI_CLASS_OH_TOKEN = 280,
    H_TOKEN = 281,
    AT_TOKEN = 282,
    PERCENT_TOKEN = 283,
    ATOM_OPEN_TOKEN = 284,
    ATOM_CLOSE_TOKEN = 285,
    NOT_TOKEN = 286,
    AND_TOKEN = 287,
    OR_TOKEN = 288,
    SEMI_TOKEN = 289,
    BEGIN_RECURSE = 290,
    END_RECURSE = 291,
    COLON_TOKEN = 292,
    UNDERSCORE_TOKEN = 293,
    BOND_TOKEN = 294,
    EOS_TOKEN = 295
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 51 "smarts.yy" /* yacc.c:1909  */

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;

#line 102 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.hpp" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner);

#endif /* !YY_YYSMARTS_MNT_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED  */
