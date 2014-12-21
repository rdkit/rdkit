/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

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

#ifndef YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
# define YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
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
    RINGSIZE_ATOM_QUERY_TOKEN = 263,
    RINGBOND_ATOM_QUERY_TOKEN = 264,
    IMPLICIT_H_ATOM_QUERY_TOKEN = 265,
    HYB_TOKEN = 266,
    ZERO_TOKEN = 267,
    NONZERO_DIGIT_TOKEN = 268,
    GROUP_OPEN_TOKEN = 269,
    GROUP_CLOSE_TOKEN = 270,
    SEPARATOR_TOKEN = 271,
    HASH_TOKEN = 272,
    MINUS_TOKEN = 273,
    PLUS_TOKEN = 274,
    CHIRAL_MARKER_TOKEN = 275,
    CHI_CLASS_TOKEN = 276,
    CHI_CLASS_OH_TOKEN = 277,
    H_TOKEN = 278,
    AT_TOKEN = 279,
    PERCENT_TOKEN = 280,
    ATOM_OPEN_TOKEN = 281,
    ATOM_CLOSE_TOKEN = 282,
    NOT_TOKEN = 283,
    AND_TOKEN = 284,
    OR_TOKEN = 285,
    SEMI_TOKEN = 286,
    BEGIN_RECURSE = 287,
    END_RECURSE = 288,
    COLON_TOKEN = 289,
    UNDERSCORE_TOKEN = 290,
    BOND_TOKEN = 291,
    EOS_TOKEN = 292
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 51 "smarts.yy" /* yacc.c:1909  */

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;

#line 99 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.hpp" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner);

#endif /* !YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED  */
