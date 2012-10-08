/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison interface for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2011 Free Software Foundation, Inc.
   
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


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     AROMATIC_ATOM_TOKEN = 258,
     ORGANIC_ATOM_TOKEN = 259,
     ATOM_TOKEN = 260,
     SIMPLE_ATOM_QUERY_TOKEN = 261,
     COMPLEX_ATOM_QUERY_TOKEN = 262,
     RINGSIZE_ATOM_QUERY_TOKEN = 263,
     HYB_TOKEN = 264,
     ZERO_TOKEN = 265,
     NONZERO_DIGIT_TOKEN = 266,
     GROUP_OPEN_TOKEN = 267,
     GROUP_CLOSE_TOKEN = 268,
     SEPARATOR_TOKEN = 269,
     HASH_TOKEN = 270,
     MINUS_TOKEN = 271,
     PLUS_TOKEN = 272,
     CHIRAL_MARKER_TOKEN = 273,
     CHI_CLASS_TOKEN = 274,
     CHI_CLASS_OH_TOKEN = 275,
     H_TOKEN = 276,
     AT_TOKEN = 277,
     PERCENT_TOKEN = 278,
     ATOM_OPEN_TOKEN = 279,
     ATOM_CLOSE_TOKEN = 280,
     NOT_TOKEN = 281,
     AND_TOKEN = 282,
     OR_TOKEN = 283,
     SEMI_TOKEN = 284,
     BEGIN_RECURSE = 285,
     END_RECURSE = 286,
     COLON_TOKEN = 287,
     UNDERSCORE_TOKEN = 288,
     BOND_TOKEN = 289,
     EOS_TOKEN = 290
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 2068 of yacc.c  */
#line 53 "smarts.yy"

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;



/* Line 2068 of yacc.c  */
#line 94 "/scratch/RDKit_trunk/Code/GraphMol/SmilesParse/smarts.tab.hpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif




