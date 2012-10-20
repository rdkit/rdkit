/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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
/* Tokens.  */
#define AROMATIC_ATOM_TOKEN 258
#define ORGANIC_ATOM_TOKEN 259
#define ATOM_TOKEN 260
#define SIMPLE_ATOM_QUERY_TOKEN 261
#define COMPLEX_ATOM_QUERY_TOKEN 262
#define RINGSIZE_ATOM_QUERY_TOKEN 263
#define HYB_TOKEN 264
#define ZERO_TOKEN 265
#define NONZERO_DIGIT_TOKEN 266
#define GROUP_OPEN_TOKEN 267
#define GROUP_CLOSE_TOKEN 268
#define SEPARATOR_TOKEN 269
#define HASH_TOKEN 270
#define MINUS_TOKEN 271
#define PLUS_TOKEN 272
#define CHIRAL_MARKER_TOKEN 273
#define CHI_CLASS_TOKEN 274
#define CHI_CLASS_OH_TOKEN 275
#define H_TOKEN 276
#define AT_TOKEN 277
#define PERCENT_TOKEN 278
#define ATOM_OPEN_TOKEN 279
#define ATOM_CLOSE_TOKEN 280
#define NOT_TOKEN 281
#define AND_TOKEN 282
#define OR_TOKEN 283
#define SEMI_TOKEN 284
#define BEGIN_RECURSE 285
#define END_RECURSE 286
#define COLON_TOKEN 287
#define UNDERSCORE_TOKEN 288
#define BOND_TOKEN 289
#define EOS_TOKEN 290




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 53 "smarts.yy"
{
  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;
}
/* Line 1529 of yacc.c.  */
#line 126 "/Users/landrgr1/RDKit_trunk/Code/GraphMol/SmilesParse/smarts.tab.hpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



