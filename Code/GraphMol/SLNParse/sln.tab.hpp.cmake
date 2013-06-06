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
     TEXT_BLOCK = 258,
     CHAR_TOKEN = 259,
     DIGIT_TOKEN = 260,
     H_TOKEN = 261,
     H_BRACKET_TOKEN = 262,
     H_ASTERIX_TOKEN = 263,
     AT_TOKEN = 264,
     ATOM_TOKEN = 265,
     COMPARE_TOKEN = 266,
     OPEN_PAREN_TOKEN = 267,
     CLOSE_PAREN_TOKEN = 268,
     OPEN_BRACKET_TOKEN = 269,
     CLOSE_BRACKET_TOKEN = 270,
     OPEN_ANGLE_TOKEN = 271,
     CLOSE_ANGLE_TOKEN = 272,
     SEPARATOR_TOKEN = 273,
     ASTERIX_TOKEN = 274,
     EOS_TOKEN = 275,
     PLUS_TOKEN = 276,
     MINUS_TOKEN = 277,
     COLON_TOKEN = 278,
     EQUALS_TOKEN = 279,
     TILDE_TOKEN = 280,
     HASH_TOKEN = 281,
     COMMA_TOKEN = 282,
     NOT_TOKEN = 283,
     AND_TOKEN = 284,
     OR_TOKEN = 285,
     SEMI_TOKEN = 286,
     CARET_EQUALS_TOKEN = 287,
     COLON_EQUALS_TOKEN = 288,
     RECURSE_TOKEN = 289,
     NEG_RECURSE_TOKEN = 290,
     ERROR_TOKEN = 291
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 2068 of yacc.c  */
#line 88 "sln.yy"

  int                      mol_T;
  RDKit::Atom *            atom_T;
  RDKit::Bond *            bond_T;
  int                      ival_T;
  std::string*             text_T;
  char                     char_T;
  RDKit::SLNParse::AttribType       *attrib_T;
  RDKit::SLNParse::AttribListType   *attriblist_T;



/* Line 2068 of yacc.c  */
#line 99 "/scratch/RDKit_trunk/Code/GraphMol/SLNParse/sln.tab.hpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif




