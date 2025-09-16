/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_YYSLN_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SLNPARSE_SLN_TAB_HPP_INCLUDED
# define YY_YYSLN_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SLNPARSE_SLN_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yysln_debug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    TEXT_BLOCK = 258,              /* TEXT_BLOCK  */
    CHAR_TOKEN = 259,              /* CHAR_TOKEN  */
    DIGIT_TOKEN = 260,             /* DIGIT_TOKEN  */
    H_TOKEN = 261,                 /* H_TOKEN  */
    H_BRACKET_TOKEN = 262,         /* H_BRACKET_TOKEN  */
    H_ASTERIX_TOKEN = 263,         /* H_ASTERIX_TOKEN  */
    AT_TOKEN = 264,                /* AT_TOKEN  */
    ATOM_TOKEN = 265,              /* ATOM_TOKEN  */
    COMPARE_TOKEN = 266,           /* COMPARE_TOKEN  */
    OPEN_PAREN_TOKEN = 267,        /* OPEN_PAREN_TOKEN  */
    CLOSE_PAREN_TOKEN = 268,       /* CLOSE_PAREN_TOKEN  */
    OPEN_BRACKET_TOKEN = 269,      /* OPEN_BRACKET_TOKEN  */
    CLOSE_BRACKET_TOKEN = 270,     /* CLOSE_BRACKET_TOKEN  */
    OPEN_ANGLE_TOKEN = 271,        /* OPEN_ANGLE_TOKEN  */
    CLOSE_ANGLE_TOKEN = 272,       /* CLOSE_ANGLE_TOKEN  */
    SEPARATOR_TOKEN = 273,         /* SEPARATOR_TOKEN  */
    ASTERIX_TOKEN = 274,           /* ASTERIX_TOKEN  */
    EOS_TOKEN = 275,               /* EOS_TOKEN  */
    PLUS_TOKEN = 276,              /* PLUS_TOKEN  */
    MINUS_TOKEN = 277,             /* MINUS_TOKEN  */
    COLON_TOKEN = 278,             /* COLON_TOKEN  */
    EQUALS_TOKEN = 279,            /* EQUALS_TOKEN  */
    TILDE_TOKEN = 280,             /* TILDE_TOKEN  */
    HASH_TOKEN = 281,              /* HASH_TOKEN  */
    COMMA_TOKEN = 282,             /* COMMA_TOKEN  */
    NOT_TOKEN = 283,               /* NOT_TOKEN  */
    AND_TOKEN = 284,               /* AND_TOKEN  */
    OR_TOKEN = 285,                /* OR_TOKEN  */
    SEMI_TOKEN = 286,              /* SEMI_TOKEN  */
    CARET_EQUALS_TOKEN = 287,      /* CARET_EQUALS_TOKEN  */
    COLON_EQUALS_TOKEN = 288,      /* COLON_EQUALS_TOKEN  */
    RECURSE_TOKEN = 289,           /* RECURSE_TOKEN  */
    NEG_RECURSE_TOKEN = 290,       /* NEG_RECURSE_TOKEN  */
    ERROR_TOKEN = 291              /* ERROR_TOKEN  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 95 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"

  int                      mol_T;
  RDKit::Atom *            atom_T;
  RDKit::Bond *            bond_T;
  int                      ival_T;
  std::string*             text_T;
  char                     char_T;
  RDKit::SLNParse::AttribType       *attrib_T;
  RDKit::SLNParse::AttribListType   *attriblist_T;

#line 111 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.hpp"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif




int yysln_parse (const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner);


#endif /* !YY_YYSLN_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SLNPARSE_SLN_TAB_HPP_INCLUDED  */
