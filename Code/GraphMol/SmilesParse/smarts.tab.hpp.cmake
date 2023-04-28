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

#ifndef YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
#define YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#if YYDEBUG
extern int yysmarts_debug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
#define YYTOKENTYPE
enum yytokentype {
  YYEMPTY = -2,
  YYEOF = 0,                             /* "end of file"  */
  YYerror = 256,                         /* error  */
  YYUNDEF = 257,                         /* "invalid token"  */
  START_MOL = 258,                       /* START_MOL  */
  START_ATOM = 259,                      /* START_ATOM  */
  START_BOND = 260,                      /* START_BOND  */
  AROMATIC_ATOM_TOKEN = 261,             /* AROMATIC_ATOM_TOKEN  */
  ORGANIC_ATOM_TOKEN = 262,              /* ORGANIC_ATOM_TOKEN  */
  ATOM_TOKEN = 263,                      /* ATOM_TOKEN  */
  SIMPLE_ATOM_QUERY_TOKEN = 264,         /* SIMPLE_ATOM_QUERY_TOKEN  */
  COMPLEX_ATOM_QUERY_TOKEN = 265,        /* COMPLEX_ATOM_QUERY_TOKEN  */
  RINGSIZE_ATOM_QUERY_TOKEN = 266,       /* RINGSIZE_ATOM_QUERY_TOKEN  */
  RINGBOND_ATOM_QUERY_TOKEN = 267,       /* RINGBOND_ATOM_QUERY_TOKEN  */
  IMPLICIT_H_ATOM_QUERY_TOKEN = 268,     /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
  HYB_TOKEN = 269,                       /* HYB_TOKEN  */
  HETERONEIGHBOR_ATOM_QUERY_TOKEN = 270, /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
  ALIPHATIC = 271,                       /* ALIPHATIC  */
  ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN =
      272,                   /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
  ZERO_TOKEN = 273,          /* ZERO_TOKEN  */
  NONZERO_DIGIT_TOKEN = 274, /* NONZERO_DIGIT_TOKEN  */
  GROUP_OPEN_TOKEN = 275,    /* GROUP_OPEN_TOKEN  */
  GROUP_CLOSE_TOKEN = 276,   /* GROUP_CLOSE_TOKEN  */
  SEPARATOR_TOKEN = 277,     /* SEPARATOR_TOKEN  */
  RANGE_OPEN_TOKEN = 278,    /* RANGE_OPEN_TOKEN  */
  RANGE_CLOSE_TOKEN = 279,   /* RANGE_CLOSE_TOKEN  */
  HASH_TOKEN = 280,          /* HASH_TOKEN  */
  MINUS_TOKEN = 281,         /* MINUS_TOKEN  */
  PLUS_TOKEN = 282,          /* PLUS_TOKEN  */
  H_TOKEN = 283,             /* H_TOKEN  */
  AT_TOKEN = 284,            /* AT_TOKEN  */
  PERCENT_TOKEN = 285,       /* PERCENT_TOKEN  */
  ATOM_OPEN_TOKEN = 286,     /* ATOM_OPEN_TOKEN  */
  ATOM_CLOSE_TOKEN = 287,    /* ATOM_CLOSE_TOKEN  */
  NOT_TOKEN = 288,           /* NOT_TOKEN  */
  AND_TOKEN = 289,           /* AND_TOKEN  */
  OR_TOKEN = 290,            /* OR_TOKEN  */
  SEMI_TOKEN = 291,          /* SEMI_TOKEN  */
  BEGIN_RECURSE = 292,       /* BEGIN_RECURSE  */
  END_RECURSE = 293,         /* END_RECURSE  */
  COLON_TOKEN = 294,         /* COLON_TOKEN  */
  UNDERSCORE_TOKEN = 295,    /* UNDERSCORE_TOKEN  */
  BOND_TOKEN = 296,          /* BOND_TOKEN  */
  CHI_CLASS_TOKEN = 297,     /* CHI_CLASS_TOKEN  */
  EOS_TOKEN = 298            /* EOS_TOKEN  */
};
typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if !defined YYSTYPE && !defined YYSTYPE_IS_DECLARED
union YYSTYPE {
#line 79 "smarts.yy"

  int moli;
  RDKit::QueryAtom *atom;
  RDKit::QueryBond *bond;
  RDKit::Atom::ChiralType chiraltype;
  int ival;

#line 115 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.hpp"
};
typedef union YYSTYPE YYSTYPE;
#define YYSTYPE_IS_TRIVIAL 1
#define YYSTYPE_IS_DECLARED 1
#endif

int yysmarts_parse(const char *input, std::vector<RDKit::RWMol *> *molList,
                   RDKit::Atom *&lastAtom, RDKit::Bond *&lastBond,
                   unsigned &numAtomsParsed, unsigned &numBondsParsed,
                   std::list<unsigned int> *branchPoints, void *scanner,
                   int &start_token);

/* "%code provides" blocks.  */
#line 74 "smarts.yy"

#define YY_DECL \
  int yylex(YYSTYPE *yylval_param, yyscan_t yyscanner, int &start_token)

#line 134 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.hpp"

#endif /* !YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED \
        */
