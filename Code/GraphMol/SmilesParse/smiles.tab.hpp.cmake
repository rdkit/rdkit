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

#ifndef YY_YYSMILES_USERS_FAARA_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED
# define YY_YYSMILES_USERS_FAARA_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yysmiles_debug;
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
    START_MOL = 258,               /* START_MOL  */
    START_ATOM = 259,              /* START_ATOM  */
    START_BOND = 260,              /* START_BOND  */
    AROMATIC_ATOM_TOKEN = 261,     /* AROMATIC_ATOM_TOKEN  */
    ATOM_TOKEN = 262,              /* ATOM_TOKEN  */
    ORGANIC_ATOM_TOKEN = 263,      /* ORGANIC_ATOM_TOKEN  */
    NONZERO_DIGIT_TOKEN = 264,     /* NONZERO_DIGIT_TOKEN  */
    ZERO_TOKEN = 265,              /* ZERO_TOKEN  */
    GROUP_OPEN_TOKEN = 266,        /* GROUP_OPEN_TOKEN  */
    GROUP_CLOSE_TOKEN = 267,       /* GROUP_CLOSE_TOKEN  */
    SEPARATOR_TOKEN = 268,         /* SEPARATOR_TOKEN  */
    LOOP_CONNECTOR_TOKEN = 269,    /* LOOP_CONNECTOR_TOKEN  */
    MINUS_TOKEN = 270,             /* MINUS_TOKEN  */
    PLUS_TOKEN = 271,              /* PLUS_TOKEN  */
    H_TOKEN = 272,                 /* H_TOKEN  */
    AT_TOKEN = 273,                /* AT_TOKEN  */
    PERCENT_TOKEN = 274,           /* PERCENT_TOKEN  */
    COLON_TOKEN = 275,             /* COLON_TOKEN  */
    HASH_TOKEN = 276,              /* HASH_TOKEN  */
    BOND_TOKEN = 277,              /* BOND_TOKEN  */
    CHI_CLASS_TOKEN = 278,         /* CHI_CLASS_TOKEN  */
    ATOM_OPEN_TOKEN = 279,         /* ATOM_OPEN_TOKEN  */
    ATOM_CLOSE_TOKEN = 280,        /* ATOM_CLOSE_TOKEN  */
    EOS_TOKEN = 281                /* EOS_TOKEN  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{

  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  RDKit::Atom::ChiralType chiraltype;
  int                      ival;


};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif




int yysmiles_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::vector<std::pair<unsigned int, unsigned int>>& branchPoints, void *scanner, int& start_token, unsigned int& current_token_position);

/* "%code provides" blocks.  */

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token, unsigned int& current_token_position)


#endif /* !YY_YYSMILES_USERS_FAARA_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED  */
