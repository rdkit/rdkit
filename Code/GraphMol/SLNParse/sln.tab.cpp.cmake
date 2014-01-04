/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse yysln_parse
#define yylex   yysln_lex
#define yyerror yysln_error
#define yylval  yysln_lval
#define yychar  yysln_char
#define yydebug yysln_debug
#define yynerrs yysln_nerrs


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
/* Tokens.  */
#define TEXT_BLOCK 258
#define CHAR_TOKEN 259
#define DIGIT_TOKEN 260
#define H_TOKEN 261
#define H_BRACKET_TOKEN 262
#define H_ASTERIX_TOKEN 263
#define AT_TOKEN 264
#define ATOM_TOKEN 265
#define COMPARE_TOKEN 266
#define OPEN_PAREN_TOKEN 267
#define CLOSE_PAREN_TOKEN 268
#define OPEN_BRACKET_TOKEN 269
#define CLOSE_BRACKET_TOKEN 270
#define OPEN_ANGLE_TOKEN 271
#define CLOSE_ANGLE_TOKEN 272
#define SEPARATOR_TOKEN 273
#define ASTERIX_TOKEN 274
#define EOS_TOKEN 275
#define PLUS_TOKEN 276
#define MINUS_TOKEN 277
#define COLON_TOKEN 278
#define EQUALS_TOKEN 279
#define TILDE_TOKEN 280
#define HASH_TOKEN 281
#define COMMA_TOKEN 282
#define NOT_TOKEN 283
#define AND_TOKEN 284
#define OR_TOKEN 285
#define SEMI_TOKEN 286
#define CARET_EQUALS_TOKEN 287
#define COLON_EQUALS_TOKEN 288
#define RECURSE_TOKEN 289
#define NEG_RECURSE_TOKEN 290
#define ERROR_TOKEN 291




/* Copy the first part of user declarations.  */
#line 3 "sln.yy"


  // $Id$
  //
  //  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
  //  All rights reserved.
  // 
  // Redistribution and use in source and binary forms, with or without
  // modification, are permitted provided that the following conditions are
  // met: 
  //
  //     * Redistributions of source code must retain the above copyright 
  //       notice, this list of conditions and the following disclaimer.
  //     * Redistributions in binary form must reproduce the above
  //       copyright notice, this list of conditions and the following 
  //       disclaimer in the documentation and/or other materials provided 
  //       with the distribution.
  //     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
  //       nor the names of its contributors may be used to endorse or promote 
  //       products derived from this software without specific prior
  //       written permission.
  //
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  //
  // Created by Greg Landrum, September 2006
  //

#include <cstring>
#include <cstdio>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SLNParse/SLNParseOps.h>
#include <GraphMol/SLNParse/SLNAttribs.h>
#include <GraphMol/SLNParse/SLNParse.h>
#include <RDGeneral/RDLog.h>
#include "sln.tab.hpp"

int yysln_lex(YYSTYPE *,void *);

#define YYDEBUG 1

void
yysln_error( const char *input,
             std::vector<RDKit::RWMol *> *ms,bool doQ,
	     void *scanner,const char * msg )
{
  BOOST_LOG(rdErrorLog)<<"SLN Parse Error: "<<msg<<" while parsing: "<<input<<std::endl;
}

 namespace SLNParse = RDKit::SLNParse;

#define YYPRINT(file, type, value)   yyprint (file, type, value)

static void
yyprint (FILE *file, int type, YYSTYPE value)
{
  if (type == TEXT_BLOCK)
    fprintf (file, " %s", value.text_T->c_str());
  else fprintf (file, " %d", type);
}




/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 88 "sln.yy"
{
  int                      mol_T;
  RDKit::Atom *            atom_T;
  RDKit::Bond *            bond_T;
  int                      ival_T;
  std::string*             text_T;
  char                     char_T;
  RDKit::SLNParse::AttribType       *attrib_T;
  RDKit::SLNParse::AttribListType   *attriblist_T;
}
/* Line 193 of yacc.c.  */
#line 265 "/Users/landrgr1/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 278 "/Users/landrgr1/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  25
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   219

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  16
/* YYNRULES -- Number of rules.  */
#define YYNRULES  72
/* YYNRULES -- Number of states.  */
#define YYNSTATES  116

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   291

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,    10,    14,    17,    20,    22,    26,
      28,    30,    33,    37,    41,    46,    51,    57,    63,    70,
      74,    79,    85,    91,    98,   105,   113,   115,   117,   120,
     124,   126,   130,   136,   140,   142,   145,   150,   157,   162,
     164,   169,   171,   176,   178,   181,   183,   185,   187,   189,
     193,   197,   201,   203,   205,   209,   211,   214,   218,   220,
     223,   225,   228,   230,   232,   235,   238,   242,   244,   248,
     252,   256,   258
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      38,     0,    -1,    39,    -1,    38,    16,    48,    17,    -1,
      38,     1,    20,    -1,    38,    20,    -1,     1,    20,    -1,
      40,    -1,    39,    18,    40,    -1,     6,    -1,    41,    -1,
      40,    41,    -1,    40,    44,    41,    -1,    40,     9,    52,
      -1,    40,    44,     9,    52,    -1,    40,    12,    40,    13,
      -1,    40,    12,    44,    40,    13,    -1,    40,    12,     9,
      52,    13,    -1,    40,    12,    44,     9,    52,    13,    -1,
      40,    44,     6,    -1,    40,     9,    52,     6,    -1,    40,
      44,     9,    52,     6,    -1,    40,    12,    40,    13,     6,
      -1,    40,    12,    44,    40,    13,     6,    -1,    40,    12,
       9,    52,    13,     6,    -1,    40,    12,    44,     9,    52,
      13,     6,    -1,    42,    -1,    43,    -1,    43,     6,    -1,
      43,     6,     5,    -1,     8,    -1,     7,    52,    15,    -1,
       7,    52,    23,    47,    15,    -1,     7,    47,    15,    -1,
      10,    -1,    43,    19,    -1,    43,    14,    52,    15,    -1,
      43,    14,    52,    23,    47,    15,    -1,    43,    14,    47,
      15,    -1,    45,    -1,    45,    14,    47,    15,    -1,    25,
      -1,    25,    14,    47,    15,    -1,    46,    -1,    45,    46,
      -1,    22,    -1,    24,    -1,    26,    -1,    23,    -1,    47,
      29,    49,    -1,    47,    30,    49,    -1,    47,    31,    49,
      -1,    49,    -1,    51,    -1,    48,    31,    51,    -1,     3,
      -1,    28,    49,    -1,     3,    11,     3,    -1,    21,    -1,
      21,     5,    -1,    22,    -1,    22,     5,    -1,    19,    -1,
      50,    -1,    34,    38,    -1,    35,    38,    -1,    50,    27,
      38,    -1,     3,    -1,     3,    24,     3,    -1,     3,    33,
       3,    -1,     3,    32,     3,    -1,     5,    -1,    52,     5,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   133,   133,   134,   142,   154,   157,   171,   172,   177,
     187,   190,   194,   198,   202,   207,   211,   217,   221,   225,
     236,   248,   261,   273,   287,   299,   313,   314,   315,   319,
     325,   333,   341,   351,   362,   363,   366,   370,   376,   384,
     385,   392,   397,   406,   407,   422,   431,   440,   450,   462,
     466,   470,   474,   481,   486,   492,   500,   504,   515,   521,
     527,   533,   539,   545,   550,   564,   580,   594,   602,   612,
     622,   634,   635
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "TEXT_BLOCK", "CHAR_TOKEN",
  "DIGIT_TOKEN", "H_TOKEN", "H_BRACKET_TOKEN", "H_ASTERIX_TOKEN",
  "AT_TOKEN", "ATOM_TOKEN", "COMPARE_TOKEN", "OPEN_PAREN_TOKEN",
  "CLOSE_PAREN_TOKEN", "OPEN_BRACKET_TOKEN", "CLOSE_BRACKET_TOKEN",
  "OPEN_ANGLE_TOKEN", "CLOSE_ANGLE_TOKEN", "SEPARATOR_TOKEN",
  "ASTERIX_TOKEN", "EOS_TOKEN", "PLUS_TOKEN", "MINUS_TOKEN", "COLON_TOKEN",
  "EQUALS_TOKEN", "TILDE_TOKEN", "HASH_TOKEN", "COMMA_TOKEN", "NOT_TOKEN",
  "AND_TOKEN", "OR_TOKEN", "SEMI_TOKEN", "CARET_EQUALS_TOKEN",
  "COLON_EQUALS_TOKEN", "RECURSE_TOKEN", "NEG_RECURSE_TOKEN",
  "ERROR_TOKEN", "$accept", "cmpd", "mol", "primmol", "atom", "hatom",
  "primatom", "bond", "primbond", "onebond", "attriblist",
  "ctabattriblist", "attrib", "recursivequery", "ctabattrib", "number", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    37,    38,    38,    38,    38,    38,    39,    39,    40,
      40,    40,    40,    40,    40,    40,    40,    40,    40,    40,
      40,    40,    40,    40,    40,    40,    41,    41,    41,    41,
      42,    42,    42,    42,    43,    43,    43,    43,    43,    44,
      44,    44,    44,    45,    45,    46,    46,    46,    46,    47,
      47,    47,    47,    48,    48,    49,    49,    49,    49,    49,
      49,    49,    49,    49,    50,    50,    50,    51,    51,    51,
      51,    52,    52
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     4,     3,     2,     2,     1,     3,     1,
       1,     2,     3,     3,     4,     4,     5,     5,     6,     3,
       4,     5,     5,     6,     6,     7,     1,     1,     2,     3,
       1,     3,     5,     3,     1,     2,     4,     6,     4,     1,
       4,     1,     4,     1,     2,     1,     1,     1,     1,     3,
       3,     3,     1,     1,     3,     1,     2,     3,     1,     2,
       1,     2,     1,     1,     2,     2,     3,     1,     3,     3,
       3,     1,     2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     9,     0,    30,    34,     0,     2,     7,    10,
      26,    27,     6,    55,    71,    62,    58,    60,     0,     0,
       0,     0,    52,    63,     0,     1,     0,     0,     5,     0,
       0,     0,    45,    48,    46,    41,    47,    11,     0,    39,
      43,    28,     0,    35,     0,    59,    61,    56,     0,     0,
      33,     0,     0,     0,     0,    72,    31,     0,     4,    67,
       0,    53,     8,    13,     0,     0,     0,     0,    19,     0,
      12,     0,    44,    29,     0,     0,    57,    49,    50,    51,
       0,     0,     0,     0,     0,     3,     0,    20,     0,    15,
       0,     0,     0,    14,     0,    38,    36,     0,    32,    68,
      70,    69,    54,    17,    22,     0,    16,    42,    21,    40,
       0,    24,    18,    23,    37,    25
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     6,     7,     8,     9,    10,    11,    38,    39,    40,
      21,    60,    22,    23,    61,    24
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -40
static const yytype_int16 yypact[] =
{
     194,   -13,   -40,    28,   -40,   -40,     5,    -6,   154,   -40,
     -40,    46,   -40,     3,   -40,   -40,    12,    15,    45,   194,
     194,     7,   -40,    17,     4,   -40,    26,    48,   -40,   209,
      50,   107,   -40,   -40,   -40,    39,   -40,   -40,   199,   170,
     -40,    54,    28,   -40,    65,   -40,   -40,   -40,    69,    75,
     -40,    45,    45,    45,   194,   -40,   -40,    45,   -40,   141,
      -2,   -40,   154,    18,    50,   127,   204,    45,   -40,    50,
     -40,    45,   -40,   -40,    57,    11,   -40,   -40,   -40,   -40,
      92,    95,    68,    71,    72,   -40,    48,   -40,    64,    86,
      50,   134,   139,    76,   152,   -40,   -40,    45,   -40,   -40,
     -40,   -40,   -40,    88,   -40,    96,    91,   -40,   -40,   -40,
     156,   -40,    97,   -40,   -40,   -40
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -40,    -9,   -40,   -27,    -8,   -40,   -40,    47,   -40,    79,
     -39,   -40,   -10,   -40,    25,   -29
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -67
static const yytype_int8 yytable[] =
{
      37,    63,    62,    74,    65,    25,    26,    12,    47,    55,
      48,    49,    29,    75,    44,    85,    55,    45,    81,    56,
      46,    27,    50,    55,    87,    28,    96,    57,    92,    86,
      70,    13,    94,    14,    97,    88,    51,    52,    53,    91,
      93,    77,    78,    79,    54,    80,    58,    15,    13,    16,
      17,    59,    41,    67,    37,    14,    18,    37,   110,    73,
      42,   105,    19,    20,    15,    43,    16,    17,    76,    55,
      26,    99,    95,    18,   100,   101,    26,   103,    66,    19,
      20,    55,   108,    37,   -64,    27,    51,    52,    53,    28,
     -65,    27,   104,    26,   111,    28,   -64,   113,   -64,   -64,
     -64,    55,   -65,   115,   -65,   -65,   -65,   -66,    27,   112,
      98,   102,    28,     2,     3,     4,    64,     5,    72,   -66,
       0,   -66,   -66,   -66,    51,    52,    53,     0,     0,    32,
      33,    34,    35,    36,     3,     4,    30,     5,     0,    31,
      89,     3,     4,    30,     5,     0,    31,   106,     0,    32,
      33,    34,    35,    36,   107,     0,    32,    33,    34,    35,
      36,     3,     4,    30,     5,    82,    31,   109,    51,    52,
      53,   114,     0,    83,    84,     0,    32,    33,    34,    35,
      36,    51,    52,    53,    71,    51,    52,    53,     0,     0,
       0,     0,    32,    33,    34,     1,    36,     0,     0,     0,
       2,     3,     4,     0,     5,    68,     3,     4,    69,     5,
       2,     3,     4,    90,     5,     2,     3,     4,     0,     5
};

static const yytype_int8 yycheck[] =
{
       8,    30,    29,    42,    31,     0,     1,    20,    18,     5,
      19,    20,    18,    42,    11,    17,     5,     5,    57,    15,
       5,    16,    15,     5,     6,    20,    15,    23,    67,    31,
      38,     3,    71,     5,    23,    64,    29,    30,    31,    66,
      69,    51,    52,    53,    27,    54,    20,    19,     3,    21,
      22,     3,     6,    14,    62,     5,    28,    65,    97,     5,
      14,    90,    34,    35,    19,    19,    21,    22,     3,     5,
       1,     3,    15,    28,     3,     3,     1,    13,    31,    34,
      35,     5,     6,    91,    15,    16,    29,    30,    31,    20,
      15,    16,     6,     1,     6,    20,    27,     6,    29,    30,
      31,     5,    27,     6,    29,    30,    31,    15,    16,    13,
      15,    86,    20,     6,     7,     8,     9,    10,    39,    27,
      -1,    29,    30,    31,    29,    30,    31,    -1,    -1,    22,
      23,    24,    25,    26,     7,     8,     9,    10,    -1,    12,
      13,     7,     8,     9,    10,    -1,    12,    13,    -1,    22,
      23,    24,    25,    26,    15,    -1,    22,    23,    24,    25,
      26,     7,     8,     9,    10,    24,    12,    15,    29,    30,
      31,    15,    -1,    32,    33,    -1,    22,    23,    24,    25,
      26,    29,    30,    31,    14,    29,    30,    31,    -1,    -1,
      -1,    -1,    22,    23,    24,     1,    26,    -1,    -1,    -1,
       6,     7,     8,    -1,    10,     6,     7,     8,     9,    10,
       6,     7,     8,     9,    10,     6,     7,     8,    -1,    10
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     6,     7,     8,    10,    38,    39,    40,    41,
      42,    43,    20,     3,     5,    19,    21,    22,    28,    34,
      35,    47,    49,    50,    52,     0,     1,    16,    20,    18,
       9,    12,    22,    23,    24,    25,    26,    41,    44,    45,
      46,     6,    14,    19,    11,     5,     5,    49,    38,    38,
      15,    29,    30,    31,    27,     5,    15,    23,    20,     3,
      48,    51,    40,    52,     9,    40,    44,    14,     6,     9,
      41,    14,    46,     5,    47,    52,     3,    49,    49,    49,
      38,    47,    24,    32,    33,    17,    31,     6,    52,    13,
       9,    40,    47,    52,    47,    15,    15,    23,    15,     3,
       3,     3,    51,    13,     6,    52,    13,    15,     6,    15,
      47,     6,    13,     6,    15,     6
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (input, molList, doQueries, scanner, YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (&yylval, YYLEX_PARAM)
#else
# define YYLEX yylex (&yylval, scanner)
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value, input, molList, doQueries, scanner); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, doQueries, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    bool doQueries;
    void *scanner;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (input);
  YYUSE (molList);
  YYUSE (doQueries);
  YYUSE (scanner);
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, input, molList, doQueries, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    bool doQueries;
    void *scanner;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, doQueries, scanner);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
#else
static void
yy_reduce_print (yyvsp, yyrule, input, molList, doQueries, scanner)
    YYSTYPE *yyvsp;
    int yyrule;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    bool doQueries;
    void *scanner;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       , input, molList, doQueries, scanner);
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule, input, molList, doQueries, scanner); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, input, molList, doQueries, scanner)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    bool doQueries;
    void *scanner;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (doQueries);
  YYUSE (scanner);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */






/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
#else
int
yyparse (input, molList, doQueries, scanner)
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    bool doQueries;
    void *scanner;
#endif
#endif
{
  /* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;

  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
#line 134 "sln.yy"
    {
  // allowing mol<attrs><attrs> seems to be a NIBR thing, I don't
  // think it's standard SLN
  RDKit::ROMol *mol=(*molList)[(yyvsp[(1) - (4)].mol_T)];
  SLNParse::parseMolAttribs(mol,*(yyvsp[(3) - (4)].attriblist_T));
  delete (yyvsp[(3) - (4)].attriblist_T);
  (yyval.mol_T)=(yyvsp[(1) - (4)].mol_T);
;}
    break;

  case 4:
#line 142 "sln.yy"
    {
  yyclearin;
  yyerrok;
  for(std::vector<RDKit::RWMol *>::iterator iter=molList->begin();
      iter!=molList->end();++iter){
    SLNParse::CleanupAfterParseError(*iter);
    delete *iter;
  }
  molList->clear();
  molList->resize(0);
  YYABORT;
;}
    break;

  case 5:
#line 154 "sln.yy"
    {
  YYACCEPT;
;}
    break;

  case 6:
#line 157 "sln.yy"
    {
  yyclearin;
  yyerrok;
  for(std::vector<RDKit::RWMol *>::iterator iter=molList->begin();
      iter!=molList->end();++iter){
    SLNParse::CleanupAfterParseError(*iter);
    delete *iter;
  }
  molList->clear();
  molList->resize(0);
  YYABORT;
;}
    break;

  case 8:
#line 172 "sln.yy"
    {
  (yyval.mol_T)=SLNParse::addFragToMol(*molList,(yyvsp[(1) - (3)].mol_T),(yyvsp[(3) - (3)].mol_T));
;}
    break;

  case 9:
#line 177 "sln.yy"
    {
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  
  (yyval.mol_T)=SLNParse::startMol(*molList,newAtom,doQueries);
;}
    break;

  case 10:
#line 187 "sln.yy"
    {
  (yyval.mol_T)=SLNParse::startMol(*molList,(yyvsp[(1) - (1)].atom_T),doQueries);
;}
    break;

  case 11:
#line 190 "sln.yy"
    {
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),(yyvsp[(2) - (2)].atom_T),doQueries);
  (yyval.mol_T)=(yyvsp[(1) - (2)].mol_T);
;}
    break;

  case 12:
#line 194 "sln.yy"
    {
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),(yyvsp[(3) - (3)].atom_T),(yyvsp[(2) - (3)].bond_T),doQueries);
  (yyval.mol_T)=(yyvsp[(1) - (3)].mol_T);
;}
    break;

  case 13:
#line 198 "sln.yy"
    {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(3) - (3)].ival_T));
  (yyval.mol_T)=(yyvsp[(1) - (3)].mol_T);
;}
    break;

  case 14:
#line 202 "sln.yy"
    {
  // closeRingBond() takes ownership of the bond
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(4) - (4)].ival_T),(yyvsp[(2) - (4)].bond_T));
  (yyval.mol_T)=(yyvsp[(1) - (4)].mol_T);
;}
    break;

  case 15:
#line 207 "sln.yy"
    {
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[(3) - (4)].mol_T));
  (yyval.mol_T)=(yyvsp[(1) - (4)].mol_T);
;}
    break;

  case 16:
#line 211 "sln.yy"
    {
  // addBranchToMol() takes ownership of the bond and deletes the
  // branch, so no leaks here'
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[(4) - (5)].mol_T),(yyvsp[(3) - (5)].bond_T));
  (yyval.mol_T)=(yyvsp[(1) - (5)].mol_T);
;}
    break;

  case 17:
#line 217 "sln.yy"
    {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(4) - (5)].ival_T));
  (yyval.mol_T)=(yyvsp[(1) - (5)].mol_T);
;}
    break;

  case 18:
#line 221 "sln.yy"
    {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(5) - (6)].ival_T),(yyvsp[(3) - (6)].bond_T));
  (yyval.mol_T)=(yyvsp[(1) - (6)].mol_T);
;}
    break;

  case 19:
#line 225 "sln.yy"
    {
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,(yyvsp[(2) - (3)].bond_T),doQueries);
  (yyval.mol_T)=(yyvsp[(1) - (3)].mol_T);
;}
    break;

  case 20:
#line 236 "sln.yy"
    {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(3) - (4)].ival_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);

  (yyval.mol_T)=(yyvsp[(1) - (4)].mol_T);
;}
    break;

  case 21:
#line 248 "sln.yy"
    {
  // closeRingBond() takes ownership of the bond
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(4) - (5)].ival_T),(yyvsp[(2) - (5)].bond_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);
  
  (yyval.mol_T)=(yyvsp[(1) - (5)].mol_T);
;}
    break;

  case 22:
#line 261 "sln.yy"
    {
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[(3) - (5)].mol_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);
  
  (yyval.mol_T)=(yyvsp[(1) - (5)].mol_T);
;}
    break;

  case 23:
#line 273 "sln.yy"
    {
  // addBranchToMol() takes ownership of the bond and deletes the
  // branch, so no leaks here'
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[(4) - (6)].mol_T),(yyvsp[(3) - (6)].bond_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);
  
  (yyval.mol_T)=(yyvsp[(1) - (6)].mol_T);
;}
    break;

  case 24:
#line 287 "sln.yy"
    {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(4) - (6)].ival_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);
  
  (yyval.mol_T)=(yyvsp[(1) - (6)].mol_T);
;}
    break;

  case 25:
#line 299 "sln.yy"
    {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[(5) - (7)].ival_T),(yyvsp[(3) - (7)].bond_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);
  
  (yyval.mol_T)=(yyvsp[(1) - (7)].mol_T);
;}
    break;

  case 28:
#line 315 "sln.yy"
    {
  (yyvsp[(1) - (2)].atom_T)->setNumExplicitHs(1);
  (yyval.atom_T)=(yyvsp[(1) - (2)].atom_T);
;}
    break;

  case 29:
#line 319 "sln.yy"
    {
  (yyvsp[(1) - (3)].atom_T)->setNumExplicitHs((yyvsp[(3) - (3)].ival_T));
  (yyval.atom_T)=(yyvsp[(1) - (3)].atom_T);
;}
    break;

  case 30:
#line 325 "sln.yy"
    {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  (yyval.atom_T)->setProp("_starred",1,true);
;}
    break;

  case 31:
#line 333 "sln.yy"
    {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  (yyval.atom_T)->setProp("_AtomID",static_cast<unsigned int>((yyvsp[(2) - (3)].ival_T)));
;}
    break;

  case 32:
#line 341 "sln.yy"
    {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  (yyval.atom_T)->setProp("_AtomID",static_cast<unsigned int>((yyvsp[(2) - (5)].ival_T)));
  SLNParse::parseAtomAttribs((yyval.atom_T),*(yyvsp[(4) - (5)].attriblist_T),doQueries);
  delete (yyvsp[(4) - (5)].attriblist_T);
;}
    break;

  case 33:
#line 351 "sln.yy"
    {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  SLNParse::parseAtomAttribs((yyval.atom_T),*(yyvsp[(2) - (3)].attriblist_T),doQueries);
  delete (yyvsp[(2) - (3)].attriblist_T);
;}
    break;

  case 35:
#line 363 "sln.yy"
    {
  (yyval.atom_T)->setProp("_starred",1,true);
;}
    break;

  case 36:
#line 366 "sln.yy"
    {
  (yyvsp[(1) - (4)].atom_T)->setProp("_AtomID",static_cast<unsigned int>((yyvsp[(3) - (4)].ival_T)));
  (yyval.atom_T)=(yyvsp[(1) - (4)].atom_T);
;}
    break;

  case 37:
#line 370 "sln.yy"
    {
  (yyvsp[(1) - (6)].atom_T)->setProp("_AtomID",static_cast<unsigned int>((yyvsp[(3) - (6)].ival_T)));
  SLNParse::parseAtomAttribs((yyvsp[(1) - (6)].atom_T),*(yyvsp[(5) - (6)].attriblist_T),doQueries);
  delete (yyvsp[(5) - (6)].attriblist_T);
  (yyval.atom_T)=(yyvsp[(1) - (6)].atom_T);
;}
    break;

  case 38:
#line 376 "sln.yy"
    {
  SLNParse::parseAtomAttribs((yyvsp[(1) - (4)].atom_T),*(yyvsp[(3) - (4)].attriblist_T),doQueries);
  delete (yyvsp[(3) - (4)].attriblist_T);
  (yyval.atom_T)=(yyvsp[(1) - (4)].atom_T);
;}
    break;

  case 40:
#line 385 "sln.yy"
    {
  SLNParse::parseBondAttribs((yyvsp[(1) - (4)].bond_T),*(yyvsp[(3) - (4)].attriblist_T),doQueries);
  delete (yyvsp[(3) - (4)].attriblist_T);
  (yyval.bond_T) = (yyvsp[(1) - (4)].bond_T);
;}
    break;

  case 41:
#line 392 "sln.yy"
    { 
  RDKit::Bond *bond=new RDKit::QueryBond();
  bond->setQuery(RDKit::makeBondNullQuery());   
  (yyval.bond_T) = bond;
;}
    break;

  case 42:
#line 397 "sln.yy"
    {
  RDKit::Bond *bond=new RDKit::QueryBond();
  bond->setQuery(RDKit::makeBondNullQuery());   
  SLNParse::parseBondAttribs(bond,*(yyvsp[(3) - (4)].attriblist_T),doQueries);
  delete (yyvsp[(3) - (4)].attriblist_T);
  (yyval.bond_T) = bond;
;}
    break;

  case 44:
#line 407 "sln.yy"
    {
	if(!doQueries){
        yysln_error(input,molList,doQueries,0,"sequential bonds not allowed in non-queries");
    molList->clear();
    molList->resize(0);
    YYABORT;
	} else {
	  RDKit::QueryBond *b1=static_cast<RDKit::QueryBond *>((yyvsp[(1) - (2)].bond_T));
	  RDKit::QueryBond *b2=static_cast<RDKit::QueryBond *>((yyvsp[(2) - (2)].bond_T));
	  b1->expandQuery(b2->getQuery()->copy(),Queries::COMPOSITE_OR,true);
		delete b2;
	}
;}
    break;

  case 45:
#line 422 "sln.yy"
    {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::SINGLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::SINGLE);
  }
  (yyval.bond_T) = bond;
;}
    break;

  case 46:
#line 431 "sln.yy"
    {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::DOUBLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::DOUBLE);
  }
  (yyval.bond_T) = bond;
;}
    break;

  case 47:
#line 440 "sln.yy"
    {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::TRIPLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::TRIPLE);
  }
  (yyval.bond_T) = bond;

;}
    break;

  case 48:
#line 450 "sln.yy"
    {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::AROMATIC);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::AROMATIC);
  }
  (yyval.bond_T) = bond;
;}
    break;

  case 49:
#line 462 "sln.yy"
    {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[(3) - (3)].attrib_T))));
;}
    break;

  case 50:
#line 466 "sln.yy"
    {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribOr,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[(3) - (3)].attrib_T))));
;}
    break;

  case 51:
#line 470 "sln.yy"
    {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribLowPriAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[(3) - (3)].attrib_T))));
;}
    break;

  case 52:
#line 474 "sln.yy"
    {
  (yyval.attriblist_T) = new SLNParse::AttribListType();
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribLowPriAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[(1) - (1)].attrib_T))));
;}
    break;

  case 53:
#line 481 "sln.yy"
    {
  (yyval.attriblist_T) = new SLNParse::AttribListType();
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[(1) - (1)].attrib_T))));
;}
    break;

  case 54:
#line 486 "sln.yy"
    {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[(3) - (3)].attrib_T))));
;}
    break;

  case 55:
#line 492 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[(1) - (1)].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  (yyval.attrib_T)->op = "";
  (yyval.attrib_T)->second = "";
  delete (yyvsp[(1) - (1)].text_T);
;}
    break;

  case 56:
#line 500 "sln.yy"
    {
  (yyvsp[(2) - (2)].attrib_T)->negated=true;
  (yyval.attrib_T)=(yyvsp[(2) - (2)].attrib_T);
;}
    break;

  case 57:
#line 504 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[(1) - (3)].text_T);
  (yyval.attrib_T)->op = *(yyvsp[(2) - (3)].text_T);
  (yyval.attrib_T)->second = *(yyvsp[(3) - (3)].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[(1) - (3)].text_T);
  delete (yyvsp[(2) - (3)].text_T);
  delete (yyvsp[(3) - (3)].text_T);
;}
    break;

  case 58:
#line 515 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = "+1";
;}
    break;

  case 59:
#line 521 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = SLNParse::convertToString((yyvsp[(2) - (2)].ival_T));
;}
    break;

  case 60:
#line 527 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = "-1";
;}
    break;

  case 61:
#line 533 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = SLNParse::convertToString(-(yyvsp[(2) - (2)].ival_T));
;}
    break;

  case 62:
#line 539 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "spin";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = "d";
;}
    break;

  case 63:
#line 545 "sln.yy"
    {
  (yyval.attrib_T) = (yyvsp[(1) - (1)].attrib_T);
;}
    break;

  case 64:
#line 550 "sln.yy"
    {
   int sz = molList->size();
   RDKit::ROMol *mol=(*molList)[(yyvsp[(2) - (2)].mol_T)];
   molList->resize( sz-1 );
   SLNParse::finalizeQueryMol(mol,true);
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(mol);
   RDKit::ATOM_OR_QUERY *orq=new RDKit::ATOM_OR_QUERY();
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   (yyval.attrib_T) = new SLNParse::AttribType();
   (yyval.attrib_T)->first="is";
   (yyval.attrib_T)->op = "=";
   (yyval.attrib_T)->second = "";
   (yyval.attrib_T)->structQuery=static_cast<void *>(orq);
;}
    break;

  case 65:
#line 564 "sln.yy"
    {
   int sz = molList->size();
   RDKit::ROMol *mol=(*molList)[(yyvsp[(2) - (2)].mol_T)];
   molList->resize( sz-1 );
   SLNParse::finalizeQueryMol(mol,true);
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(mol);
   RDKit::ATOM_OR_QUERY *orq=new RDKit::ATOM_OR_QUERY();
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   orq->setNegation(true);

   (yyval.attrib_T) = new SLNParse::AttribType();
   (yyval.attrib_T)->first="is";
   (yyval.attrib_T)->op = "=";
   (yyval.attrib_T)->second = "";
   (yyval.attrib_T)->structQuery=static_cast<void *>(orq);
;}
    break;

  case 66:
#line 580 "sln.yy"
    {
   int sz = molList->size();
   RDKit::ROMol *mol=(*molList)[(yyvsp[(3) - (3)].mol_T)];
   molList->resize( sz-1 );
   SLNParse::finalizeQueryMol(mol,true);
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(mol);

   RDKit::ATOM_OR_QUERY *orq=static_cast<RDKit::ATOM_OR_QUERY *>((yyvsp[(1) - (3)].attrib_T)->structQuery);
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   (yyval.attrib_T)=(yyvsp[(1) - (3)].attrib_T);
;}
    break;

  case 67:
#line 594 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[(1) - (1)].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  (yyval.attrib_T)->op = "";
  (yyval.attrib_T)->second = "";
  delete (yyvsp[(1) - (1)].text_T);
;}
    break;

  case 68:
#line 602 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[(1) - (3)].text_T);
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = *(yyvsp[(3) - (3)].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[(1) - (3)].text_T);
  delete (yyvsp[(3) - (3)].text_T);
;}
    break;

  case 69:
#line 612 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[(1) - (3)].text_T);
  (yyval.attrib_T)->op = ":=";
  (yyval.attrib_T)->second = *(yyvsp[(3) - (3)].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[(1) - (3)].text_T);
  delete (yyvsp[(3) - (3)].text_T);
;}
    break;

  case 70:
#line 622 "sln.yy"
    {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[(1) - (3)].text_T);
  (yyval.attrib_T)->op = "^=";
  (yyval.attrib_T)->second = *(yyvsp[(3) - (3)].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[(1) - (3)].text_T);
  delete (yyvsp[(3) - (3)].text_T);
;}
    break;

  case 72:
#line 635 "sln.yy"
    { (yyval.ival_T) = (yyvsp[(1) - (2)].ival_T)*10 + (yyvsp[(2) - (2)].ival_T); ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2346 "/Users/landrgr1/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (input, molList, doQueries, scanner, YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (input, molList, doQueries, scanner, yymsg);
	  }
	else
	  {
	    yyerror (input, molList, doQueries, scanner, YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval, input, molList, doQueries, scanner);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp, input, molList, doQueries, scanner);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (input, molList, doQueries, scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval, input, molList, doQueries, scanner);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, input, molList, doQueries, scanner);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 637 "sln.yy"




