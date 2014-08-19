
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
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
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse         yysmarts_parse
#define yylex           yysmarts_lex
#define yyerror         yysmarts_error
#define yylval          yysmarts_lval
#define yychar          yysmarts_char
#define yydebug         yysmarts_debug
#define yynerrs         yysmarts_nerrs


/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 1 "smarts.yy"


  // $Id$
  //
  //  Copyright (C) 2003-2011 Greg Landrum and Rational Discovery LLC
  //
  //   @@ All Rights Reserved  @@
  //
#include <cstring>
#include <iostream>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>  
#include <GraphMol/SmilesParse/SmilesParseOps.h>  
#include <RDGeneral/RDLog.h>
#include "smarts.tab.hpp"

extern int yysmarts_lex(YYSTYPE *,void *);


#define YYDEBUG 1

void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
		void *scanner,const char * msg )
{
  throw RDKit::SmilesParseException(msg);
}

using namespace RDKit;
namespace {
 void yyErrorCleanup(std::vector<RDKit::RWMol *> *molList){
  for(std::vector<RDKit::RWMol *>::iterator iter=molList->begin();
      iter != molList->end(); ++iter){
     delete *iter;
  }
  molList->clear();
  molList->resize(0);
 }
}


/* Line 189 of yacc.c  */
#line 127 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
     IMPLICIT_H_ATOM_QUERY_TOKEN = 264,
     HYB_TOKEN = 265,
     ZERO_TOKEN = 266,
     NONZERO_DIGIT_TOKEN = 267,
     GROUP_OPEN_TOKEN = 268,
     GROUP_CLOSE_TOKEN = 269,
     SEPARATOR_TOKEN = 270,
     HASH_TOKEN = 271,
     MINUS_TOKEN = 272,
     PLUS_TOKEN = 273,
     CHIRAL_MARKER_TOKEN = 274,
     CHI_CLASS_TOKEN = 275,
     CHI_CLASS_OH_TOKEN = 276,
     H_TOKEN = 277,
     AT_TOKEN = 278,
     PERCENT_TOKEN = 279,
     ATOM_OPEN_TOKEN = 280,
     ATOM_CLOSE_TOKEN = 281,
     NOT_TOKEN = 282,
     AND_TOKEN = 283,
     OR_TOKEN = 284,
     SEMI_TOKEN = 285,
     BEGIN_RECURSE = 286,
     END_RECURSE = 287,
     COLON_TOKEN = 288,
     UNDERSCORE_TOKEN = 289,
     BOND_TOKEN = 290,
     EOS_TOKEN = 291
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 52 "smarts.yy"

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;



/* Line 214 of yacc.c  */
#line 208 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 220 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
# if YYENABLE_NLS
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
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
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
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
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
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  33
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   375

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  18
/* YYNRULES -- Number of rules.  */
#define YYNRULES  77
/* YYNRULES -- Number of states.  */
#define YYNSTATES  109

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
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     9,    12,    15,    17,    20,    24,
      27,    31,    34,    38,    42,    47,    49,    53,    59,    63,
      69,    73,    77,    81,    84,    86,    89,    91,    93,    97,
     103,   105,   108,   110,   113,   116,   120,   122,   125,   127,
     130,   132,   135,   138,   140,   142,   145,   147,   149,   151,
     153,   155,   157,   161,   165,   169,   171,   173,   176,   178,
     180,   182,   184,   186,   189,   192,   195,   197,   200,   203,
     205,   207,   211,   213,   215,   217,   220,   222
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      38,     0,    -1,    39,    -1,    38,     1,    36,    -1,    38,
      36,    -1,     1,    36,    -1,    41,    -1,    39,    41,    -1,
      39,    47,    41,    -1,    39,    51,    -1,    39,    47,    51,
      -1,    39,    40,    -1,    39,    15,    39,    -1,    13,    39,
      14,    -1,    13,    47,    39,    14,    -1,    46,    -1,    25,
      22,    26,    -1,    25,    22,    33,    52,    26,    -1,    25,
      42,    26,    -1,    25,    42,    33,    52,    26,    -1,    42,
      28,    42,    -1,    42,    29,    42,    -1,    42,    30,    42,
      -1,    42,    43,    -1,    43,    -1,    27,    43,    -1,    44,
      -1,    45,    -1,    31,    39,    32,    -1,    31,    39,    32,
      34,    53,    -1,    46,    -1,    52,    46,    -1,     5,    -1,
      52,     5,    -1,    16,    52,    -1,    52,    16,    52,    -1,
       7,    -1,     7,    52,    -1,     8,    -1,     8,    52,    -1,
       9,    -1,     9,    52,    -1,    22,    52,    -1,    22,    -1,
      50,    -1,    23,    23,    -1,    23,    -1,    10,    -1,    52,
      -1,     4,    -1,     3,    -1,     6,    -1,    47,    28,    47,
      -1,    47,    29,    47,    -1,    47,    30,    47,    -1,    48,
      -1,    49,    -1,    48,    49,    -1,    35,    -1,    17,    -1,
      16,    -1,    33,    -1,    23,    -1,    27,    49,    -1,    18,
      18,    -1,    18,    52,    -1,    18,    -1,    17,    17,    -1,
      17,    52,    -1,    17,    -1,    54,    -1,    24,    12,    54,
      -1,    11,    -1,    53,    -1,    12,    -1,    53,    54,    -1,
      12,    -1,    11,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    87,    87,    88,    94,    97,   108,   117,   146,   167,
     200,   218,   230,   243,   244,   255,   256,   260,   265,   269,
     277,   282,   287,   292,   297,   300,   304,   305,   309,   326,
     350,   351,   355,   356,   360,   361,   365,   366,   369,   370,
     374,   375,   379,   384,   389,   394,   400,   406,   407,   415,
     427,   432,   437,   441,   445,   449,   452,   453,   460,   461,
     467,   473,   479,   484,   491,   492,   493,   494,   495,   496,
     500,   501,   506,   507,   511,   512,   515,   516
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "AROMATIC_ATOM_TOKEN",
  "ORGANIC_ATOM_TOKEN", "ATOM_TOKEN", "SIMPLE_ATOM_QUERY_TOKEN",
  "COMPLEX_ATOM_QUERY_TOKEN", "RINGSIZE_ATOM_QUERY_TOKEN",
  "IMPLICIT_H_ATOM_QUERY_TOKEN", "HYB_TOKEN", "ZERO_TOKEN",
  "NONZERO_DIGIT_TOKEN", "GROUP_OPEN_TOKEN", "GROUP_CLOSE_TOKEN",
  "SEPARATOR_TOKEN", "HASH_TOKEN", "MINUS_TOKEN", "PLUS_TOKEN",
  "CHIRAL_MARKER_TOKEN", "CHI_CLASS_TOKEN", "CHI_CLASS_OH_TOKEN",
  "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN", "ATOM_OPEN_TOKEN",
  "ATOM_CLOSE_TOKEN", "NOT_TOKEN", "AND_TOKEN", "OR_TOKEN", "SEMI_TOKEN",
  "BEGIN_RECURSE", "END_RECURSE", "COLON_TOKEN", "UNDERSCORE_TOKEN",
  "BOND_TOKEN", "EOS_TOKEN", "$accept", "cmpd", "mol", "branch", "atomd",
  "atom_expr", "point_query", "recursive_query", "atom_query",
  "simple_atom", "bond_expr", "bond_query", "bondd", "charge_spec",
  "ring_number", "number", "nonzero_number", "digit", 0
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
       0,    37,    38,    38,    38,    38,    39,    39,    39,    39,
      39,    39,    39,    40,    40,    41,    41,    41,    41,    41,
      42,    42,    42,    42,    42,    43,    43,    43,    44,    44,
      45,    45,    45,    45,    45,    45,    45,    45,    45,    45,
      45,    45,    45,    45,    45,    45,    45,    45,    45,    46,
      46,    46,    47,    47,    47,    47,    48,    48,    49,    49,
      49,    49,    49,    49,    50,    50,    50,    50,    50,    50,
      51,    51,    52,    52,    53,    53,    54,    54
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     3,     2,     2,     1,     2,     3,     2,
       3,     2,     3,     3,     4,     1,     3,     5,     3,     5,
       3,     3,     3,     2,     1,     2,     1,     1,     3,     5,
       1,     2,     1,     2,     2,     3,     1,     2,     1,     2,
       1,     2,     2,     1,     1,     2,     1,     1,     1,     1,
       1,     1,     3,     3,     3,     1,     1,     2,     1,     1,
       1,     1,     1,     2,     2,     2,     1,     2,     2,     1,
       1,     3,     1,     1,     1,     2,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    50,    49,    51,     0,     0,     2,     6,    15,
       5,    32,    36,    38,    40,    47,    72,    74,     0,    69,
      66,    43,    46,     0,     0,     0,    24,    26,    27,    30,
      44,    48,    73,     1,     0,     4,    77,    76,     0,     0,
      60,    59,    62,     0,     0,    61,    58,    11,     7,     0,
      55,    56,     9,    70,    37,    39,    41,    34,    67,    68,
      64,    65,    16,     0,    42,    45,    43,    25,     0,    18,
       0,     0,     0,     0,    23,    33,     0,    31,    75,     3,
       0,     0,    12,     0,    63,     0,     0,     0,     8,    10,
      57,     0,    28,    20,    21,    22,     0,    35,    13,     0,
      71,    52,    53,    54,    17,     0,    19,    14,    29
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     6,     7,    47,     8,    25,    26,    27,    28,     9,
      49,    50,    51,    30,    52,    31,    32,    53
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -54
static const yytype_int16 yypact[] =
{
      40,   -26,   -54,   -54,   -54,   294,     4,   157,   -54,   -54,
     -54,   -54,    60,    60,    60,   -54,   -54,   -54,    60,    47,
      51,    24,   -17,   323,    45,   205,   -54,   -54,   -54,   -54,
     -54,   123,    71,   -54,   -20,   -54,   -54,   -54,   172,    45,
     -54,   -54,   -54,     8,   320,   -54,   -54,   -54,   -54,   345,
     320,   -54,   -54,   -54,   -54,   -54,   -54,   -54,   -54,   -54,
     -54,   -54,   -54,    60,   -54,   -54,    60,   -54,    82,   -54,
     323,   323,   323,    60,   -54,   -54,    60,   -54,   -54,   -54,
     107,     9,   157,    71,   -54,   320,   320,   320,   -54,   -54,
     -54,     0,   -12,   323,   265,   236,     2,   -54,   -54,   132,
     -54,   -54,     5,    75,   -54,    33,   -54,   -54,    71
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -54,   -54,   -21,   -54,    -7,    30,   -14,   -54,   -54,    -4,
     -31,   -54,   -36,   -54,    -2,    11,   -53,   -30
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      48,    29,    78,    68,    33,    34,    65,    81,    84,    67,
      10,    74,     2,     3,    90,     4,    79,    80,    82,    29,
      83,    29,   105,    54,    55,    56,   104,    77,   106,    57,
      59,    61,    64,    85,     5,    16,    17,    85,    86,    87,
      35,     1,    88,     2,     3,    17,     4,    89,     2,     3,
      62,     4,   108,   100,   101,   102,   103,    63,    16,    17,
      99,    48,    16,    17,    58,     5,    29,    29,    29,    60,
       5,    16,    17,    48,    91,    48,     0,    64,    78,    74,
      74,    74,    36,    37,    96,     2,     3,    97,     4,    29,
      29,    29,    48,    36,    37,    38,     0,    39,    40,    41,
      93,    94,    95,    85,    86,    42,    43,     5,     0,    44,
       2,     3,     0,     4,    92,    45,     0,    46,    36,    37,
      38,    98,    39,    40,    41,     0,     2,     3,    75,     4,
      42,    43,     5,     0,    44,     2,     3,     0,     4,    76,
      45,     0,    46,    36,    37,    38,   107,    39,    40,    41,
       0,     0,     0,     0,     0,    42,    43,     5,     0,    44,
       2,     3,     0,     4,     0,    45,     0,    46,    36,    37,
      38,     0,    39,    40,    41,     2,     3,     0,     4,     0,
      42,    43,     5,     0,    44,     0,     0,     0,    40,    41,
      45,     0,    46,     0,     0,    42,     0,     5,     0,    44,
       0,     0,     0,     0,     0,    45,     0,    46,     2,     3,
      11,     4,    12,    13,    14,    15,    16,    17,     0,     0,
       0,    18,    19,    20,     0,     0,     0,    66,    22,     0,
       0,    69,    23,    70,    71,    72,    24,     0,    73,     2,
       3,    11,     4,    12,    13,    14,    15,    16,    17,     0,
       0,     0,    18,    19,    20,     0,     0,     0,    66,    22,
       0,     0,     0,    23,    70,    71,     0,    24,     2,     3,
      11,     4,    12,    13,    14,    15,    16,    17,     0,     0,
       0,    18,    19,    20,     0,     0,     0,    66,    22,     0,
       0,     0,    23,    70,     0,     0,    24,     2,     3,    11,
       4,    12,    13,    14,    15,    16,    17,     0,     0,     0,
      18,    19,    20,     0,     0,     0,    21,    22,     0,     0,
       0,    23,     0,     0,     0,    24,     2,     3,    11,     4,
      12,    13,    14,    15,    16,    17,    40,    41,     0,    18,
      19,    20,     0,    42,     0,    66,    22,    44,     2,     3,
      23,     4,     0,    45,    24,    46,    36,    37,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    43,
       5,     0,     0,    85,    86,    87
};

static const yytype_int8 yycheck[] =
{
       7,     5,    32,    24,     0,     1,    23,    38,    44,    23,
      36,    25,     3,     4,    50,     6,    36,    38,    39,    23,
      12,    25,    34,    12,    13,    14,    26,    31,    26,    18,
      19,    20,    21,    28,    25,    11,    12,    28,    29,    30,
      36,     1,    49,     3,     4,    12,     6,    49,     3,     4,
      26,     6,   105,    83,    85,    86,    87,    33,    11,    12,
      81,    68,    11,    12,    17,    25,    70,    71,    72,    18,
      25,    11,    12,    80,    63,    82,    -1,    66,   108,    93,
      94,    95,    11,    12,    73,     3,     4,    76,     6,    93,
      94,    95,    99,    11,    12,    13,    -1,    15,    16,    17,
      70,    71,    72,    28,    29,    23,    24,    25,    -1,    27,
       3,     4,    -1,     6,    32,    33,    -1,    35,    11,    12,
      13,    14,    15,    16,    17,    -1,     3,     4,     5,     6,
      23,    24,    25,    -1,    27,     3,     4,    -1,     6,    16,
      33,    -1,    35,    11,    12,    13,    14,    15,    16,    17,
      -1,    -1,    -1,    -1,    -1,    23,    24,    25,    -1,    27,
       3,     4,    -1,     6,    -1,    33,    -1,    35,    11,    12,
      13,    -1,    15,    16,    17,     3,     4,    -1,     6,    -1,
      23,    24,    25,    -1,    27,    -1,    -1,    -1,    16,    17,
      33,    -1,    35,    -1,    -1,    23,    -1,    25,    -1,    27,
      -1,    -1,    -1,    -1,    -1,    33,    -1,    35,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    -1,    -1,
      -1,    16,    17,    18,    -1,    -1,    -1,    22,    23,    -1,
      -1,    26,    27,    28,    29,    30,    31,    -1,    33,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    -1,
      -1,    -1,    16,    17,    18,    -1,    -1,    -1,    22,    23,
      -1,    -1,    -1,    27,    28,    29,    -1,    31,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    -1,    -1,
      -1,    16,    17,    18,    -1,    -1,    -1,    22,    23,    -1,
      -1,    -1,    27,    28,    -1,    -1,    31,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    -1,    -1,    -1,
      16,    17,    18,    -1,    -1,    -1,    22,    23,    -1,    -1,
      -1,    27,    -1,    -1,    -1,    31,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    16,    17,    -1,    16,
      17,    18,    -1,    23,    -1,    22,    23,    27,     3,     4,
      27,     6,    -1,    33,    31,    35,    11,    12,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    24,
      25,    -1,    -1,    28,    29,    30
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     6,    25,    38,    39,    41,    46,
      36,     5,     7,     8,     9,    10,    11,    12,    16,    17,
      18,    22,    23,    27,    31,    42,    43,    44,    45,    46,
      50,    52,    53,     0,     1,    36,    11,    12,    13,    15,
      16,    17,    23,    24,    27,    33,    35,    40,    41,    47,
      48,    49,    51,    54,    52,    52,    52,    52,    17,    52,
      18,    52,    26,    33,    52,    23,    22,    43,    39,    26,
      28,    29,    30,    33,    43,     5,    16,    46,    54,    36,
      39,    47,    39,    12,    49,    28,    29,    30,    41,    51,
      49,    52,    32,    42,    42,    42,    52,    52,    14,    39,
      54,    47,    47,    47,    26,    34,    26,    14,    53
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
      yyerror (input, molList, scanner, YY_("syntax error: cannot back up")); \
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
# if YYLTYPE_IS_TRIVIAL
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
		  Type, Value, input, molList, scanner); \
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
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    void *scanner;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (input);
  YYUSE (molList);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, input, molList, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    void *scanner;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, scanner);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
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
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
#else
static void
yy_reduce_print (yyvsp, yyrule, input, molList, scanner)
    YYSTYPE *yyvsp;
    int yyrule;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
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
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       , input, molList, scanner);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule, input, molList, scanner); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, input, molList, scanner)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    void *scanner;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
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
int yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */





/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

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
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
#else
int
yyparse (input, molList, scanner)
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    void *scanner;
#endif
#endif
{
/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

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
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
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

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
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

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
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

/* Line 1455 of yacc.c  */
#line 88 "smarts.yy"
    {
  yyclearin;
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
;}
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 94 "smarts.yy"
    {
  YYACCEPT;
;}
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 97 "smarts.yy"
    {
  yyclearin;
  yyerrok;
  
  yyErrorCleanup(molList);
  YYABORT;
;}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 108 "smarts.yy"
    {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (yyvsp[(1) - (1)].atom)->setProp("_SmilesStart",1);
  (*molList)[ sz ]->addAtom((yyvsp[(1) - (1)].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 117 "smarts.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[(2) - (2)].atom),true,true);

  QueryBond *newB;
  // this is a bit of a hack to try and get nicer "SMILES" from 
  // a SMARTS molecule:
  if(!(a1->getIsAromatic() && (yyvsp[(2) - (2)].atom)->getIsAromatic())){
    newB = new QueryBond(Bond::SINGLE);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
  			    Queries::COMPOSITE_OR,
  			    true);
  } else {
    newB = new QueryBond(Bond::AROMATIC);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::SINGLE),
  			    Queries::COMPOSITE_OR,
  			    true);
  }
  newB->setProp("_unspecifiedOrder",1);
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  mp->addBond(newB);
  delete newB;
  //delete $2;
;}
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 146 "smarts.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[(3) - (3)].atom),true,true);
  if( (yyvsp[(2) - (3)].bond)->getBondType() == Bond::DATIVER ){
    (yyvsp[(2) - (3)].bond)->setBeginAtomIdx(atomIdx1);
    (yyvsp[(2) - (3)].bond)->setEndAtomIdx(atomIdx2);
    (yyvsp[(2) - (3)].bond)->setBondType(Bond::DATIVE);
  }else if ( (yyvsp[(2) - (3)].bond)->getBondType() == Bond::DATIVEL ){
    (yyvsp[(2) - (3)].bond)->setBeginAtomIdx(atomIdx2);
    (yyvsp[(2) - (3)].bond)->setEndAtomIdx(atomIdx1);
    (yyvsp[(2) - (3)].bond)->setBondType(Bond::DATIVE);
  } else {
    (yyvsp[(2) - (3)].bond)->setBeginAtomIdx(atomIdx1);
    (yyvsp[(2) - (3)].bond)->setEndAtomIdx(atomIdx2);
  }
  mp->addBond((yyvsp[(2) - (3)].bond));
  delete (yyvsp[(2) - (3)].bond);

;}
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 167 "smarts.yy"
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();

  // this is a bit of a hack to try and get nicer "SMILES" from 
  // a SMARTS molecule:
  QueryBond * newB;
  if(!atom->getIsAromatic()){
    newB = new QueryBond(Bond::SINGLE);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::AROMATIC),
  			    Queries::COMPOSITE_OR,
  			    true);
  } else {
    newB = new QueryBond(Bond::AROMATIC);
    newB->expandQuery(makeBondOrderEqualsQuery(Bond::SINGLE),
  			    Queries::COMPOSITE_OR,
  			    true);
  }
  newB->setProp("_unspecifiedOrder",1);
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atom->getIdx());
  mp->setBondBookmark(newB,(yyvsp[(2) - (2)].ival));

  mp->setAtomBookmark(atom,(yyvsp[(2) - (2)].ival));
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-((yyvsp[(2) - (2)].ival)+1));
  atom->setProp("_RingClosures",tmp);

;}
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 200 "smarts.yy"
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();

  mp->setBondBookmark((yyvsp[(2) - (3)].bond),(yyvsp[(3) - (3)].ival));
  (yyvsp[(2) - (3)].bond)->setOwningMol(mp);
  (yyvsp[(2) - (3)].bond)->setBeginAtomIdx(atom->getIdx());

  mp->setAtomBookmark(atom,(yyvsp[(3) - (3)].ival));
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-((yyvsp[(3) - (3)].ival)+1));
  atom->setProp("_RingClosures",tmp);

;}
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 218 "smarts.yy"
    {
  RWMol *m1_p = (*molList)[(yyval.moli)],*m2_p=(*molList)[(yyvsp[(2) - (2)].moli)];
  m2_p->getAtomWithIdx(0)->setProp("_SmilesStart",1);
  // FIX: handle generic bonds here
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED,Bond::NONE,false,true);
  delete m2_p;
  int sz = molList->size();
  if ( sz==(yyvsp[(2) - (2)].moli)+1) {
    molList->resize( sz-1 );
  }
;}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 230 "smarts.yy"
    {
  RWMol *m1_p = (*molList)[(yyvsp[(1) - (3)].moli)],*m2_p=(*molList)[(yyvsp[(3) - (3)].moli)];
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::IONIC,Bond::NONE,true);
  delete m2_p;
  int sz = molList->size();
  if ( sz==(yyvsp[(3) - (3)].moli)+1) {
    molList->resize( sz-1 );
  }
;}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 243 "smarts.yy"
    { (yyval.moli) = (yyvsp[(2) - (3)].moli); ;}
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 244 "smarts.yy"
    {
  // FIX: this needs to handle arbitrary bond_exprs
  (yyval.moli) = (yyvsp[(3) - (4)].moli);
  int sz     = molList->size();
  (yyvsp[(2) - (4)].bond)->setOwningMol((*molList)[ sz-1 ]);
  (yyvsp[(2) - (4)].bond)->setBeginAtomIdx(0);
  (*molList)[ sz-1 ]->setBondBookmark((yyvsp[(2) - (4)].bond),ci_LEADING_BOND);
;}
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 257 "smarts.yy"
    {
  (yyval.atom) = new QueryAtom(1);
;}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 261 "smarts.yy"
    {
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp("molAtomMapNumber",(yyvsp[(4) - (5)].ival));
;}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 266 "smarts.yy"
    {
  (yyval.atom) = (yyvsp[(2) - (3)].atom);
;}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 270 "smarts.yy"
    {
  (yyval.atom) = (yyvsp[(2) - (5)].atom);
  (yyval.atom)->setProp("molAtomMapNumber",(yyvsp[(4) - (5)].ival));
;}
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 277 "smarts.yy"
    {
  (yyvsp[(1) - (3)].atom)->expandQuery((yyvsp[(3) - (3)].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[(1) - (3)].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[(1) - (3)].atom)->setChiralTag((yyvsp[(3) - (3)].atom)->getChiralTag());
  delete (yyvsp[(3) - (3)].atom);
;}
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 282 "smarts.yy"
    {
  (yyvsp[(1) - (3)].atom)->expandQuery((yyvsp[(3) - (3)].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[(1) - (3)].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[(1) - (3)].atom)->setChiralTag((yyvsp[(3) - (3)].atom)->getChiralTag());
  delete (yyvsp[(3) - (3)].atom);
;}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 287 "smarts.yy"
    {
  (yyvsp[(1) - (3)].atom)->expandQuery((yyvsp[(3) - (3)].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[(1) - (3)].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[(1) - (3)].atom)->setChiralTag((yyvsp[(3) - (3)].atom)->getChiralTag());
  delete (yyvsp[(3) - (3)].atom);
;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 292 "smarts.yy"
    {
  (yyvsp[(1) - (2)].atom)->expandQuery((yyvsp[(2) - (2)].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[(1) - (2)].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[(1) - (2)].atom)->setChiralTag((yyvsp[(2) - (2)].atom)->getChiralTag());
  delete (yyvsp[(2) - (2)].atom);
;}
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 300 "smarts.yy"
    {
  (yyvsp[(2) - (2)].atom)->getQuery()->setNegation(!((yyvsp[(2) - (2)].atom)->getQuery()->getNegation()));
  (yyval.atom) = (yyvsp[(2) - (2)].atom);
;}
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 309 "smarts.yy"
    {
  // this is a recursive SMARTS expression
  QueryAtom *qA = new QueryAtom();
  //  FIX: there's maybe a leak here
  RWMol *molP = (*molList)[(yyvsp[(2) - (3)].moli)];
  // close any rings in the molecule:
  SmilesParseOps::CloseMolRings(molP,0);

  //molP->debugMol(std::cout);
  qA->setQuery(new RecursiveStructureQuery(molP));
  //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
  int sz = molList->size();
  if ( sz==(yyvsp[(2) - (3)].moli)+1) {
    molList->resize( sz-1 );
  }
  (yyval.atom) = qA;
;}
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 326 "smarts.yy"
    {
  // UNDOCUMENTED EXTENSION:
  // this is a recursive SMARTS expression with a serial number
  // please don't write your own SMARTS that include this extension:
  // the RDKit smarts parsing code will automatically insert serial
  // numbers for recursive smarts patterns.
  QueryAtom *qA = new QueryAtom();
  //  FIX: there's maybe a leak here
  RWMol *molP = (*molList)[(yyvsp[(2) - (5)].moli)];
  // close any rings in the molecule:
  SmilesParseOps::CloseMolRings(molP,0);

  //molP->debugMol(std::cout);
  qA->setQuery(new RecursiveStructureQuery(molP,(yyvsp[(5) - (5)].ival)));
  //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
  int sz = molList->size();
  if ( sz==(yyvsp[(2) - (5)].moli)+1) {
    molList->resize( sz-1 );
  }
  (yyval.atom) = qA;
;}
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 351 "smarts.yy"
    {
  (yyvsp[(2) - (2)].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[(1) - (2)].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[(2) - (2)].atom);
;}
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 356 "smarts.yy"
    {
  (yyvsp[(2) - (2)].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[(1) - (2)].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[(2) - (2)].atom);
;}
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 360 "smarts.yy"
    { (yyval.atom) = new QueryAtom((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 361 "smarts.yy"
    {
  (yyval.atom) = new QueryAtom((yyvsp[(3) - (3)].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[(1) - (3)].ival)),Queries::COMPOSITE_AND,true);
;}
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 366 "smarts.yy"
    {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[(1) - (2)].atom)->getQuery())->setVal((yyvsp[(2) - (2)].ival));
;}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 370 "smarts.yy"
    {
  delete (yyvsp[(1) - (2)].atom)->getQuery();
  (yyvsp[(1) - (2)].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[(2) - (2)].ival)));
;}
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 375 "smarts.yy"
    {
  delete (yyvsp[(1) - (2)].atom)->getQuery();
  (yyvsp[(1) - (2)].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[(2) - (2)].ival)));
;}
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 379 "smarts.yy"
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[(2) - (2)].ival)));
  (yyval.atom)=newQ;
;}
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 384 "smarts.yy"
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  (yyval.atom)=newQ;
;}
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 389 "smarts.yy"
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[(1) - (1)].ival)));
  (yyval.atom)=newQ;
;}
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 394 "smarts.yy"
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
;}
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 400 "smarts.yy"
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
;}
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 407 "smarts.yy"
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[(1) - (1)].ival)));
  (yyval.atom)=newQ;
;}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 415 "smarts.yy"
    {
  //
  // This construction (and some others) may seem odd, but the
  // SMARTS definition requires that an atom which is aliphatic on
  // input (i.e. something in the "organic subset" that is given with
  // a capital letter) only match aliphatic atoms.
  //
  // The following rule applies a similar logic to aromatic atoms.
  //
  (yyval.atom) = new QueryAtom((yyvsp[(1) - (1)].ival));
  (yyval.atom)->expandQuery(makeAtomAliphaticQuery(),Queries::COMPOSITE_AND);
;}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 427 "smarts.yy"
    {
  (yyval.atom) = new QueryAtom((yyvsp[(1) - (1)].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->expandQuery(makeAtomAromaticQuery(),Queries::COMPOSITE_AND);
;}
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 437 "smarts.yy"
    {
  (yyvsp[(1) - (3)].bond)->expandQuery((yyvsp[(3) - (3)].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[(3) - (3)].bond);
;}
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 441 "smarts.yy"
    {
  (yyvsp[(1) - (3)].bond)->expandQuery((yyvsp[(3) - (3)].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[(3) - (3)].bond);
;}
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 445 "smarts.yy"
    {
  (yyvsp[(1) - (3)].bond)->expandQuery((yyvsp[(3) - (3)].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[(3) - (3)].bond);
;}
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 453 "smarts.yy"
    {
  (yyvsp[(1) - (2)].bond)->expandQuery((yyvsp[(2) - (2)].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[(2) - (2)].bond);
;}
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 461 "smarts.yy"
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
;}
    break;

  case 60:

/* Line 1455 of yacc.c  */
#line 467 "smarts.yy"
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
;}
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 473 "smarts.yy"
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
;}
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 479 "smarts.yy"
    {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
;}
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 484 "smarts.yy"
    {
  (yyvsp[(2) - (2)].bond)->getQuery()->setNegation(!((yyvsp[(2) - (2)].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[(2) - (2)].bond);
;}
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 491 "smarts.yy"
    { (yyval.ival)=2; ;}
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 492 "smarts.yy"
    { (yyval.ival)=(yyvsp[(2) - (2)].ival); ;}
    break;

  case 66:

/* Line 1455 of yacc.c  */
#line 493 "smarts.yy"
    { (yyval.ival)=1; ;}
    break;

  case 67:

/* Line 1455 of yacc.c  */
#line 494 "smarts.yy"
    { (yyval.ival)=-2; ;}
    break;

  case 68:

/* Line 1455 of yacc.c  */
#line 495 "smarts.yy"
    { (yyval.ival)=-(yyvsp[(2) - (2)].ival); ;}
    break;

  case 69:

/* Line 1455 of yacc.c  */
#line 496 "smarts.yy"
    { (yyval.ival)=-1; ;}
    break;

  case 71:

/* Line 1455 of yacc.c  */
#line 501 "smarts.yy"
    { (yyval.ival) = (yyvsp[(2) - (3)].ival)*10+(yyvsp[(3) - (3)].ival); ;}
    break;

  case 75:

/* Line 1455 of yacc.c  */
#line 512 "smarts.yy"
    { (yyval.ival) = (yyvsp[(1) - (2)].ival)*10 + (yyvsp[(2) - (2)].ival); ;}
    break;



/* Line 1455 of yacc.c  */
#line 2272 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
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
      yyerror (input, molList, scanner, YY_("syntax error"));
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
	    yyerror (input, molList, scanner, yymsg);
	  }
	else
	  {
	    yyerror (input, molList, scanner, YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
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
		      yytoken, &yylval, input, molList, scanner);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
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
		  yystos[yystate], yyvsp, input, molList, scanner);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

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

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (input, molList, scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval, input, molList, scanner);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, input, molList, scanner);
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



/* Line 1675 of yacc.c  */
#line 519 "smarts.yy"




