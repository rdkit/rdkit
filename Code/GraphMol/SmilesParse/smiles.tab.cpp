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
#define yyparse yysmiles_parse
#define yylex   yysmiles_lex
#define yyerror yysmiles_error
#define yylval  yysmiles_lval
#define yychar  yysmiles_char
#define yydebug yysmiles_debug
#define yynerrs yysmiles_nerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     AROMATIC_ATOM_TOKEN = 258,
     ATOM_TOKEN = 259,
     ORGANIC_ATOM_TOKEN = 260,
     NONZERO_DIGIT_TOKEN = 261,
     ZERO_TOKEN = 262,
     GROUP_OPEN_TOKEN = 263,
     GROUP_CLOSE_TOKEN = 264,
     SEPARATOR_TOKEN = 265,
     LOOP_CONNECTOR_TOKEN = 266,
     MINUS_TOKEN = 267,
     PLUS_TOKEN = 268,
     CHIRAL_MARKER_TOKEN = 269,
     CHI_CLASS_TOKEN = 270,
     CHI_CLASS_OH_TOKEN = 271,
     H_TOKEN = 272,
     AT_TOKEN = 273,
     PERCENT_TOKEN = 274,
     COLON_TOKEN = 275,
     BOND_TOKEN = 276,
     ATOM_OPEN_TOKEN = 277,
     ATOM_CLOSE_TOKEN = 278,
     EOS_TOKEN = 279
   };
#endif
/* Tokens.  */
#define AROMATIC_ATOM_TOKEN 258
#define ATOM_TOKEN 259
#define ORGANIC_ATOM_TOKEN 260
#define NONZERO_DIGIT_TOKEN 261
#define ZERO_TOKEN 262
#define GROUP_OPEN_TOKEN 263
#define GROUP_CLOSE_TOKEN 264
#define SEPARATOR_TOKEN 265
#define LOOP_CONNECTOR_TOKEN 266
#define MINUS_TOKEN 267
#define PLUS_TOKEN 268
#define CHIRAL_MARKER_TOKEN 269
#define CHI_CLASS_TOKEN 270
#define CHI_CLASS_OH_TOKEN 271
#define H_TOKEN 272
#define AT_TOKEN 273
#define PERCENT_TOKEN 274
#define COLON_TOKEN 275
#define BOND_TOKEN 276
#define ATOM_OPEN_TOKEN 277
#define ATOM_CLOSE_TOKEN 278
#define EOS_TOKEN 279




/* Copy the first part of user declarations.  */
#line 3 "smiles.yy"


  // $Id: smiles.yy 2159 2012-08-28 04:01:49Z glandrum $
  //
  //  Copyright (C) 2001-2010 Randal Henne, Greg Landrum and Rational Discovery LLC
  //
  //   @@ All Rights Reserved  @@
  //

#include <cstring>
#include <iostream>
#include <vector>
#include <list>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>  
#include <GraphMol/SmilesParse/SmilesParseOps.h>  
#include <RDGeneral/RDLog.h>
#include "smiles.tab.hpp"

extern int yysmiles_lex(YYSTYPE *,void *);

#define YYDEBUG 1
#define YYLEX_PARAM scanner

void
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::list<unsigned int> *branchPoints,
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
#line 56 "smiles.yy"
{
  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  int                      ival;
}
/* Line 193 of yacc.c.  */
#line 207 "/Users/landrgr1/RDKit_trunk/Code/GraphMol/SmilesParse/smiles.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 220 "/Users/landrgr1/RDKit_trunk/Code/GraphMol/SmilesParse/smiles.tab.cpp"

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
#define YYFINAL  21
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   86

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  25
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  13
/* YYNRULES -- Number of rules.  */
#define YYNRULES  53
/* YYNRULES -- Number of states.  */
#define YYNSTATES  70

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   279

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
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     5,     9,    13,    16,    19,    21,    24,
      28,    32,    35,    39,    43,    47,    52,    57,    60,    62,
      68,    72,    74,    77,    81,    85,    88,    92,    96,    98,
     101,   104,   108,   112,   117,   119,   122,   126,   128,   131,
     135,   137,   140,   142,   145,   147,   149,   151,   155,   157,
     159,   161,   164,   166
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      26,     0,    -1,    27,    -1,    26,    10,    27,    -1,    26,
       1,    24,    -1,    26,    24,    -1,     1,    24,    -1,    28,
      -1,    27,    28,    -1,    27,    21,    28,    -1,    27,    12,
      28,    -1,    27,    34,    -1,    27,    21,    34,    -1,    27,
      12,    34,    -1,    27,     8,    28,    -1,    27,     8,    21,
      28,    -1,    27,     8,    12,    28,    -1,    27,     9,    -1,
      33,    -1,    22,    29,    20,    35,    23,    -1,    22,    29,
      23,    -1,    30,    -1,    30,    13,    -1,    30,    13,    13,
      -1,    30,    13,    35,    -1,    30,    12,    -1,    30,    12,
      12,    -1,    30,    12,    35,    -1,    17,    -1,    35,    17,
      -1,    17,    17,    -1,    35,    17,    17,    -1,    17,    17,
      35,    -1,    35,    17,    17,    35,    -1,    31,    -1,    31,
      17,    -1,    31,    17,    35,    -1,    32,    -1,    32,    18,
      -1,    32,    18,    18,    -1,    33,    -1,    35,    33,    -1,
       4,    -1,    35,     4,    -1,     5,    -1,     3,    -1,    37,
      -1,    19,     6,    37,    -1,     7,    -1,    36,    -1,     6,
      -1,    36,    37,    -1,     6,    -1,     7,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    78,    78,    79,    88,    94,    97,   107,   118,   128,
     148,   156,   173,   190,   205,   215,   235,   243,   253,   254,
     261,   269,   270,   271,   272,   273,   274,   275,   279,   280,
     281,   282,   283,   284,   285,   286,   287,   291,   292,   293,
     297,   298,   299,   300,   304,   305,   309,   310,   314,   315,
     319,   320,   323,   324
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "AROMATIC_ATOM_TOKEN", "ATOM_TOKEN",
  "ORGANIC_ATOM_TOKEN", "NONZERO_DIGIT_TOKEN", "ZERO_TOKEN",
  "GROUP_OPEN_TOKEN", "GROUP_CLOSE_TOKEN", "SEPARATOR_TOKEN",
  "LOOP_CONNECTOR_TOKEN", "MINUS_TOKEN", "PLUS_TOKEN",
  "CHIRAL_MARKER_TOKEN", "CHI_CLASS_TOKEN", "CHI_CLASS_OH_TOKEN",
  "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN", "COLON_TOKEN", "BOND_TOKEN",
  "ATOM_OPEN_TOKEN", "ATOM_CLOSE_TOKEN", "EOS_TOKEN", "$accept", "cmpd",
  "mol", "atomd", "charge_element", "h_element", "chiral_element",
  "element", "simple_atom", "ring_number", "number", "nonzero_number",
  "digit", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    25,    26,    26,    26,    26,    26,    27,    27,    27,
      27,    27,    27,    27,    27,    27,    27,    27,    28,    28,
      28,    29,    29,    29,    29,    29,    29,    29,    30,    30,
      30,    30,    30,    30,    30,    30,    30,    31,    31,    31,
      32,    32,    32,    32,    33,    33,    34,    34,    35,    35,
      36,    36,    37,    37
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     3,     3,     2,     2,     1,     2,     3,
       3,     2,     3,     3,     3,     4,     4,     2,     1,     5,
       3,     1,     2,     3,     3,     2,     3,     3,     1,     2,
       2,     3,     3,     4,     1,     2,     3,     1,     2,     3,
       1,     2,     1,     2,     1,     1,     1,     3,     1,     1,
       1,     2,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    45,    44,     0,     0,     2,     7,    18,     6,
      42,    50,    48,    28,     0,    21,    34,    37,    40,     0,
      49,     1,     0,     0,     5,    52,    53,     0,    17,     0,
       0,     0,     8,    11,    46,    30,     0,    20,    25,    22,
      35,    38,    43,    29,    41,    51,     4,     3,     0,     0,
      14,    10,    13,     0,     9,    12,    32,     0,    26,    27,
      23,    24,    36,    39,    31,    16,    15,    47,    19,    33
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     5,     6,     7,    14,    15,    16,    17,     8,    33,
      19,    20,    34
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -34
static const yytype_int8 yypact[] =
{
      11,   -20,   -34,   -34,    62,     8,    41,   -34,   -34,   -34,
     -34,   -34,   -34,    -2,   -10,    15,    28,     6,   -34,    67,
      23,   -34,    16,    14,   -34,   -34,   -34,    52,   -34,    32,
      46,    32,   -34,   -34,   -34,    74,    74,   -34,    71,    69,
      74,    38,   -34,    42,   -34,   -34,   -34,    41,    14,    14,
     -34,   -34,   -34,    23,   -34,   -34,   -34,    35,   -34,   -34,
     -34,   -34,   -34,   -34,    74,   -34,   -34,   -34,   -34,   -34
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -34,   -34,    63,    -6,   -34,   -34,   -34,   -34,     7,    -9,
     -33,   -34,   -19
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      32,    45,    56,    57,     9,    59,    61,    62,    21,    22,
      36,    18,     1,    37,     2,    35,     3,     2,    23,     3,
      52,    50,    55,    51,    41,    54,    44,    38,    39,    25,
      26,    69,    24,     4,    67,     2,     4,     3,    25,    26,
      46,    32,    65,    66,     2,    40,     3,    25,    26,    27,
      28,    30,    53,    29,     4,     2,    63,     3,    68,    64,
      30,     0,    31,     4,    48,     2,    10,     3,    11,    12,
       2,    42,     3,    49,     4,    11,    12,    11,    12,    13,
      11,    12,    60,    58,    43,     0,    47
};

static const yytype_int8 yycheck[] =
{
       6,    20,    35,    36,    24,    38,    39,    40,     0,     1,
      20,     4,     1,    23,     3,    17,     5,     3,    10,     5,
      29,    27,    31,    29,    18,    31,    19,    12,    13,     6,
       7,    64,    24,    22,    53,     3,    22,     5,     6,     7,
      24,    47,    48,    49,     3,    17,     5,     6,     7,     8,
       9,    19,     6,    12,    22,     3,    18,     5,    23,    17,
      19,    -1,    21,    22,    12,     3,     4,     5,     6,     7,
       3,     4,     5,    21,    22,     6,     7,     6,     7,    17,
       6,     7,    13,    12,    17,    -1,    23
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     5,    22,    26,    27,    28,    33,    24,
       4,     6,     7,    17,    29,    30,    31,    32,    33,    35,
      36,     0,     1,    10,    24,     6,     7,     8,     9,    12,
      19,    21,    28,    34,    37,    17,    20,    23,    12,    13,
      17,    18,     4,    17,    33,    37,    24,    27,    12,    21,
      28,    28,    34,     6,    28,    34,    35,    35,    12,    35,
      13,    35,    35,    18,    17,    28,    28,    37,    23,    35
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
      yyerror (input, molList, branchPoints, scanner, YY_("syntax error: cannot back up")); \
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
# define YYLEX yylex (&yylval)
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
		  Type, Value, input, molList, branchPoints, scanner); \
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
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, std::list<unsigned int> *branchPoints, void *scanner)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, branchPoints, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    std::list<unsigned int> *branchPoints;
    void *scanner;
#endif
{
  if (!yyvaluep)
    return;
  YYUSE (input);
  YYUSE (molList);
  YYUSE (branchPoints);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, std::list<unsigned int> *branchPoints, void *scanner)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep, input, molList, branchPoints, scanner)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    std::list<unsigned int> *branchPoints;
    void *scanner;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, branchPoints, scanner);
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
yy_reduce_print (YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, std::list<unsigned int> *branchPoints, void *scanner)
#else
static void
yy_reduce_print (yyvsp, yyrule, input, molList, branchPoints, scanner)
    YYSTYPE *yyvsp;
    int yyrule;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    std::list<unsigned int> *branchPoints;
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
		       		       , input, molList, branchPoints, scanner);
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule, input, molList, branchPoints, scanner); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, std::list<unsigned int> *branchPoints, void *scanner)
#else
static void
yydestruct (yymsg, yytype, yyvaluep, input, molList, branchPoints, scanner)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    std::list<unsigned int> *branchPoints;
    void *scanner;
#endif
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (branchPoints);
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
int yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, std::list<unsigned int> *branchPoints, void *scanner);
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
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, std::list<unsigned int> *branchPoints, void *scanner)
#else
int
yyparse (input, molList, branchPoints, scanner)
    const char *input;
    std::vector<RDKit::RWMol *> *molList;
    std::list<unsigned int> *branchPoints;
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
#line 79 "smiles.yy"
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

  case 4:
#line 88 "smiles.yy"
    {
  yyclearin;
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
;}
    break;

  case 5:
#line 94 "smiles.yy"
    {
  YYACCEPT;
;}
    break;

  case 6:
#line 97 "smiles.yy"
    {
  yyclearin;
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
;}
    break;

  case 7:
#line 107 "smiles.yy"
    {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  RDKit::RWMol *curMol = (*molList)[ sz ];
  (yyvsp[(1) - (1)].atom)->setProp("_SmilesStart",1);
  curMol->addAtom((yyvsp[(1) - (1)].atom));
  delete (yyvsp[(1) - (1)].atom);
  (yyval.moli) = sz;
;}
    break;

  case 8:
#line 118 "smiles.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[(2) - (2)].atom));
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  delete (yyvsp[(2) - (2)].atom);
;}
    break;

  case 9:
#line 128 "smiles.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[(3) - (3)].atom));
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
  mp->addBond((yyvsp[(2) - (3)].bond),true);
  delete (yyvsp[(3) - (3)].atom);
;}
    break;

  case 10:
#line 148 "smiles.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[(3) - (3)].atom));
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  delete (yyvsp[(3) - (3)].atom);
;}
    break;

  case 11:
#line 156 "smiles.yy"
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  mp->setAtomBookmark(atom,(yyvsp[(2) - (2)].ival));

  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::UNSPECIFIED);
  mp->setBondBookmark(newB,(yyvsp[(2) - (2)].ival));
  newB->setProp("_unspecifiedOrder",1);
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-((yyvsp[(2) - (2)].ival)+1));
  atom->setProp("_RingClosures",tmp);
;}
    break;

  case 12:
#line 173 "smiles.yy"
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     (yyvsp[(2) - (3)].bond)->getBondType());
  newB->setBondDir((yyvsp[(2) - (3)].bond)->getBondDir());
  mp->setAtomBookmark(atom,(yyvsp[(3) - (3)].ival));
  mp->setBondBookmark(newB,(yyvsp[(3) - (3)].ival));
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-((yyvsp[(3) - (3)].ival)+1));
  atom->setProp("_RingClosures",tmp);
  delete (yyvsp[(2) - (3)].bond);
;}
    break;

  case 13:
#line 190 "smiles.yy"
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::SINGLE);
  mp->setAtomBookmark(atom,(yyvsp[(3) - (3)].ival));
  mp->setBondBookmark(newB,(yyvsp[(3) - (3)].ival));
  INT_VECT tmp;
  if(atom->hasProp("_RingClosures")){
    atom->getProp("_RingClosures",tmp);
  }
  tmp.push_back(-((yyvsp[(3) - (3)].ival)+1));
  atom->setProp("_RingClosures",tmp);
;}
    break;

  case 14:
#line 205 "smiles.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[(3) - (3)].atom));
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  delete (yyvsp[(3) - (3)].atom);
  branchPoints->push_back(atomIdx1);
;}
    break;

  case 15:
#line 215 "smiles.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[(4) - (4)].atom));
  if( (yyvsp[(3) - (4)].bond)->getBondType() == Bond::DATIVER ){
    (yyvsp[(3) - (4)].bond)->setBeginAtomIdx(atomIdx1);
    (yyvsp[(3) - (4)].bond)->setEndAtomIdx(atomIdx2);
    (yyvsp[(3) - (4)].bond)->setBondType(Bond::DATIVE);
  }else if ( (yyvsp[(3) - (4)].bond)->getBondType() == Bond::DATIVEL ){
    (yyvsp[(3) - (4)].bond)->setBeginAtomIdx(atomIdx2);
    (yyvsp[(3) - (4)].bond)->setEndAtomIdx(atomIdx1);
    (yyvsp[(3) - (4)].bond)->setBondType(Bond::DATIVE);
  } else {
    (yyvsp[(3) - (4)].bond)->setBeginAtomIdx(atomIdx1);
    (yyvsp[(3) - (4)].bond)->setEndAtomIdx(atomIdx2);
  }
  mp->addBond((yyvsp[(3) - (4)].bond),true);
  delete (yyvsp[(4) - (4)].atom);
  branchPoints->push_back(atomIdx1);
;}
    break;

  case 16:
#line 235 "smiles.yy"
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[(4) - (4)].atom));
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  delete (yyvsp[(4) - (4)].atom);
  branchPoints->push_back(atomIdx1);
;}
    break;

  case 17:
#line 243 "smiles.yy"
    {
  if(branchPoints->empty()) yyerror(input,molList,branchPoints,scanner,"extra close parentheses");
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->setActiveAtom(branchPoints->back());
  branchPoints->pop_back();
;}
    break;

  case 19:
#line 255 "smiles.yy"
    {
  (yyval.atom) = (yyvsp[(2) - (5)].atom);
  (yyval.atom)->setNoImplicit(true);
  (yyval.atom)->setProp("molAtomMapNumber",(yyvsp[(4) - (5)].ival));
;}
    break;

  case 20:
#line 262 "smiles.yy"
    {
  (yyval.atom) = (yyvsp[(2) - (3)].atom);
  (yyvsp[(2) - (3)].atom)->setNoImplicit(true);
;}
    break;

  case 22:
#line 270 "smiles.yy"
    { (yyvsp[(1) - (2)].atom)->setFormalCharge(1); ;}
    break;

  case 23:
#line 271 "smiles.yy"
    { (yyvsp[(1) - (3)].atom)->setFormalCharge(2); ;}
    break;

  case 24:
#line 272 "smiles.yy"
    { (yyvsp[(1) - (3)].atom)->setFormalCharge((yyvsp[(3) - (3)].ival)); ;}
    break;

  case 25:
#line 273 "smiles.yy"
    { (yyvsp[(1) - (2)].atom)->setFormalCharge(-1); ;}
    break;

  case 26:
#line 274 "smiles.yy"
    { (yyvsp[(1) - (3)].atom)->setFormalCharge(-2); ;}
    break;

  case 27:
#line 275 "smiles.yy"
    { (yyvsp[(1) - (3)].atom)->setFormalCharge(-(yyvsp[(3) - (3)].ival)); ;}
    break;

  case 28:
#line 279 "smiles.yy"
    { (yyval.atom) = new Atom(1); ;}
    break;

  case 29:
#line 280 "smiles.yy"
    { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[(1) - (2)].ival)); ;}
    break;

  case 30:
#line 281 "smiles.yy"
    { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs(1); ;}
    break;

  case 31:
#line 282 "smiles.yy"
    { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[(1) - (3)].ival)); (yyval.atom)->setNumExplicitHs(1);;}
    break;

  case 32:
#line 283 "smiles.yy"
    { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs((yyvsp[(3) - (3)].ival)); ;}
    break;

  case 33:
#line 284 "smiles.yy"
    { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[(1) - (4)].ival)); (yyval.atom)->setNumExplicitHs((yyvsp[(4) - (4)].ival));;}
    break;

  case 35:
#line 286 "smiles.yy"
    { (yyval.atom) = (yyvsp[(1) - (2)].atom); (yyvsp[(1) - (2)].atom)->setNumExplicitHs(1);;}
    break;

  case 36:
#line 287 "smiles.yy"
    { (yyval.atom) = (yyvsp[(1) - (3)].atom); (yyvsp[(1) - (3)].atom)->setNumExplicitHs((yyvsp[(3) - (3)].ival));;}
    break;

  case 38:
#line 292 "smiles.yy"
    { (yyvsp[(1) - (2)].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); ;}
    break;

  case 39:
#line 293 "smiles.yy"
    { (yyvsp[(1) - (3)].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); ;}
    break;

  case 41:
#line 298 "smiles.yy"
    { (yyvsp[(2) - (2)].atom)->setIsotope( (yyvsp[(1) - (2)].ival) ); (yyval.atom) = (yyvsp[(2) - (2)].atom); ;}
    break;

  case 43:
#line 300 "smiles.yy"
    { (yyvsp[(2) - (2)].atom)->setIsotope( (yyvsp[(1) - (2)].ival) ); (yyval.atom) = (yyvsp[(2) - (2)].atom); ;}
    break;

  case 47:
#line 310 "smiles.yy"
    { (yyval.ival) = (yyvsp[(2) - (3)].ival)*10+(yyvsp[(3) - (3)].ival); ;}
    break;

  case 51:
#line 320 "smiles.yy"
    { (yyval.ival) = (yyvsp[(1) - (2)].ival)*10 + (yyvsp[(2) - (2)].ival); ;}
    break;


/* Line 1267 of yacc.c.  */
#line 1846 "/Users/landrgr1/RDKit_trunk/Code/GraphMol/SmilesParse/smiles.tab.cpp"
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
      yyerror (input, molList, branchPoints, scanner, YY_("syntax error"));
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
	    yyerror (input, molList, branchPoints, scanner, yymsg);
	  }
	else
	  {
	    yyerror (input, molList, branchPoints, scanner, YY_("syntax error"));
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
		      yytoken, &yylval, input, molList, branchPoints, scanner);
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
		  yystos[yystate], yyvsp, input, molList, branchPoints, scanner);
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
  yyerror (input, molList, branchPoints, scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval, input, molList, branchPoints, scanner);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp, input, molList, branchPoints, scanner);
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


#line 334 "smiles.yy"




