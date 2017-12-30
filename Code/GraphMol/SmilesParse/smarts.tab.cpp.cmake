/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

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
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         yysmarts_parse
#define yylex           yysmarts_lex
#define yyerror         yysmarts_error
#define yydebug         yysmarts_debug
#define yynerrs         yysmarts_nerrs


/* Copy the first part of user declarations.  */
#line 1 "smarts.yy" /* yacc.c:339  */


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

#define YYDEBUG 1
#include "smarts.tab.hpp"

extern int yysmarts_lex(YYSTYPE *,void *);

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

#line 116 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "smarts.tab.hpp".  */
#ifndef YY_YYSMARTS_MNT_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
# define YY_YYSMARTS_MNT_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
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
    RANGE_OPEN_TOKEN = 272,
    RANGE_CLOSE_TOKEN = 273,
    HASH_TOKEN = 274,
    MINUS_TOKEN = 275,
    PLUS_TOKEN = 276,
    CHIRAL_MARKER_TOKEN = 277,
    CHI_CLASS_TOKEN = 278,
    CHI_CLASS_OH_TOKEN = 279,
    H_TOKEN = 280,
    AT_TOKEN = 281,
    PERCENT_TOKEN = 282,
    ATOM_OPEN_TOKEN = 283,
    ATOM_CLOSE_TOKEN = 284,
    NOT_TOKEN = 285,
    AND_TOKEN = 286,
    OR_TOKEN = 287,
    SEMI_TOKEN = 288,
    BEGIN_RECURSE = 289,
    END_RECURSE = 290,
    COLON_TOKEN = 291,
    UNDERSCORE_TOKEN = 292,
    BOND_TOKEN = 293,
    EOS_TOKEN = 294
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 51 "smarts.yy" /* yacc.c:355  */

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;

#line 203 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner);

#endif /* !YY_YYSMARTS_MNT_C_USERS_GLANDRUM_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 219 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:358  */

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
#else
typedef signed char yytype_int8;
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
# elif ! defined YYSIZE_T
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
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
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
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
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

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  35
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   430

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  40
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  94
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  146

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   294

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
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
      35,    36,    37,    38,    39
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    87,    87,    88,    94,    97,   108,   117,   146,   166,
     172,   207,   227,   242,   243,   254,   255,   256,   260,   283,
     287,   292,   298,   307,   312,   319,   326,   343,   348,   353,
     358,   363,   366,   370,   371,   375,   392,   416,   417,   421,
     422,   426,   427,   431,   432,   435,   442,   443,   446,   447,
     450,   451,   454,   457,   460,   465,   470,   475,   481,   487,
     488,   496,   508,   513,   518,   522,   526,   530,   533,   534,
     541,   542,   548,   554,   560,   565,   572,   573,   574,   575,
     576,   577,   581,   582,   583,   584,   585,   586,   587,   592,
     593,   597,   598,   601,   602
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "AROMATIC_ATOM_TOKEN",
  "ORGANIC_ATOM_TOKEN", "ATOM_TOKEN", "SIMPLE_ATOM_QUERY_TOKEN",
  "COMPLEX_ATOM_QUERY_TOKEN", "RINGSIZE_ATOM_QUERY_TOKEN",
  "RINGBOND_ATOM_QUERY_TOKEN", "IMPLICIT_H_ATOM_QUERY_TOKEN", "HYB_TOKEN",
  "ZERO_TOKEN", "NONZERO_DIGIT_TOKEN", "GROUP_OPEN_TOKEN",
  "GROUP_CLOSE_TOKEN", "SEPARATOR_TOKEN", "RANGE_OPEN_TOKEN",
  "RANGE_CLOSE_TOKEN", "HASH_TOKEN", "MINUS_TOKEN", "PLUS_TOKEN",
  "CHIRAL_MARKER_TOKEN", "CHI_CLASS_TOKEN", "CHI_CLASS_OH_TOKEN",
  "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN", "ATOM_OPEN_TOKEN",
  "ATOM_CLOSE_TOKEN", "NOT_TOKEN", "AND_TOKEN", "OR_TOKEN", "SEMI_TOKEN",
  "BEGIN_RECURSE", "END_RECURSE", "COLON_TOKEN", "UNDERSCORE_TOKEN",
  "BOND_TOKEN", "EOS_TOKEN", "$accept", "cmpd", "mol", "branch", "atomd",
  "hydrogen_atom", "atom_expr", "point_query", "recursive_query",
  "atom_query", "simple_atom", "bond_expr", "bond_query", "bondd",
  "charge_spec", "ring_number", "number", "nonzero_number", "digit", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294
};
# endif

#define YYPACT_NINF -37

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-37)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
      25,   -36,   -37,   -37,   -37,   339,    13,   205,   -37,   -37,
     -37,   -37,   -37,   105,    88,    88,    88,   -37,   -37,   -37,
      88,    37,    68,   387,     1,   371,    30,   241,   -37,   -37,
     -37,   -18,   -37,   125,   102,   -37,    -9,   -37,   -37,   -37,
      18,    30,   -37,   -37,   -37,    29,   392,   -37,   -37,   -37,
     -37,   382,   392,   -37,   -37,   -37,    88,   -37,   -37,   -37,
     -37,   -37,   -37,   -37,   -37,   -37,   -37,    88,   -24,   -37,
     -37,    88,   -37,   167,   121,   -37,   371,   371,   371,    88,
     -37,    88,   -37,    88,    43,   -37,   -37,   -37,   149,    80,
     -37,   102,   102,   -37,   392,   392,   392,   -37,   -37,   -37,
      12,    38,   -37,    88,    34,   371,   307,   275,    44,   -37,
     -37,   -37,    88,    41,   -37,   177,   -37,   130,   -37,    57,
     107,    88,   -37,    46,    84,   -37,    64,   -37,    88,   -37,
     -37,   186,   108,   -37,   102,   -37,   103,   -37,   214,   -37,
     -37,   -37,   224,   -37,   131,   -37
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    62,    61,    63,     0,     0,     2,     6,    16,
      15,     5,    39,    43,    46,    48,    50,    59,    89,    91,
       0,    81,    78,    55,    58,     0,     0,     0,    31,    33,
      34,    37,    56,    60,    90,     1,     0,     4,    94,    93,
       0,     0,    72,    71,    74,     0,     0,    73,    70,    12,
       7,     0,    67,    68,    10,    82,     0,    44,    47,    49,
      51,    41,    79,    80,    76,    77,    19,     0,     0,    54,
      57,    55,    32,    60,     0,    17,     0,     0,     0,     0,
      30,    53,    40,     0,     0,    38,    92,     3,     0,     0,
       9,     0,     0,    75,     0,     0,     0,     8,    11,    69,
       0,     0,    23,     0,    35,    27,    28,    29,     0,    52,
      42,    21,     0,     0,    13,     0,    83,     0,    64,    65,
      66,     0,    20,     0,     0,    18,     0,    25,     0,    14,
      84,     0,     0,    24,    36,    22,     0,    85,     0,    45,
      26,    86,     0,    87,     0,    88
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -37,   -37,   -20,   -37,    -6,   -37,   179,    -2,   -37,   -37,
      14,     0,   -37,     9,   -19,    85,    -5,    36,   -32
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     6,     7,    49,     8,     9,    27,    28,    29,    30,
      10,    51,    52,    53,    32,    54,    73,    34,    55
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      33,    50,    86,    11,    68,   102,    74,    81,    57,    58,
      59,    60,   103,    35,    36,    61,    63,    65,    69,    31,
      88,     2,     3,    72,     4,    80,     1,    70,     2,     3,
      87,     4,   121,     2,     3,    90,     4,    42,    43,    31,
      89,    31,    91,    92,    44,    97,     5,    85,    46,    18,
      19,   100,    37,     5,    47,    93,    48,    62,     5,   116,
     117,    99,   101,    21,    22,   113,    69,   122,    50,   115,
     127,   124,   111,   125,   108,   133,   109,   128,   110,   112,
      18,    19,    50,     2,     3,   131,     4,    85,    94,    64,
      31,    31,    31,   135,   118,   119,   120,    19,   123,   138,
      18,    19,    86,    80,    80,    80,   142,   126,     5,    50,
     144,    94,    95,    96,    38,    39,   132,    18,    19,    31,
      31,    31,    56,   136,     2,     3,   139,     4,     2,     3,
      82,     4,   140,    38,    39,    40,    98,    41,    94,    95,
      42,    43,    38,    39,    83,   130,   145,    44,    45,     5,
      84,    46,     2,     3,     0,     4,   104,    47,     0,    48,
     134,    38,    39,    40,   114,    41,     0,     0,    42,    43,
       2,     3,    82,     4,     0,    44,    45,     5,     0,    46,
       2,     3,     0,     4,     0,    47,    83,    48,     0,    38,
      39,    40,   129,    41,     0,     0,    42,    43,    38,    39,
       0,   137,     0,    44,    45,     5,     0,    46,     2,     3,
       0,     4,     0,    47,     0,    48,     0,    38,    39,    40,
       0,    41,     0,     0,    42,    43,    38,    39,     0,   141,
       0,    44,    45,     5,     0,    46,    38,    39,     0,   143,
       0,    47,     0,    48,     2,     3,    12,     4,    13,    14,
      15,    16,    17,    18,    19,   105,   106,   107,     0,     0,
      20,    21,    22,     0,     0,     0,    71,    24,     0,     0,
      75,    25,    76,    77,    78,    26,     0,    79,     2,     3,
      12,     4,    13,    14,    15,    16,    17,    18,    19,     0,
       0,     0,     0,     0,    20,    21,    22,     0,     0,     0,
      71,    24,     0,     0,     0,    25,    76,    77,     0,    26,
       2,     3,    12,     4,    13,    14,    15,    16,    17,    18,
      19,     0,     0,     0,     0,     0,    20,    21,    22,     0,
       0,     0,    71,    24,     0,     0,     0,    25,    76,     0,
       0,    26,     2,     3,    12,     4,    13,    14,    15,    16,
      17,    18,    19,     0,     0,     0,     0,     0,    20,    21,
      22,     0,     0,     0,    23,    24,     0,     0,     0,    25,
       0,     0,     0,    26,     2,     3,    12,     4,    13,    14,
      15,    16,    17,    18,    19,     2,     3,     0,     4,     0,
      20,    21,    22,     0,    38,    39,    71,    24,     0,    18,
      19,    25,     0,     0,     0,    26,     0,    21,    22,    45,
       5,    42,    43,    94,    95,    96,    66,     0,    44,     0,
       0,     0,    46,    67,     0,     0,     0,     0,    47,     0,
      48
};

static const yytype_int16 yycheck[] =
{
       5,     7,    34,    39,    23,    29,    26,    25,    13,    14,
      15,    16,    36,     0,     1,    20,    21,    22,    23,     5,
      40,     3,     4,    25,     6,    27,     1,    26,     3,     4,
      39,     6,    20,     3,     4,    41,     6,    19,    20,    25,
      40,    27,    13,    14,    26,    51,    28,    33,    30,    12,
      13,    56,    39,    28,    36,    46,    38,    20,    28,    91,
      92,    52,    67,    20,    21,    84,    71,    29,    74,    89,
      29,    37,    29,    29,    79,    29,    81,    36,    83,    36,
      12,    13,    88,     3,     4,   117,     6,    73,    31,    21,
      76,    77,    78,    29,    94,    95,    96,    13,   103,   131,
      12,    13,   134,   105,   106,   107,   138,   112,    28,   115,
     142,    31,    32,    33,    12,    13,   121,    12,    13,   105,
     106,   107,    17,   128,     3,     4,    18,     6,     3,     4,
       5,     6,    29,    12,    13,    14,    51,    16,    31,    32,
      19,    20,    12,    13,    19,    15,    15,    26,    27,    28,
      25,    30,     3,     4,    -1,     6,    35,    36,    -1,    38,
     124,    12,    13,    14,    15,    16,    -1,    -1,    19,    20,
       3,     4,     5,     6,    -1,    26,    27,    28,    -1,    30,
       3,     4,    -1,     6,    -1,    36,    19,    38,    -1,    12,
      13,    14,    15,    16,    -1,    -1,    19,    20,    12,    13,
      -1,    15,    -1,    26,    27,    28,    -1,    30,     3,     4,
      -1,     6,    -1,    36,    -1,    38,    -1,    12,    13,    14,
      -1,    16,    -1,    -1,    19,    20,    12,    13,    -1,    15,
      -1,    26,    27,    28,    -1,    30,    12,    13,    -1,    15,
      -1,    36,    -1,    38,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    76,    77,    78,    -1,    -1,
      19,    20,    21,    -1,    -1,    -1,    25,    26,    -1,    -1,
      29,    30,    31,    32,    33,    34,    -1,    36,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    -1,
      -1,    -1,    -1,    -1,    19,    20,    21,    -1,    -1,    -1,
      25,    26,    -1,    -1,    -1,    30,    31,    32,    -1,    34,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    -1,    -1,    -1,    -1,    -1,    19,    20,    21,    -1,
      -1,    -1,    25,    26,    -1,    -1,    -1,    30,    31,    -1,
      -1,    34,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    -1,    -1,    -1,    -1,    -1,    19,    20,
      21,    -1,    -1,    -1,    25,    26,    -1,    -1,    -1,    30,
      -1,    -1,    -1,    34,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,     3,     4,    -1,     6,    -1,
      19,    20,    21,    -1,    12,    13,    25,    26,    -1,    12,
      13,    30,    -1,    -1,    -1,    34,    -1,    20,    21,    27,
      28,    19,    20,    31,    32,    33,    29,    -1,    26,    -1,
      -1,    -1,    30,    36,    -1,    -1,    -1,    -1,    36,    -1,
      38
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     6,    28,    41,    42,    44,    45,
      50,    39,     5,     7,     8,     9,    10,    11,    12,    13,
      19,    20,    21,    25,    26,    30,    34,    46,    47,    48,
      49,    50,    54,    56,    57,     0,     1,    39,    12,    13,
      14,    16,    19,    20,    26,    27,    30,    36,    38,    43,
      44,    51,    52,    53,    55,    58,    17,    56,    56,    56,
      56,    56,    20,    56,    21,    56,    29,    36,    54,    56,
      26,    25,    47,    56,    42,    29,    31,    32,    33,    36,
      47,    25,     5,    19,    25,    50,    58,    39,    42,    51,
      44,    13,    14,    53,    31,    32,    33,    44,    55,    53,
      56,    56,    29,    36,    35,    46,    46,    46,    56,    56,
      56,    29,    36,    54,    15,    42,    58,    58,    51,    51,
      51,    20,    29,    56,    37,    29,    56,    29,    36,    15,
      15,    58,    56,    29,    57,    29,    56,    15,    58,    18,
      29,    15,    58,    15,    58,    15
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    40,    41,    41,    41,    41,    42,    42,    42,    42,
      42,    42,    42,    43,    43,    44,    44,    44,    44,    45,
      45,    45,    45,    45,    45,    45,    45,    46,    46,    46,
      46,    46,    47,    47,    47,    48,    48,    49,    49,    49,
      49,    49,    49,    49,    49,    49,    49,    49,    49,    49,
      49,    49,    49,    49,    49,    49,    49,    49,    49,    49,
      49,    50,    50,    50,    51,    51,    51,    51,    52,    52,
      53,    53,    53,    53,    53,    53,    54,    54,    54,    54,
      54,    54,    55,    55,    55,    55,    55,    55,    55,    56,
      56,    57,    57,    58,    58
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     3,     2,     2,     1,     2,     3,     3,
       2,     3,     2,     3,     4,     1,     1,     3,     5,     3,
       5,     4,     6,     4,     6,     5,     7,     3,     3,     3,
       2,     1,     2,     1,     1,     3,     5,     1,     2,     1,
       2,     2,     3,     1,     2,     6,     1,     2,     1,     2,
       1,     2,     3,     2,     2,     1,     1,     2,     1,     1,
       1,     1,     1,     1,     3,     3,     3,     1,     1,     2,
       1,     1,     1,     1,     1,     2,     2,     2,     1,     2,
       2,     1,     1,     3,     4,     5,     6,     7,     8,     1,
       1,     1,     2,     1,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (input, molList, scanner, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, input, molList, scanner); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (scanner);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, scanner);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              , input, molList, scanner);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, molList, scanner); \
} while (0)

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
#ifndef YYINITDEPTH
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
static YYSIZE_T
yystrlen (const char *yystr)
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
static char *
yystpcpy (char *yydest, const char *yysrc)
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

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (scanner);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
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
  int yytoken = 0;
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

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
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
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (&yylval, scanner);
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
      if (yytable_value_is_error (yyn))
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
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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
     '$$ = $1'.

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
#line 88 "smarts.yy" /* yacc.c:1646  */
    {
  yyclearin;
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1476 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 4:
#line 94 "smarts.yy" /* yacc.c:1646  */
    {
  YYACCEPT;
}
#line 1484 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 5:
#line 97 "smarts.yy" /* yacc.c:1646  */
    {
  yyclearin;
  yyerrok;

  yyErrorCleanup(molList);
  YYABORT;
}
#line 1496 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 108 "smarts.yy" /* yacc.c:1646  */
    {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1);
  (*molList)[ sz ]->addAtom((yyvsp[0].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
}
#line 1510 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 117 "smarts.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);

  QueryBond *newB;
  // this is a bit of a hack to try and get nicer "SMILES" from
  // a SMARTS molecule:
  if(!(a1->getIsAromatic() && (yyvsp[0].atom)->getIsAromatic())){
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
  newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  mp->addBond(newB);
  delete newB;
  //delete $2;
}
#line 1543 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 8:
#line 146 "smarts.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[0].atom),true,true);
  if( (yyvsp[-1].bond)->getBondType() == Bond::DATIVER ){
    (yyvsp[-1].bond)->setBeginAtomIdx(atomIdx1);
    (yyvsp[-1].bond)->setEndAtomIdx(atomIdx2);
    (yyvsp[-1].bond)->setBondType(Bond::DATIVE);
  }else if ( (yyvsp[-1].bond)->getBondType() == Bond::DATIVEL ){
    (yyvsp[-1].bond)->setBeginAtomIdx(atomIdx2);
    (yyvsp[-1].bond)->setEndAtomIdx(atomIdx1);
    (yyvsp[-1].bond)->setBondType(Bond::DATIVE);
  } else {
    (yyvsp[-1].bond)->setBeginAtomIdx(atomIdx1);
    (yyvsp[-1].bond)->setEndAtomIdx(atomIdx2);
  }
  mp->addBond((yyvsp[-1].bond));
  delete (yyvsp[-1].bond);
}
#line 1567 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 9:
#line 166 "smarts.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1577 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 10:
#line 172 "smarts.yy" /* yacc.c:1646  */
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
  newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atom->getIdx());
  mp->setBondBookmark(newB,(yyvsp[0].ival));
  mp->setAtomBookmark(atom,(yyvsp[0].ival));

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  if(atom->hasProp(RDKit::common_properties::_RingClosures)){
    atom->getProp(RDKit::common_properties::_RingClosures,tmp);
  }
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);

}
#line 1616 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 11:
#line 207 "smarts.yy" /* yacc.c:1646  */
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();

  mp->setBondBookmark((yyvsp[-1].bond),(yyvsp[0].ival));
  (yyvsp[-1].bond)->setOwningMol(mp);
  (yyvsp[-1].bond)->setBeginAtomIdx(atom->getIdx());
  mp->setAtomBookmark(atom,(yyvsp[0].ival));

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  if(atom->hasProp(RDKit::common_properties::_RingClosures)){
    atom->getProp(RDKit::common_properties::_RingClosures,tmp);
  }
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);

}
#line 1640 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 12:
#line 227 "smarts.yy" /* yacc.c:1646  */
    {
  RWMol *m1_p = (*molList)[(yyval.moli)],*m2_p=(*molList)[(yyvsp[0].moli)];
  m2_p->getAtomWithIdx(0)->setProp(RDKit::common_properties::_SmilesStart,1);
  // FIX: handle generic bonds here
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED,Bond::NONE,false,true);
  delete m2_p;
  int sz = molList->size();
  if ( sz==(yyvsp[0].moli)+1) {
    molList->resize( sz-1 );
  }
}
#line 1656 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 13:
#line 242 "smarts.yy" /* yacc.c:1646  */
    { (yyval.moli) = (yyvsp[-1].moli); }
#line 1662 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 14:
#line 243 "smarts.yy" /* yacc.c:1646  */
    {
  // FIX: this needs to handle arbitrary bond_exprs
  (yyval.moli) = (yyvsp[-1].moli);
  int sz     = molList->size();
  (yyvsp[-2].bond)->setOwningMol((*molList)[ sz-1 ]);
  (yyvsp[-2].bond)->setBeginAtomIdx(0);
  (*molList)[ sz-1 ]->setBondBookmark((yyvsp[-2].bond),ci_LEADING_BOND);
}
#line 1675 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 17:
#line 257 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = (yyvsp[-1].atom);
}
#line 1683 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 18:
#line 261 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1692 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 19:
#line 284 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom(1);
}
#line 1700 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 20:
#line 288 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1709 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 21:
#line 292 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1720 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 22:
#line 298 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-4].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-4].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1733 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 23:
#line 307 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1743 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 24:
#line 312 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1755 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 25:
#line 319 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1767 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 26:
#line 326 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-5].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-5].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1781 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 27:
#line 343 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  delete (yyvsp[0].atom);
}
#line 1791 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 28:
#line 348 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  delete (yyvsp[0].atom);
}
#line 1801 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 29:
#line 353 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  delete (yyvsp[0].atom);
}
#line 1811 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 30:
#line 358 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-1].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-1].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  delete (yyvsp[0].atom);
}
#line 1821 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 32:
#line 366 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->getQuery()->setNegation(!((yyvsp[0].atom)->getQuery()->getNegation()));
  (yyval.atom) = (yyvsp[0].atom);
}
#line 1830 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 35:
#line 375 "smarts.yy" /* yacc.c:1646  */
    {
  // this is a recursive SMARTS expression
  QueryAtom *qA = new QueryAtom();
  //  FIX: there's maybe a leak here
  RWMol *molP = (*molList)[(yyvsp[-1].moli)];
  // close any rings in the molecule:
  SmilesParseOps::CloseMolRings(molP,0);

  //molP->debugMol(std::cout);
  qA->setQuery(new RecursiveStructureQuery(molP));
  //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
  int sz = molList->size();
  if ( sz==(yyvsp[-1].moli)+1) {
    molList->resize( sz-1 );
  }
  (yyval.atom) = qA;
}
#line 1852 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 36:
#line 392 "smarts.yy" /* yacc.c:1646  */
    {
  // UNDOCUMENTED EXTENSION:
  // this is a recursive SMARTS expression with a serial number
  // please don't write your own SMARTS that include this extension:
  // the RDKit smarts parsing code will automatically insert serial
  // numbers for recursive smarts patterns.
  QueryAtom *qA = new QueryAtom();
  //  FIX: there's maybe a leak here
  RWMol *molP = (*molList)[(yyvsp[-3].moli)];
  // close any rings in the molecule:
  SmilesParseOps::CloseMolRings(molP,0);

  //molP->debugMol(std::cout);
  qA->setQuery(new RecursiveStructureQuery(molP,(yyvsp[0].ival)));
  //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
  int sz = molList->size();
  if ( sz==(yyvsp[-3].moli)+1) {
    molList->resize( sz-1 );
  }
  (yyval.atom) = qA;
}
#line 1878 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 38:
#line 417 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 1887 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 422 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 1896 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 41:
#line 426 "smarts.yy" /* yacc.c:1646  */
    { (yyval.atom) = new QueryAtom((yyvsp[0].ival)); }
#line 1902 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 42:
#line 427 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
}
#line 1911 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 44:
#line 432 "smarts.yy" /* yacc.c:1646  */
    {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
}
#line 1919 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 45:
#line 435 "smarts.yy" /* yacc.c:1646  */
    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
}
#line 1931 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 47:
#line 443 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
}
#line 1939 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 49:
#line 447 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
}
#line 1947 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 51:
#line 451 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
}
#line 1955 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 52:
#line 454 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom)->expandQuery(makeAtomHCountQuery((yyvsp[0].ival)),Queries::COMPOSITE_AND);
}
#line 1963 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 53:
#line 457 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom)->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND);
}
#line 1971 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 54:
#line 460 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 1981 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 55:
#line 465 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  (yyval.atom)=newQ;
}
#line 1991 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 56:
#line 470 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 2001 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 57:
#line 475 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
#line 2012 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 58:
#line 481 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
#line 2023 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 60:
#line 488 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 2033 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 61:
#line 496 "smarts.yy" /* yacc.c:1646  */
    {
  //
  // This construction (and some others) may seem odd, but the
  // SMARTS definition requires that an atom which is aliphatic on
  // input (i.e. something in the "organic subset" that is given with
  // a capital letter) only match aliphatic atoms.
  //
  // The following rule applies a similar logic to aromatic atoms.
  //
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->expandQuery(makeAtomAliphaticQuery(),Queries::COMPOSITE_AND);
}
#line 2050 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 62:
#line 508 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->expandQuery(makeAtomAromaticQuery(),Queries::COMPOSITE_AND);
}
#line 2060 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 64:
#line 518 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2069 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 65:
#line 522 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
}
#line 2078 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 66:
#line 526 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2087 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 69:
#line 534 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2096 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 71:
#line 542 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
#line 2107 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 72:
#line 548 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
#line 2118 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 73:
#line 554 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
#line 2129 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 74:
#line 560 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
#line 2139 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 75:
#line 565 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
#line 2148 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 76:
#line 572 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=2; }
#line 2154 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 77:
#line 573 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=(yyvsp[0].ival); }
#line 2160 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 78:
#line 574 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=1; }
#line 2166 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 79:
#line 575 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=-2; }
#line 2172 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 80:
#line 576 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=-(yyvsp[0].ival); }
#line 2178 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 81:
#line 577 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=-1; }
#line 2184 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 83:
#line 582 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 2190 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 84:
#line 583 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 2196 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 85:
#line 584 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2202 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 86:
#line 585 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2208 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 87:
#line 586 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2214 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 88:
#line 587 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2220 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 92:
#line 598 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
#line 2226 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;


#line 2230 "/mnt/c/Users/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (input, molList, scanner, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (input, molList, scanner, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
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

  /* Do not reclaim the symbols of the rule whose action triggered
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
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
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

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


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

#if !defined yyoverflow || YYERROR_VERBOSE
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
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, input, molList, scanner);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
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
  return yyresult;
}
#line 605 "smarts.yy" /* yacc.c:1906  */

