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


  //
  //  Copyright (C) 2003-2018 Greg Landrum and Rational Discovery LLC
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

#line 115 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:339  */

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
#ifndef YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
# define YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
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
    HETERONEIGHBOR_ATOM_QUERY_TOKEN = 267,
    ALIPHATIC = 268,
    ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN = 269,
    ZERO_TOKEN = 270,
    NONZERO_DIGIT_TOKEN = 271,
    GROUP_OPEN_TOKEN = 272,
    GROUP_CLOSE_TOKEN = 273,
    SEPARATOR_TOKEN = 274,
    RANGE_OPEN_TOKEN = 275,
    RANGE_CLOSE_TOKEN = 276,
    HASH_TOKEN = 277,
    MINUS_TOKEN = 278,
    PLUS_TOKEN = 279,
    CHIRAL_MARKER_TOKEN = 280,
    CHI_CLASS_TOKEN = 281,
    CHI_CLASS_OH_TOKEN = 282,
    H_TOKEN = 283,
    AT_TOKEN = 284,
    PERCENT_TOKEN = 285,
    ATOM_OPEN_TOKEN = 286,
    ATOM_CLOSE_TOKEN = 287,
    NOT_TOKEN = 288,
    AND_TOKEN = 289,
    OR_TOKEN = 290,
    SEMI_TOKEN = 291,
    BEGIN_RECURSE = 292,
    END_RECURSE = 293,
    COLON_TOKEN = 294,
    UNDERSCORE_TOKEN = 295,
    BOND_TOKEN = 296,
    EOS_TOKEN = 297
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 50 "smarts.yy" /* yacc.c:355  */

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;

#line 205 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, void *scanner);

#endif /* !YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 221 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:358  */

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
#define YYFINAL  38
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   470

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  43
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  20
/* YYNRULES -- Number of rules.  */
#define YYNRULES  106
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  155

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   297

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
      35,    36,    37,    38,    39,    40,    41,    42
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    87,    87,    88,    94,    97,   108,   116,   145,   165,
     170,   205,   225,   239,   240,   251,   252,   253,   257,   280,
     284,   289,   295,   304,   310,   318,   326,   339,   345,   351,
     357,   380,   383,   389,   390,   394,   411,   435,   436,   441,
     442,   447,   448,   453,   454,   455,   456,   457,   458,   459,
     462,   465,   468,   471,   474,   477,   483,   489,   496,   504,
     512,   518,   524,   530,   536,   542,   543,   550,   551,   554,
     557,   560,   563,   569,   581,   586,   591,   595,   599,   603,
     606,   607,   614,   615,   621,   627,   633,   638,   645,   646,
     647,   648,   649,   650,   654,   655,   656,   657,   658,   659,
     660,   665,   666,   670,   671,   674,   675
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
  "HETERONEIGHBOR_ATOM_QUERY_TOKEN", "ALIPHATIC",
  "ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN", "ZERO_TOKEN",
  "NONZERO_DIGIT_TOKEN", "GROUP_OPEN_TOKEN", "GROUP_CLOSE_TOKEN",
  "SEPARATOR_TOKEN", "RANGE_OPEN_TOKEN", "RANGE_CLOSE_TOKEN", "HASH_TOKEN",
  "MINUS_TOKEN", "PLUS_TOKEN", "CHIRAL_MARKER_TOKEN", "CHI_CLASS_TOKEN",
  "CHI_CLASS_OH_TOKEN", "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN",
  "ATOM_OPEN_TOKEN", "ATOM_CLOSE_TOKEN", "NOT_TOKEN", "AND_TOKEN",
  "OR_TOKEN", "SEMI_TOKEN", "BEGIN_RECURSE", "END_RECURSE", "COLON_TOKEN",
  "UNDERSCORE_TOKEN", "BOND_TOKEN", "EOS_TOKEN", "$accept", "cmpd", "mol",
  "branch", "atomd", "hydrogen_atom", "atom_expr", "point_query",
  "recursive_query", "atom_query", "possible_range_query", "simple_atom",
  "bond_expr", "bond_query", "bondd", "charge_spec", "ring_number",
  "number", "nonzero_number", "digit", YY_NULLPTR
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
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297
};
# endif

#define YYPACT_NINF -57

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-57)))

#define YYTABLE_NINF -73

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
      27,   -26,   -57,   -57,   -57,   369,    22,   202,   -57,   -57,
     -57,   -57,   -57,   113,   162,   222,   224,   -57,   233,   235,
     -57,   -57,    69,    80,    25,    20,     5,   404,     1,   262,
     -57,   -57,   -57,    30,   -57,   -57,    51,   124,   -57,    52,
     -57,   -57,   -57,   223,     1,   -57,   -57,   -57,   185,   429,
     -57,   -57,   -57,   -57,   419,   429,   -57,   -57,   -57,   -57,
     -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,
     -57,    69,    21,   -57,   -57,    69,   -57,    71,   119,   -57,
     404,   404,   404,    69,   -57,    82,   -57,    69,   424,   -57,
     -57,   -57,   150,   140,   -57,   124,   124,   -57,   429,   429,
     429,   -57,   -57,   -57,    15,   -57,    69,    69,     2,   404,
     334,   299,    31,    69,    77,   -57,   -57,    69,    49,   -57,
     -57,   181,   -57,   146,   -57,    73,   253,   -57,    84,    98,
     -57,    97,   116,    92,   -57,    69,   -57,   -57,   177,   -57,
     124,   -57,   -57,   105,   -57,   115,   -57,   242,   -57,   -57,
     -57,   264,   -57,   127,   -57
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    74,    73,    75,     0,     0,     2,     6,    16,
      15,     5,    39,    43,    46,    47,    48,    65,    44,    45,
     101,   103,     0,    93,    90,    61,    64,     0,     0,     0,
      31,    33,    34,     0,    37,    62,    66,   102,     1,     0,
       4,   106,   105,     0,     0,    84,    83,    86,     0,     0,
      85,    82,    12,     7,     0,    79,    80,    10,    94,    49,
      52,    53,    54,    50,    51,    41,    91,    92,    88,    89,
      19,     0,     0,    60,    63,    61,    32,    66,     0,    17,
       0,     0,     0,     0,    30,     0,    40,     0,    58,    38,
     104,     3,     0,     0,     9,     0,     0,    87,     0,     0,
       0,     8,    11,    81,     0,    23,     0,    58,    35,    27,
      28,    29,     0,     0,     0,    42,    21,     0,     0,    59,
      13,     0,    95,     0,    76,    77,    78,    20,     0,     0,
      18,     0,     0,     0,    25,     0,    14,    96,     0,    24,
      36,    55,    56,     0,    22,     0,    97,     0,    57,    26,
      98,     0,    99,     0,   100
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -57,   -57,   -22,   -57,    -6,   -57,   -56,     0,   -57,   -57,
     -57,    10,   -31,   -57,    -4,   -23,   101,    -5,    34,   -34
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     6,     7,    52,     8,     9,    29,    30,    31,    32,
      33,    10,    54,    55,    56,    35,    57,    77,    37,    58
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      36,    53,    72,    90,     2,     3,    78,     4,    59,    60,
      61,    62,    93,    63,    64,    34,    11,    65,    67,    69,
      73,    92,    38,    39,   109,   110,   111,    76,     1,    84,
       2,     3,     5,     4,    74,    20,    21,    34,    94,    34,
      20,    21,   129,    23,    24,    97,    89,   127,   101,    68,
      85,   103,    70,   105,     2,     3,    86,     4,     5,    71,
     106,   122,   123,   130,    40,   118,   104,   124,   125,   126,
      73,   121,    53,    87,     2,     3,    86,     4,   112,    88,
     114,   134,   115,   119,    20,    21,    53,    89,   135,   138,
      34,    34,    34,    87,    91,    20,    21,    20,    21,   107,
     132,   128,   119,    66,   147,   113,    90,    98,   131,    84,
      84,    84,   133,   151,    21,    53,   139,   153,   141,    34,
      34,    34,     2,     3,   144,     4,   148,   143,    20,    21,
     145,    20,    21,   -67,    41,    42,    43,   142,    44,    41,
      42,    45,    46,     2,     3,   154,     4,   149,    47,    48,
       5,     0,    49,     2,     3,   102,     4,   108,    50,     0,
      51,    41,    42,   140,   137,    41,    42,    43,   120,    44,
       0,     5,    45,    46,    98,    99,   100,    20,    21,    47,
      48,     5,   -70,    49,     2,     3,     0,     4,     0,    50,
       0,    51,    41,    42,     0,   146,    41,    42,    43,   136,
      44,    95,    96,    45,    46,     2,     3,     0,     4,     0,
      47,    48,     5,     0,    49,     0,     0,    41,    42,    43,
      50,    44,    51,     0,    45,    46,     2,     3,     0,     4,
       0,    47,    48,     5,     0,    49,     0,    20,    21,    20,
      21,    50,   -71,    51,   -72,    45,    46,     0,    20,    21,
      20,    21,    47,   -68,     5,   -69,    49,    41,    42,     0,
     150,     0,    50,     0,    51,     2,     3,    12,     4,    13,
      14,    15,    16,    17,    18,     0,    19,    20,    21,    41,
      42,     0,   152,     0,    22,    23,    24,    98,    99,     0,
      75,    26,     0,     0,    79,    27,    80,    81,    82,    28,
       0,    83,     2,     3,    12,     4,    13,    14,    15,    16,
      17,    18,     0,    19,    20,    21,     0,     0,     0,     0,
       0,    22,    23,    24,     0,     0,     0,    75,    26,     0,
       0,     0,    27,    80,    81,     0,    28,     2,     3,    12,
       4,    13,    14,    15,    16,    17,    18,     0,    19,    20,
      21,     0,     0,     0,     0,     0,    22,    23,    24,     0,
       0,     0,    75,    26,     0,     0,     0,    27,    80,     0,
       0,    28,     2,     3,    12,     4,    13,    14,    15,    16,
      17,    18,     0,    19,    20,    21,     0,     0,     0,     0,
       0,    22,    23,    24,     0,     0,     0,    25,    26,     0,
       0,     0,    27,     0,     0,     0,    28,     2,     3,    12,
       4,    13,    14,    15,    16,    17,    18,     0,    19,    20,
      21,     0,     2,     3,     0,     4,    22,    23,    24,     0,
       0,     0,    75,    26,    41,    42,     0,    27,     0,    20,
      21,    28,     0,     0,     0,     0,     0,    23,    24,    48,
       5,    45,    46,    98,    99,   100,   116,     0,    47,     0,
       0,     0,    49,   117,     0,     0,     0,     0,    50,     0,
      51
};

static const yytype_int16 yycheck[] =
{
       5,     7,    25,    37,     3,     4,    28,     6,    13,    14,
      15,    16,    43,    18,    19,     5,    42,    22,    23,    24,
      25,    43,     0,     1,    80,    81,    82,    27,     1,    29,
       3,     4,    31,     6,    29,    15,    16,    27,    44,    29,
      15,    16,    40,    23,    24,    49,    36,    32,    54,    24,
      20,    55,    32,    32,     3,     4,     5,     6,    31,    39,
      39,    95,    96,    32,    42,    88,    71,    98,    99,   100,
      75,    93,    78,    22,     3,     4,     5,     6,    83,    28,
      85,    32,    87,    88,    15,    16,    92,    77,    39,   123,
      80,    81,    82,    22,    42,    15,    16,    15,    16,    28,
      23,   106,   107,    23,   138,    23,   140,    34,   113,   109,
     110,   111,   117,   147,    16,   121,    32,   151,    21,   109,
     110,   111,     3,     4,    32,     6,    21,   132,    15,    16,
     135,    15,    16,    20,    15,    16,    17,    21,    19,    15,
      16,    22,    23,     3,     4,    18,     6,    32,    29,    30,
      31,    -1,    33,     3,     4,    54,     6,    38,    39,    -1,
      41,    15,    16,   129,    18,    15,    16,    17,    18,    19,
      -1,    31,    22,    23,    34,    35,    36,    15,    16,    29,
      30,    31,    20,    33,     3,     4,    -1,     6,    -1,    39,
      -1,    41,    15,    16,    -1,    18,    15,    16,    17,    18,
      19,    16,    17,    22,    23,     3,     4,    -1,     6,    -1,
      29,    30,    31,    -1,    33,    -1,    -1,    15,    16,    17,
      39,    19,    41,    -1,    22,    23,     3,     4,    -1,     6,
      -1,    29,    30,    31,    -1,    33,    -1,    15,    16,    15,
      16,    39,    20,    41,    20,    22,    23,    -1,    15,    16,
      15,    16,    29,    20,    31,    20,    33,    15,    16,    -1,
      18,    -1,    39,    -1,    41,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    -1,    14,    15,    16,    15,
      16,    -1,    18,    -1,    22,    23,    24,    34,    35,    -1,
      28,    29,    -1,    -1,    32,    33,    34,    35,    36,    37,
      -1,    39,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    -1,    14,    15,    16,    -1,    -1,    -1,    -1,
      -1,    22,    23,    24,    -1,    -1,    -1,    28,    29,    -1,
      -1,    -1,    33,    34,    35,    -1,    37,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    -1,    14,    15,
      16,    -1,    -1,    -1,    -1,    -1,    22,    23,    24,    -1,
      -1,    -1,    28,    29,    -1,    -1,    -1,    33,    34,    -1,
      -1,    37,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    -1,    14,    15,    16,    -1,    -1,    -1,    -1,
      -1,    22,    23,    24,    -1,    -1,    -1,    28,    29,    -1,
      -1,    -1,    33,    -1,    -1,    -1,    37,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    -1,    14,    15,
      16,    -1,     3,     4,    -1,     6,    22,    23,    24,    -1,
      -1,    -1,    28,    29,    15,    16,    -1,    33,    -1,    15,
      16,    37,    -1,    -1,    -1,    -1,    -1,    23,    24,    30,
      31,    22,    23,    34,    35,    36,    32,    -1,    29,    -1,
      -1,    -1,    33,    39,    -1,    -1,    -1,    -1,    39,    -1,
      41
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     6,    31,    44,    45,    47,    48,
      54,    42,     5,     7,     8,     9,    10,    11,    12,    14,
      15,    16,    22,    23,    24,    28,    29,    33,    37,    49,
      50,    51,    52,    53,    54,    58,    60,    61,     0,     1,
      42,    15,    16,    17,    19,    22,    23,    29,    30,    33,
      39,    41,    46,    47,    55,    56,    57,    59,    62,    60,
      60,    60,    60,    60,    60,    60,    23,    60,    24,    60,
      32,    39,    58,    60,    29,    28,    50,    60,    45,    32,
      34,    35,    36,    39,    50,    20,     5,    22,    28,    54,
      62,    42,    45,    55,    47,    16,    17,    57,    34,    35,
      36,    47,    59,    57,    60,    32,    39,    28,    38,    49,
      49,    49,    60,    23,    60,    60,    32,    39,    58,    60,
      18,    45,    62,    62,    55,    55,    55,    32,    60,    40,
      32,    60,    23,    60,    32,    39,    18,    18,    62,    32,
      61,    21,    21,    60,    32,    60,    18,    62,    21,    32,
      18,    62,    18,    62,    18
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    43,    44,    44,    44,    44,    45,    45,    45,    45,
      45,    45,    45,    46,    46,    47,    47,    47,    47,    48,
      48,    48,    48,    48,    48,    48,    48,    49,    49,    49,
      49,    49,    50,    50,    50,    51,    51,    52,    52,    52,
      52,    52,    52,    52,    52,    52,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    52,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    52,    52,    53,    53,    53,
      53,    53,    53,    54,    54,    54,    55,    55,    55,    55,
      56,    56,    57,    57,    57,    57,    57,    57,    58,    58,
      58,    58,    58,    58,    59,    59,    59,    59,    59,    59,
      59,    60,    60,    61,    61,    62,    62
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     3,     2,     2,     1,     2,     3,     3,
       2,     3,     2,     3,     4,     1,     1,     3,     5,     3,
       5,     4,     6,     4,     6,     5,     7,     3,     3,     3,
       2,     1,     2,     1,     1,     3,     5,     1,     2,     1,
       2,     2,     3,     1,     1,     1,     1,     1,     1,     2,
       2,     2,     2,     2,     2,     5,     5,     6,     2,     3,
       2,     1,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     3,     3,     1,
       1,     2,     1,     1,     1,     1,     1,     2,     2,     2,
       1,     2,     2,     1,     1,     3,     4,     5,     6,     7,
       8,     1,     1,     1,     2,     1,     1
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
#line 1495 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 4:
#line 94 "smarts.yy" /* yacc.c:1646  */
    {
  YYACCEPT;
}
#line 1503 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 5:
#line 97 "smarts.yy" /* yacc.c:1646  */
    {
  yyclearin;
  yyerrok;

  yyErrorCleanup(molList);
  YYABORT;
}
#line 1515 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 108 "smarts.yy" /* yacc.c:1646  */
    {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (*molList)[ sz ]->addAtom((yyvsp[0].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
}
#line 1528 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 116 "smarts.yy" /* yacc.c:1646  */
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
#line 1561 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 8:
#line 145 "smarts.yy" /* yacc.c:1646  */
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
#line 1585 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 9:
#line 165 "smarts.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1594 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 10:
#line 170 "smarts.yy" /* yacc.c:1646  */
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
#line 1633 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 11:
#line 205 "smarts.yy" /* yacc.c:1646  */
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
#line 1657 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 12:
#line 225 "smarts.yy" /* yacc.c:1646  */
    {
  RWMol *m1_p = (*molList)[(yyval.moli)],*m2_p=(*molList)[(yyvsp[0].moli)];
  // FIX: handle generic bonds here
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED,Bond::NONE,false,true);
  delete m2_p;
  int sz = molList->size();
  if ( sz==(yyvsp[0].moli)+1) {
    molList->resize( sz-1 );
  }
}
#line 1672 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 13:
#line 239 "smarts.yy" /* yacc.c:1646  */
    { (yyval.moli) = (yyvsp[-1].moli); }
#line 1678 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 14:
#line 240 "smarts.yy" /* yacc.c:1646  */
    {
  // FIX: this needs to handle arbitrary bond_exprs
  (yyval.moli) = (yyvsp[-1].moli);
  int sz     = molList->size();
  (yyvsp[-2].bond)->setOwningMol((*molList)[ sz-1 ]);
  (yyvsp[-2].bond)->setBeginAtomIdx(0);
  (*molList)[ sz-1 ]->setBondBookmark((yyvsp[-2].bond),ci_LEADING_BOND);
}
#line 1691 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 17:
#line 254 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = (yyvsp[-1].atom);
}
#line 1699 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 18:
#line 258 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1708 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 19:
#line 281 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom(1);
}
#line 1716 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 20:
#line 285 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1725 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 21:
#line 289 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1736 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 22:
#line 295 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-4].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-4].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1749 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 23:
#line 304 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1760 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 24:
#line 310 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1773 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 25:
#line 318 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-3].ival));
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1786 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 26:
#line 326 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-5].ival));
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-5].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1801 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 27:
#line 339 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 1812 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 28:
#line 345 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 1823 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 29:
#line 351 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 1834 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 30:
#line 357 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-1].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-1].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  if((yyvsp[0].atom)->getNumExplicitHs()){
    if(!(yyvsp[-1].atom)->getNumExplicitHs()){
      (yyvsp[-1].atom)->setNumExplicitHs((yyvsp[0].atom)->getNumExplicitHs());
      (yyvsp[-1].atom)->setNoImplicit(true);
    } else if((yyvsp[-1].atom)->getNumExplicitHs()!=(yyvsp[0].atom)->getNumExplicitHs()){
      // conflicting queries...
      (yyvsp[-1].atom)->setNumExplicitHs(0);
      (yyvsp[-1].atom)->setNoImplicit(false);
    }
  }
  if((yyvsp[0].atom)->getFormalCharge()){
    if(!(yyvsp[-1].atom)->getFormalCharge()){
      (yyvsp[-1].atom)->setFormalCharge((yyvsp[0].atom)->getFormalCharge());
    } else if((yyvsp[-1].atom)->getFormalCharge()!=(yyvsp[0].atom)->getFormalCharge()){
      // conflicting queries...
      (yyvsp[-1].atom)->setFormalCharge(0);
    }
  }
  delete (yyvsp[0].atom);
}
#line 1862 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 32:
#line 383 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->getQuery()->setNegation(!((yyvsp[0].atom)->getQuery()->getNegation()));
  (yyvsp[0].atom)->setAtomicNum(0);
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[0].atom));
  (yyval.atom) = (yyvsp[0].atom);
}
#line 1873 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 35:
#line 394 "smarts.yy" /* yacc.c:1646  */
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
#line 1895 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 36:
#line 411 "smarts.yy" /* yacc.c:1646  */
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
#line 1921 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 38:
#line 436 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 1931 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 442 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 1941 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 41:
#line 447 "smarts.yy" /* yacc.c:1646  */
    { (yyval.atom) = new QueryAtom((yyvsp[0].ival)); }
#line 1947 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 42:
#line 448 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsotope((yyvsp[-2].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
}
#line 1957 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 49:
#line 459 "smarts.yy" /* yacc.c:1646  */
    {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
}
#line 1965 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 50:
#line 462 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 1973 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 51:
#line 465 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 1981 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 52:
#line 468 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
}
#line 1989 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 53:
#line 471 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
}
#line 1997 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 54:
#line 474 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
}
#line 2005 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 55:
#line 477 "smarts.yy" /* yacc.c:1646  */
    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>((yyvsp[-1].ival),oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2016 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 56:
#line 483 "smarts.yy" /* yacc.c:1646  */
    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>((yyvsp[-2].ival),oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2027 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 57:
#line 489 "smarts.yy" /* yacc.c:1646  */
    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
}
#line 2039 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 58:
#line 496 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)));
  newQ->setIsotope((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2052 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 59:
#line 504 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)));
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomHCountQuery((yyvsp[0].ival)),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2065 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 60:
#line 512 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2076 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 61:
#line 518 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2087 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 62:
#line 524 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  newQ->setFormalCharge((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2098 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 63:
#line 530 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
#line 2109 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 64:
#line 536 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
#line 2120 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 66:
#line 543 "smarts.yy" /* yacc.c:1646  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 2130 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 68:
#line 551 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
}
#line 2138 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 69:
#line 554 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
}
#line 2146 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 70:
#line 557 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
}
#line 2154 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 71:
#line 560 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomRingBondCountQuery(0));
}
#line 2162 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 72:
#line 563 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomImplicitHCountQuery(0));
}
#line 2170 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 73:
#line 569 "smarts.yy" /* yacc.c:1646  */
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
#line 2187 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 74:
#line 581 "smarts.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->expandQuery(makeAtomAromaticQuery(),Queries::COMPOSITE_AND);
}
#line 2197 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 76:
#line 591 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2206 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 77:
#line 595 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
}
#line 2215 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 78:
#line 599 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2224 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 81:
#line 607 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2233 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 83:
#line 615 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
#line 2244 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 84:
#line 621 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
#line 2255 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 85:
#line 627 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
#line 2266 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 86:
#line 633 "smarts.yy" /* yacc.c:1646  */
    {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
#line 2276 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 87:
#line 638 "smarts.yy" /* yacc.c:1646  */
    {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
#line 2285 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 88:
#line 645 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=2; }
#line 2291 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 89:
#line 646 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=(yyvsp[0].ival); }
#line 2297 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 90:
#line 647 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=1; }
#line 2303 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 91:
#line 648 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=-2; }
#line 2309 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 92:
#line 649 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=-(yyvsp[0].ival); }
#line 2315 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 93:
#line 650 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival)=-1; }
#line 2321 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 95:
#line 655 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 2327 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 96:
#line 656 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 2333 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 97:
#line 657 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2339 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 98:
#line 658 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2345 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 99:
#line 659 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2351 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 100:
#line 660 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2357 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;

  case 104:
#line 671 "smarts.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
#line 2363 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
    break;


#line 2367 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1646  */
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
#line 678 "smarts.yy" /* yacc.c:1906  */

