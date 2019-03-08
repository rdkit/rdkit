/* A Bison parser, made by GNU Bison 3.1.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018 Free Software Foundation, Inc.

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
#define YYBISON_VERSION "3.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 2

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

extern int yysmarts_lex(YYSTYPE *,void *, int &);

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
void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                RDKit::Atom* &lastAtom,
                RDKit::Bond* &lastBond,
		void *scanner,int start_token, const char * msg )
{
  RDUNUSED_PARAM(input);
  RDUNUSED_PARAM(lastAtom);
  RDUNUSED_PARAM(lastBond);
  RDUNUSED_PARAM(scanner);
  RDUNUSED_PARAM(start_token);
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error: " << msg << " while parsing: " << input << std::endl;
}

#line 118 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:339  */

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
#ifndef YY_YYSMARTS_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
# define YY_YYSMARTS_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED
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
    START_MOL = 258,
    START_ATOM = 259,
    START_BOND = 260,
    AROMATIC_ATOM_TOKEN = 261,
    ORGANIC_ATOM_TOKEN = 262,
    ATOM_TOKEN = 263,
    SIMPLE_ATOM_QUERY_TOKEN = 264,
    COMPLEX_ATOM_QUERY_TOKEN = 265,
    RINGSIZE_ATOM_QUERY_TOKEN = 266,
    RINGBOND_ATOM_QUERY_TOKEN = 267,
    IMPLICIT_H_ATOM_QUERY_TOKEN = 268,
    HYB_TOKEN = 269,
    HETERONEIGHBOR_ATOM_QUERY_TOKEN = 270,
    ALIPHATIC = 271,
    ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN = 272,
    ZERO_TOKEN = 273,
    NONZERO_DIGIT_TOKEN = 274,
    GROUP_OPEN_TOKEN = 275,
    GROUP_CLOSE_TOKEN = 276,
    SEPARATOR_TOKEN = 277,
    RANGE_OPEN_TOKEN = 278,
    RANGE_CLOSE_TOKEN = 279,
    HASH_TOKEN = 280,
    MINUS_TOKEN = 281,
    PLUS_TOKEN = 282,
    CHIRAL_MARKER_TOKEN = 283,
    CHI_CLASS_TOKEN = 284,
    CHI_CLASS_OH_TOKEN = 285,
    H_TOKEN = 286,
    AT_TOKEN = 287,
    PERCENT_TOKEN = 288,
    ATOM_OPEN_TOKEN = 289,
    ATOM_CLOSE_TOKEN = 290,
    NOT_TOKEN = 291,
    AND_TOKEN = 292,
    OR_TOKEN = 293,
    SEMI_TOKEN = 294,
    BEGIN_RECURSE = 295,
    END_RECURSE = 296,
    COLON_TOKEN = 297,
    UNDERSCORE_TOKEN = 298,
    BOND_TOKEN = 299,
    EOS_TOKEN = 300
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 62 "smarts.yy" /* yacc.c:355  */

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;

#line 211 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token);
/* "%code provides" blocks.  */
#line 57 "smarts.yy" /* yacc.c:355  */

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)

#line 228 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:355  */

#endif /* !YY_YYSMARTS_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 234 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:358  */

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
typedef unsigned short yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short yytype_int16;
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
#  define YYSIZE_T unsigned
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

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
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
#define YYFINAL  55
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   639

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  46
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  21
/* YYNRULES -- Number of rules.  */
#define YYNRULES  116
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  173

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   300

#define YYTRANSLATE(YYX)                                                \
  ((unsigned) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

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
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   106,   106,   109,   113,   116,   119,   123,   127,   130,
     135,   138,   146,   147,   148,   149,   157,   165,   190,   210,
     215,   246,   266,   279,   280,   291,   292,   293,   297,   320,
     324,   329,   335,   344,   350,   358,   366,   379,   385,   392,
     398,   421,   424,   430,   431,   435,   452,   476,   477,   482,
     483,   488,   489,   494,   495,   496,   497,   498,   499,   500,
     503,   506,   509,   512,   515,   518,   524,   530,   537,   545,
     553,   559,   565,   571,   577,   583,   584,   591,   592,   595,
     598,   601,   604,   610,   622,   627,   632,   636,   640,   644,
     647,   648,   655,   656,   662,   668,   674,   679,   686,   687,
     688,   689,   690,   691,   695,   696,   697,   698,   699,   700,
     701,   706,   707,   711,   712,   715,   716
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "START_MOL", "START_ATOM", "START_BOND",
  "AROMATIC_ATOM_TOKEN", "ORGANIC_ATOM_TOKEN", "ATOM_TOKEN",
  "SIMPLE_ATOM_QUERY_TOKEN", "COMPLEX_ATOM_QUERY_TOKEN",
  "RINGSIZE_ATOM_QUERY_TOKEN", "RINGBOND_ATOM_QUERY_TOKEN",
  "IMPLICIT_H_ATOM_QUERY_TOKEN", "HYB_TOKEN",
  "HETERONEIGHBOR_ATOM_QUERY_TOKEN", "ALIPHATIC",
  "ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN", "ZERO_TOKEN",
  "NONZERO_DIGIT_TOKEN", "GROUP_OPEN_TOKEN", "GROUP_CLOSE_TOKEN",
  "SEPARATOR_TOKEN", "RANGE_OPEN_TOKEN", "RANGE_CLOSE_TOKEN", "HASH_TOKEN",
  "MINUS_TOKEN", "PLUS_TOKEN", "CHIRAL_MARKER_TOKEN", "CHI_CLASS_TOKEN",
  "CHI_CLASS_OH_TOKEN", "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN",
  "ATOM_OPEN_TOKEN", "ATOM_CLOSE_TOKEN", "NOT_TOKEN", "AND_TOKEN",
  "OR_TOKEN", "SEMI_TOKEN", "BEGIN_RECURSE", "END_RECURSE", "COLON_TOKEN",
  "UNDERSCORE_TOKEN", "BOND_TOKEN", "EOS_TOKEN", "$accept", "meta_start",
  "bad_atom_def", "mol", "branch", "atomd", "hydrogen_atom", "atom_expr",
  "point_query", "recursive_query", "atom_query", "possible_range_query",
  "simple_atom", "bond_expr", "bond_query", "bondd", "charge_spec",
  "ring_number", "number", "nonzero_number", "digit", YY_NULLPTR
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
     295,   296,   297,   298,   299,   300
};
# endif

#define YYPACT_NINF -56

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-56)))

#define YYTABLE_NINF -83

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     106,   -13,    30,   322,   176,    23,   -56,   -56,   -56,   -56,
     538,   225,   -56,   -56,   -56,   -56,    54,    79,   104,   112,
     -56,   147,   245,   -56,   -56,   154,     3,     7,   154,    33,
     359,   396,   573,    30,   396,   -56,    12,   433,   -56,   -56,
     -56,    32,    46,   -56,   132,   178,   -56,   -56,   -56,   176,
     -56,   -56,     4,   176,   -56,   -56,    51,   -56,   593,   285,
     -56,   146,   -56,   -56,   246,    30,   184,   -56,   -56,   588,
     -56,   -56,   -56,   -56,   -56,   -56,   -56,   -56,   -56,   -56,
     -56,   -56,   -56,   -56,   -56,   396,   -56,   285,   -56,   -56,
     142,   -56,   -56,   573,   573,   573,   -56,    43,   -56,   154,
     154,   -56,   -56,   -56,   176,   176,   176,   -56,   -56,   -56,
     -56,   154,    -4,   -56,   154,   597,   173,    45,   -56,   178,
     178,   -56,   -56,    27,   573,   503,   468,   154,    52,   -56,
     -56,   -56,    83,   190,    53,   -56,   154,    89,   -56,   154,
       5,   -56,   204,   -56,   115,   109,   118,    26,   -56,   121,
     -56,   124,   -56,   154,   -56,   -56,   169,   178,   -56,   -56,
     126,   -56,   -56,   134,   -56,   255,   -56,   -56,   -56,   265,
     -56,   160,   -56
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     5,     8,     0,    11,    84,    83,    85,
       0,     2,    16,    26,    25,    49,    53,    56,    57,    58,
      75,    54,    55,   111,   113,     0,   103,   100,    71,    74,
       0,     0,     0,     0,     0,     4,     0,    15,    41,    43,
      44,     0,    47,    72,    76,   112,    94,    93,    96,     0,
      95,    92,     7,    89,    90,     1,     0,    10,    71,     0,
      47,    76,   116,   115,     0,     0,     0,    22,    17,     0,
      20,   104,    59,    62,    63,    64,    60,    61,    51,   101,
     102,    98,    99,    70,    73,     0,    12,    15,    13,    42,
       0,    14,     3,     0,     0,     0,    40,     0,    50,     0,
      68,    48,   114,    97,     0,     0,     0,     6,    91,     9,
      29,     0,     0,    27,     0,    68,     0,     0,    19,     0,
       0,    18,    21,    45,    37,    38,    39,     0,     0,    52,
      69,    86,    87,    88,     0,    33,     0,     0,    31,     0,
       0,    23,     0,   105,     0,     0,     0,     0,    30,     0,
      28,     0,    35,     0,    24,   106,     0,    46,    65,    66,
       0,    34,    32,     0,   107,     0,    67,    36,   108,     0,
     109,     0,   110
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -56,   -56,    25,   -31,   -56,    16,   -56,     0,    21,   -56,
     -56,   -56,     2,    10,   -56,    50,   -55,   116,   -10,    44,
     -44
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     5,    86,    11,    67,    12,    13,    37,    38,    39,
      40,    41,    60,    69,    53,    54,    43,    70,    44,    45,
      71
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      61,   102,    90,   112,    14,    42,    72,    73,    74,    75,
      59,    76,    77,    14,    52,    78,    80,    82,    83,    36,
      61,    23,    24,    55,    56,    23,    24,    68,    35,    79,
      87,   135,     6,   116,    81,    14,     7,     8,   136,     9,
     152,   104,   105,   106,    23,    24,   101,   153,    83,   107,
     159,     7,     8,    89,     9,    97,    88,    92,    96,    91,
     140,    23,    24,   101,    10,    84,    14,    14,    57,   127,
     145,    14,    23,    24,   117,   143,   144,   -77,   147,    10,
      96,   118,   104,   105,   106,   121,   142,   128,   148,   129,
     130,   -25,    14,   124,   125,   126,   109,    23,    24,   103,
     156,   134,   -80,   108,   137,   130,    68,     1,    96,     2,
       3,     4,   165,   102,   131,   132,   133,   146,    14,    14,
     104,   169,    23,    24,   150,   171,   149,   -81,    24,   151,
      23,    24,    68,    62,    63,   -82,   155,   160,     7,     8,
      98,     9,   158,   163,    14,    96,    96,    96,     7,     8,
     166,     9,     7,     8,    98,     9,   161,    99,    68,   162,
      62,    63,    64,   100,    65,    23,    24,    46,    47,   167,
     -78,    99,    23,    24,    48,    66,    10,   115,    49,     7,
       8,   172,     9,   123,    50,   122,    51,    62,    63,   157,
     164,    62,    63,    64,   141,    65,    62,    63,    46,    47,
       0,    46,    47,   119,   120,    48,    66,    10,    48,    49,
       7,     8,    49,     9,     0,    50,     0,    51,    50,     0,
      51,     0,    62,    63,    64,   154,    65,   104,   105,    46,
      47,     7,     8,     0,     9,     0,    48,    66,    10,     0,
      49,     0,     0,    62,    63,    64,    50,    65,    51,     0,
      46,    47,     7,     8,     0,     9,     0,    48,    66,    10,
       0,    49,     0,    23,    24,     0,     0,    50,   -79,    51,
       0,    46,    47,    62,    63,     0,   168,     0,    48,     0,
      10,     0,    49,    62,    63,     0,   170,     0,    50,     0,
      51,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,     0,    22,    23,    24,     0,     0,     0,     0,     0,
      25,    26,    27,     0,     0,     0,    28,    29,     0,     0,
     113,    32,    93,    94,    95,    33,     0,   114,     7,     8,
      15,     9,    16,    17,    18,    19,    20,    21,     0,    22,
      23,    24,     0,     0,     0,     0,     0,    25,    26,    27,
       0,     0,     0,    28,    29,     0,    30,    31,    32,     0,
       0,     0,    33,     0,    34,     7,     8,    15,     9,    16,
      17,    18,    19,    20,    21,     0,    22,    23,    24,     0,
       0,     0,     0,     0,    25,    26,    27,     0,     0,     0,
      58,    29,     0,    85,    31,    32,     0,     0,     0,    33,
       0,    34,     7,     8,    15,     9,    16,    17,    18,    19,
      20,    21,     0,    22,    23,    24,     0,     0,     0,     0,
       0,    25,    26,    27,     0,     0,     0,    28,    29,     0,
      85,    31,    32,     0,     0,     0,    33,     0,    34,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,     0,
      22,    23,    24,     0,     0,     0,     0,     0,    25,    26,
      27,     0,     0,     0,    28,    29,     0,     0,     0,    32,
      93,    94,    95,    33,     7,     8,    15,     9,    16,    17,
      18,    19,    20,    21,     0,    22,    23,    24,     0,     0,
       0,     0,     0,    25,    26,    27,     0,     0,     0,    28,
      29,     0,     0,     0,    32,    93,    94,     0,    33,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,     0,
      22,    23,    24,     0,     0,     0,     0,     0,    25,    26,
      27,     0,     0,     0,    28,    29,     0,     0,     0,    32,
      93,     0,     0,    33,     7,     8,    15,     9,    16,    17,
      18,    19,    20,    21,     0,    22,    23,    24,     0,     0,
       0,     0,     0,    25,    26,    27,     0,     0,     0,    58,
      29,     0,     0,     0,    32,     0,     0,     0,    33,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,     0,
      22,    23,    24,     0,     7,     8,     0,     9,    25,    26,
      27,     0,     0,     0,    28,    29,    62,    63,     0,    32,
       0,    23,    24,    33,     0,    23,    24,     0,     0,    26,
      27,    66,    10,    26,    27,   104,   105,   106,   110,     0,
       0,     0,   138,     0,     0,   111,     0,     0,     0,   139
};

static const yytype_int16 yycheck[] =
{
      10,    45,    33,    58,     2,     3,    16,    17,    18,    19,
      10,    21,    22,    11,     4,    25,    26,    27,    28,     3,
      30,    18,    19,     0,     1,    18,    19,    11,     3,    26,
      30,    35,    45,    64,    27,    33,     6,     7,    42,     9,
      35,    37,    38,    39,    18,    19,    44,    42,    58,    45,
      24,     6,     7,    32,     9,    23,    31,    45,    37,    34,
     115,    18,    19,    61,    34,    32,    64,    65,    45,    26,
      43,    69,    18,    19,    64,   119,   120,    23,    26,    34,
      59,    65,    37,    38,    39,    69,   117,    97,    35,    99,
     100,    45,    90,    93,    94,    95,    45,    18,    19,    49,
     144,   111,    23,    53,   114,   115,    90,     1,    87,     3,
       4,     5,   156,   157,   104,   105,   106,   127,   116,   117,
      37,   165,    18,    19,    35,   169,   136,    23,    19,   139,
      18,    19,   116,    18,    19,    23,    21,   147,     6,     7,
       8,     9,    24,   153,   142,   124,   125,   126,     6,     7,
      24,     9,     6,     7,     8,     9,    35,    25,   142,    35,
      18,    19,    20,    31,    22,    18,    19,    25,    26,    35,
      23,    25,    18,    19,    32,    33,    34,    31,    36,     6,
       7,    21,     9,    41,    42,    69,    44,    18,    19,   145,
      21,    18,    19,    20,    21,    22,    18,    19,    25,    26,
      -1,    25,    26,    19,    20,    32,    33,    34,    32,    36,
       6,     7,    36,     9,    -1,    42,    -1,    44,    42,    -1,
      44,    -1,    18,    19,    20,    21,    22,    37,    38,    25,
      26,     6,     7,    -1,     9,    -1,    32,    33,    34,    -1,
      36,    -1,    -1,    18,    19,    20,    42,    22,    44,    -1,
      25,    26,     6,     7,    -1,     9,    -1,    32,    33,    34,
      -1,    36,    -1,    18,    19,    -1,    -1,    42,    23,    44,
      -1,    25,    26,    18,    19,    -1,    21,    -1,    32,    -1,
      34,    -1,    36,    18,    19,    -1,    21,    -1,    42,    -1,
      44,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    -1,    17,    18,    19,    -1,    -1,    -1,    -1,    -1,
      25,    26,    27,    -1,    -1,    -1,    31,    32,    -1,    -1,
      35,    36,    37,    38,    39,    40,    -1,    42,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    -1,    17,
      18,    19,    -1,    -1,    -1,    -1,    -1,    25,    26,    27,
      -1,    -1,    -1,    31,    32,    -1,    34,    35,    36,    -1,
      -1,    -1,    40,    -1,    42,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    -1,    17,    18,    19,    -1,
      -1,    -1,    -1,    -1,    25,    26,    27,    -1,    -1,    -1,
      31,    32,    -1,    34,    35,    36,    -1,    -1,    -1,    40,
      -1,    42,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    -1,    17,    18,    19,    -1,    -1,    -1,    -1,
      -1,    25,    26,    27,    -1,    -1,    -1,    31,    32,    -1,
      34,    35,    36,    -1,    -1,    -1,    40,    -1,    42,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    -1,
      17,    18,    19,    -1,    -1,    -1,    -1,    -1,    25,    26,
      27,    -1,    -1,    -1,    31,    32,    -1,    -1,    -1,    36,
      37,    38,    39,    40,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    -1,    17,    18,    19,    -1,    -1,
      -1,    -1,    -1,    25,    26,    27,    -1,    -1,    -1,    31,
      32,    -1,    -1,    -1,    36,    37,    38,    -1,    40,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    -1,
      17,    18,    19,    -1,    -1,    -1,    -1,    -1,    25,    26,
      27,    -1,    -1,    -1,    31,    32,    -1,    -1,    -1,    36,
      37,    -1,    -1,    40,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    -1,    17,    18,    19,    -1,    -1,
      -1,    -1,    -1,    25,    26,    27,    -1,    -1,    -1,    31,
      32,    -1,    -1,    -1,    36,    -1,    -1,    -1,    40,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    -1,
      17,    18,    19,    -1,     6,     7,    -1,     9,    25,    26,
      27,    -1,    -1,    -1,    31,    32,    18,    19,    -1,    36,
      -1,    18,    19,    40,    -1,    18,    19,    -1,    -1,    26,
      27,    33,    34,    26,    27,    37,    38,    39,    35,    -1,
      -1,    -1,    35,    -1,    -1,    42,    -1,    -1,    -1,    42
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,    47,    45,     6,     7,     9,
      34,    49,    51,    52,    58,     8,    10,    11,    12,    13,
      14,    15,    17,    18,    19,    25,    26,    27,    31,    32,
      34,    35,    36,    40,    42,    48,    51,    53,    54,    55,
      56,    57,    58,    62,    64,    65,    25,    26,    32,    36,
      42,    44,    59,    60,    61,     0,     1,    45,    31,    53,
      58,    64,    18,    19,    20,    22,    33,    50,    51,    59,
      63,    66,    64,    64,    64,    64,    64,    64,    64,    26,
      64,    27,    64,    64,    32,    34,    48,    53,    48,    54,
      49,    48,    45,    37,    38,    39,    54,    23,     8,    25,
      31,    58,    66,    61,    37,    38,    39,    45,    61,    45,
      35,    42,    62,    35,    42,    31,    49,    59,    51,    19,
      20,    51,    63,    41,    53,    53,    53,    26,    64,    64,
      64,    59,    59,    59,    64,    35,    42,    64,    35,    42,
      62,    21,    49,    66,    66,    43,    64,    26,    35,    64,
      35,    64,    35,    42,    21,    21,    66,    65,    24,    24,
      64,    35,    35,    64,    21,    66,    24,    35,    21,    66,
      21,    66,    21
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    46,    47,    47,    47,    47,    47,    47,    47,    47,
      47,    47,    48,    48,    48,    48,    49,    49,    49,    49,
      49,    49,    49,    50,    50,    51,    51,    51,    51,    52,
      52,    52,    52,    52,    52,    52,    52,    53,    53,    53,
      53,    53,    54,    54,    54,    55,    55,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    57,    57,    57,
      57,    57,    57,    58,    58,    58,    59,    59,    59,    59,
      60,    60,    61,    61,    61,    61,    61,    61,    62,    62,
      62,    62,    62,    62,    63,    63,    63,    63,    63,    63,
      63,    64,    64,    65,    65,    66,    66
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     3,     2,     1,     3,     2,     1,     3,
       2,     2,     2,     2,     2,     1,     1,     2,     3,     3,
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
      yyerror (input, molList, lastAtom, lastBond, scanner, start_token, YY_("syntax error: cannot back up")); \
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
                  Type, Value, input, molList, lastAtom, lastBond, scanner, start_token); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (scanner);
  YYUSE (start_token);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, lastAtom, lastBond, scanner, start_token);
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
{
  unsigned long yylno = yyrline[yyrule];
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
                                              , input, molList, lastAtom, lastBond, scanner, start_token);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, molList, lastAtom, lastBond, scanner, start_token); \
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
    default: /* Avoid compiler warnings. */
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (scanner);
  YYUSE (start_token);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yytype)
    {
          case 8: /* ATOM_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1298 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 9: /* SIMPLE_ATOM_QUERY_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1304 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 10: /* COMPLEX_ATOM_QUERY_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1310 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 11: /* RINGSIZE_ATOM_QUERY_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1316 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 12: /* RINGBOND_ATOM_QUERY_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1322 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 13: /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1328 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 14: /* HYB_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1334 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 15: /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1340 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 16: /* ALIPHATIC  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1346 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 17: /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1352 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 44: /* BOND_TOKEN  */
#line 98 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).bond); }
#line 1358 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 51: /* atomd  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1364 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 52: /* hydrogen_atom  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1370 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 53: /* atom_expr  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1376 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 54: /* point_query  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1382 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 55: /* recursive_query  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1388 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 56: /* atom_query  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1394 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 57: /* possible_range_query  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1400 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 58: /* simple_atom  */
#line 97 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).atom); }
#line 1406 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 59: /* bond_expr  */
#line 98 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).bond); }
#line 1412 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 60: /* bond_query  */
#line 98 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).bond); }
#line 1418 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;

    case 61: /* bondd  */
#line 98 "smarts.yy" /* yacc.c:1258  */
      { delete ((*yyvaluep).bond); }
#line 1424 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1258  */
        break;


      default:
        break;
    }
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
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
                  (unsigned long) yystacksize));

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
      yychar = yylex (&yylval, scanner, start_token);
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
        case 2:
#line 106 "smarts.yy" /* yacc.c:1651  */
    {
// the molList has already been updated, no need to do anything
}
#line 1694 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 3:
#line 109 "smarts.yy" /* yacc.c:1651  */
    {
  lastAtom = (yyvsp[-1].atom);
  YYACCEPT;
}
#line 1703 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 4:
#line 113 "smarts.yy" /* yacc.c:1651  */
    {
  YYABORT;
}
#line 1711 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 5:
#line 116 "smarts.yy" /* yacc.c:1651  */
    {
  YYABORT;
}
#line 1719 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 6:
#line 119 "smarts.yy" /* yacc.c:1651  */
    {
  lastBond = (yyvsp[-1].bond);
  YYACCEPT;
}
#line 1728 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 7:
#line 123 "smarts.yy" /* yacc.c:1651  */
    {
  delete (yyvsp[0].bond);
  YYABORT;
}
#line 1737 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 8:
#line 127 "smarts.yy" /* yacc.c:1651  */
    {
  YYABORT;
}
#line 1745 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 9:
#line 130 "smarts.yy" /* yacc.c:1651  */
    {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1755 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 10:
#line 135 "smarts.yy" /* yacc.c:1651  */
    {
  YYACCEPT;
}
#line 1763 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 11:
#line 138 "smarts.yy" /* yacc.c:1651  */
    {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1773 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 15:
#line 149 "smarts.yy" /* yacc.c:1651  */
    {
  delete (yyvsp[0].atom);
  YYABORT;
}
#line 1782 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 16:
#line 157 "smarts.yy" /* yacc.c:1651  */
    {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (*molList)[ sz ]->addAtom((yyvsp[0].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
}
#line 1795 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 17:
#line 165 "smarts.yy" /* yacc.c:1651  */
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
    newB->setQuery(makeSingleOrAromaticBondQuery());
  } else {
    newB = new QueryBond(Bond::AROMATIC);
    newB->setQuery(makeSingleOrAromaticBondQuery());
  }
  newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  mp->addBond(newB);
  delete newB;
  //delete $2;
}
#line 1824 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 18:
#line 190 "smarts.yy" /* yacc.c:1651  */
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
#line 1848 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 19:
#line 210 "smarts.yy" /* yacc.c:1651  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1857 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 20:
#line 215 "smarts.yy" /* yacc.c:1651  */
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();

  // this is a bit of a hack to try and get nicer "SMILES" from
  // a SMARTS molecule:
  QueryBond * newB;
  if(!atom->getIsAromatic()){
    newB = new QueryBond(Bond::SINGLE);
    newB->setQuery(makeSingleOrAromaticBondQuery());
  } else {
    newB = new QueryBond(Bond::AROMATIC);
    newB->setQuery(makeSingleOrAromaticBondQuery());
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
#line 1892 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 21:
#line 246 "smarts.yy" /* yacc.c:1651  */
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
#line 1916 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 22:
#line 266 "smarts.yy" /* yacc.c:1651  */
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
#line 1931 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 23:
#line 279 "smarts.yy" /* yacc.c:1651  */
    { (yyval.moli) = (yyvsp[-1].moli); }
#line 1937 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 24:
#line 280 "smarts.yy" /* yacc.c:1651  */
    {
  // FIX: this needs to handle arbitrary bond_exprs
  (yyval.moli) = (yyvsp[-1].moli);
  int sz     = molList->size();
  (yyvsp[-2].bond)->setOwningMol((*molList)[ sz-1 ]);
  (yyvsp[-2].bond)->setBeginAtomIdx(0);
  (*molList)[ sz-1 ]->setBondBookmark((yyvsp[-2].bond),ci_LEADING_BOND);
}
#line 1950 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 27:
#line 294 "smarts.yy" /* yacc.c:1651  */
    {
  (yyval.atom) = (yyvsp[-1].atom);
}
#line 1958 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 28:
#line 298 "smarts.yy" /* yacc.c:1651  */
    {
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1967 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 29:
#line 321 "smarts.yy" /* yacc.c:1651  */
    {
  (yyval.atom) = new QueryAtom(1);
}
#line 1975 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 30:
#line 325 "smarts.yy" /* yacc.c:1651  */
    {
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1984 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 31:
#line 329 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1995 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 32:
#line 335 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-4].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-4].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2008 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 33:
#line 344 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2019 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 34:
#line 350 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2032 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 35:
#line 358 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-3].ival));
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2045 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 36:
#line 366 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-5].ival));
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-5].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2060 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 37:
#line 379 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 2071 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 38:
#line 385 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  (yyvsp[-2].atom)->setAtomicNum(0);
  delete (yyvsp[0].atom);
}
#line 2083 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 39:
#line 392 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 2094 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 40:
#line 398 "smarts.yy" /* yacc.c:1651  */
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
#line 2122 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 42:
#line 424 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->getQuery()->setNegation(!((yyvsp[0].atom)->getQuery()->getNegation()));
  (yyvsp[0].atom)->setAtomicNum(0);
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[0].atom));
  (yyval.atom) = (yyvsp[0].atom);
}
#line 2133 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 45:
#line 435 "smarts.yy" /* yacc.c:1651  */
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
#line 2155 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 46:
#line 452 "smarts.yy" /* yacc.c:1651  */
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
#line 2181 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 48:
#line 477 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2191 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 50:
#line 483 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2201 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 51:
#line 488 "smarts.yy" /* yacc.c:1651  */
    { (yyval.atom) = new QueryAtom((yyvsp[0].ival)); }
#line 2207 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 52:
#line 489 "smarts.yy" /* yacc.c:1651  */
    {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsotope((yyvsp[-2].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
}
#line 2217 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 59:
#line 500 "smarts.yy" /* yacc.c:1651  */
    {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
}
#line 2225 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 60:
#line 503 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2233 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 61:
#line 506 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2241 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 62:
#line 509 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
}
#line 2249 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 63:
#line 512 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
}
#line 2257 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 64:
#line 515 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
}
#line 2265 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 65:
#line 518 "smarts.yy" /* yacc.c:1651  */
    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>((yyvsp[-1].ival),oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2276 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 66:
#line 524 "smarts.yy" /* yacc.c:1651  */
    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>((yyvsp[-2].ival),oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2287 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 67:
#line 530 "smarts.yy" /* yacc.c:1651  */
    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
}
#line 2299 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 68:
#line 537 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)));
  newQ->setIsotope((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2312 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 69:
#line 545 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)));
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomHCountQuery((yyvsp[0].ival)),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2325 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 70:
#line 553 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2336 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 71:
#line 559 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2347 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 72:
#line 565 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  newQ->setFormalCharge((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2358 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 73:
#line 571 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
#line 2369 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 74:
#line 577 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
#line 2380 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 76:
#line 584 "smarts.yy" /* yacc.c:1651  */
    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 2390 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 78:
#line 592 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
}
#line 2398 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 79:
#line 595 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
}
#line 2406 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 80:
#line 598 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
}
#line 2414 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 81:
#line 601 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomRingBondCountQuery(0));
}
#line 2422 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 82:
#line 604 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].atom)->setQuery(makeAtomImplicitHCountQuery(0));
}
#line 2430 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 83:
#line 610 "smarts.yy" /* yacc.c:1651  */
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
  (yyval.atom)->setQuery(makeAtomTypeQuery((yyvsp[0].ival),false));
}
#line 2447 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 84:
#line 622 "smarts.yy" /* yacc.c:1651  */
    {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->setQuery(makeAtomTypeQuery((yyvsp[0].ival),true));
}
#line 2457 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 86:
#line 632 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2466 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 87:
#line 636 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
}
#line 2475 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 88:
#line 640 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2484 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 91:
#line 648 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2493 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 93:
#line 656 "smarts.yy" /* yacc.c:1651  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
#line 2504 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 94:
#line 662 "smarts.yy" /* yacc.c:1651  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
#line 2515 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 95:
#line 668 "smarts.yy" /* yacc.c:1651  */
    {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
#line 2526 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 96:
#line 674 "smarts.yy" /* yacc.c:1651  */
    {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
#line 2536 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 97:
#line 679 "smarts.yy" /* yacc.c:1651  */
    {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
#line 2545 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 98:
#line 686 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival)=2; }
#line 2551 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 99:
#line 687 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival)=(yyvsp[0].ival); }
#line 2557 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 100:
#line 688 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival)=1; }
#line 2563 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 101:
#line 689 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival)=-2; }
#line 2569 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 102:
#line 690 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival)=-(yyvsp[0].ival); }
#line 2575 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 103:
#line 691 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival)=-1; }
#line 2581 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 105:
#line 696 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 2587 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 106:
#line 697 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 2593 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 107:
#line 698 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2599 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 108:
#line 699 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2605 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 109:
#line 700 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2611 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 110:
#line 701 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2617 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;

  case 114:
#line 712 "smarts.yy" /* yacc.c:1651  */
    { (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
#line 2623 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
    break;


#line 2627 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp" /* yacc.c:1651  */
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
      yyerror (input, molList, lastAtom, lastBond, scanner, start_token, YY_("syntax error"));
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
        yyerror (input, molList, lastAtom, lastBond, scanner, start_token, yymsgp);
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
                      yytoken, &yylval, input, molList, lastAtom, lastBond, scanner, start_token);
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
                  yystos[yystate], yyvsp, input, molList, lastAtom, lastBond, scanner, start_token);
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
  yyerror (input, molList, lastAtom, lastBond, scanner, start_token, YY_("memory exhausted"));
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
                  yytoken, &yylval, input, molList, lastAtom, lastBond, scanner, start_token);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, input, molList, lastAtom, lastBond, scanner, start_token);
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
#line 719 "smarts.yy" /* yacc.c:1910  */

