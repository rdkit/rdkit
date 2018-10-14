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
#define YYPURE 2

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         yysmiles_parse
#define yylex           yysmiles_lex
#define yyerror         yysmiles_error
#define yydebug         yysmiles_debug
#define yynerrs         yysmiles_nerrs


/* Copy the first part of user declarations.  */
#line 1 "smiles.yy" /* yacc.c:339  */


  // $Id$
  //
  //  Copyright (C) 2001-2016 Randal Henne, Greg Landrum and Rational Discovery LLC
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

#define YYDEBUG 1
#include "smiles.tab.hpp"

extern int yysmiles_lex(YYSTYPE *,void *,int &);


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
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                RDKit::Atom* &lastAtom,
                RDKit::Bond* &lastBond,
                std::list<unsigned int> *branchPoints,
		void *scanner,int start_token, const char * msg )
{
  yyErrorCleanup(ms);
  throw RDKit::SmilesParseException(msg);
}

void
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::list<unsigned int> *branchPoints,
		void *scanner,int start_token, const char * msg )
{
  yyErrorCleanup(ms);
  throw RDKit::SmilesParseException(msg);
}



#line 133 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:339  */

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
   by #include "smiles.tab.hpp".  */
#ifndef YY_YYSMILES_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED
# define YY_YYSMILES_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yysmiles_debug;
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
    ATOM_TOKEN = 262,
    ORGANIC_ATOM_TOKEN = 263,
    NONZERO_DIGIT_TOKEN = 264,
    ZERO_TOKEN = 265,
    GROUP_OPEN_TOKEN = 266,
    GROUP_CLOSE_TOKEN = 267,
    SEPARATOR_TOKEN = 268,
    LOOP_CONNECTOR_TOKEN = 269,
    MINUS_TOKEN = 270,
    PLUS_TOKEN = 271,
    CHIRAL_MARKER_TOKEN = 272,
    CHI_CLASS_TOKEN = 273,
    CHI_CLASS_OH_TOKEN = 274,
    H_TOKEN = 275,
    AT_TOKEN = 276,
    PERCENT_TOKEN = 277,
    COLON_TOKEN = 278,
    HASH_TOKEN = 279,
    BOND_TOKEN = 280,
    ATOM_OPEN_TOKEN = 281,
    ATOM_CLOSE_TOKEN = 282,
    EOS_TOKEN = 283
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 78 "smiles.yy" /* yacc.c:355  */

  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  int                      ival;

#line 209 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmiles_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, std::list<unsigned int> *branchPoints, void *scanner, int& start_token);
/* "%code provides" blocks.  */
#line 73 "smiles.yy" /* yacc.c:355  */

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)

#line 226 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:355  */

#endif /* !YY_YYSMILES_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 232 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:358  */

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
#define YYFINAL  17
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   115

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  29
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  14
/* YYNRULES -- Number of rules.  */
#define YYNRULES  64
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  92

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   283

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
      25,    26,    27,    28
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   103,   103,   106,   109,   112,   118,   121,   132,   143,
     153,   173,   181,   187,   206,   224,   240,   250,   270,   278,
     287,   288,   296,   297,   304,   312,   313,   314,   315,   316,
     317,   318,   322,   323,   324,   325,   326,   327,   328,   329,
     330,   334,   335,   336,   340,   341,   342,   343,   344,   345,
     349,   350,   354,   355,   356,   357,   358,   359,   360,   364,
     365,   369,   370,   373,   374
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "START_MOL", "START_ATOM", "START_BOND",
  "AROMATIC_ATOM_TOKEN", "ATOM_TOKEN", "ORGANIC_ATOM_TOKEN",
  "NONZERO_DIGIT_TOKEN", "ZERO_TOKEN", "GROUP_OPEN_TOKEN",
  "GROUP_CLOSE_TOKEN", "SEPARATOR_TOKEN", "LOOP_CONNECTOR_TOKEN",
  "MINUS_TOKEN", "PLUS_TOKEN", "CHIRAL_MARKER_TOKEN", "CHI_CLASS_TOKEN",
  "CHI_CLASS_OH_TOKEN", "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN",
  "COLON_TOKEN", "HASH_TOKEN", "BOND_TOKEN", "ATOM_OPEN_TOKEN",
  "ATOM_CLOSE_TOKEN", "EOS_TOKEN", "$accept", "meta_start", "mol", "bondd",
  "atomd", "charge_element", "h_element", "chiral_element", "element",
  "simple_atom", "ring_number", "number", "nonzero_number", "digit", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283
};
# endif

#define YYPACT_NINF -32

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-32)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
       5,   -11,     8,     8,   -10,     2,   -32,   -32,   -32,    88,
      54,   -32,   -32,   -32,   -32,   -32,   -32,   -32,    -8,   -32,
     -32,   -32,   -32,     4,    68,    -5,    74,     9,    15,   -32,
      79,   105,   -32,   -32,    67,   -32,     8,    62,    39,    62,
     -32,   -32,   -32,   -32,    68,   -32,    68,   -32,    32,     3,
      68,    18,   -32,    25,    68,   -32,   -32,     8,     8,   -32,
     -32,   -32,   -32,   105,   105,   -32,   -32,   -32,    24,   -32,
     -32,   -32,   -32,   -32,   -32,    68,   -32,   -32,   -32,   -32,
      34,   -32,   -32,   -32,    92,   -32,    97,   -32,   101,   -32,
      49,   -32
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     7,    51,    50,     0,
       2,     8,    22,     3,    21,    20,     4,     1,     0,     6,
      46,    61,    59,    32,     0,     0,    25,    38,    41,    44,
       0,    60,    63,    64,     0,    19,     0,     0,     0,     0,
       9,    13,    52,     5,    34,    48,     0,    24,    29,    26,
      39,    42,    47,    33,     0,    45,    62,     0,     0,    16,
      12,    11,    15,     0,     0,    10,    14,    36,     0,    30,
      31,    27,    28,    40,    43,    35,    49,    18,    17,    53,
       0,    23,    37,    54,     0,    55,     0,    56,     0,    57,
       0,    58
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -32,   -32,   -32,   -32,     1,   -32,   -32,   -32,   -32,    -2,
      17,   -23,   -32,   -31
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     5,    10,    16,    11,    25,    26,    27,    28,    12,
      41,    30,    31,    42
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      56,    45,    17,    18,    13,    14,     1,    29,     2,     3,
       4,    40,    21,    22,     7,    15,     8,     6,    46,    71,
      43,    67,    47,    68,    44,    70,    72,    73,    55,    50,
      19,    76,    79,    80,     9,    59,    51,    60,    61,    74,
      65,    21,    22,    32,    33,    75,    83,    69,    63,    84,
      64,    81,    82,    86,    62,    88,    66,    90,    77,    78,
       7,    91,     8,    32,    33,    34,    35,    36,     7,    37,
       8,    32,    33,     7,     0,     8,    38,    21,    22,    39,
       9,     0,    57,     0,    38,     7,    52,     8,     9,    48,
      49,     0,    58,     9,     7,    20,     8,    21,    22,    53,
       0,    32,    33,    54,    85,     0,    32,    33,    23,    87,
      32,    33,    24,    89,    32,    33
};

static const yytype_int8 yycheck[] =
{
      31,    24,     0,     1,     3,    15,     1,     9,     3,     4,
       5,    10,     9,    10,     6,    25,     8,    28,    23,    16,
      28,    44,    27,    46,    20,    48,    49,    50,    30,    20,
      28,    54,    63,    64,    26,    34,    21,    36,    37,    21,
      39,     9,    10,     9,    10,    20,    12,    15,     9,    80,
      11,    27,    75,    84,    37,    86,    39,    88,    57,    58,
       6,    12,     8,     9,    10,    11,    12,    13,     6,    15,
       8,     9,    10,     6,    -1,     8,    22,     9,    10,    25,
      26,    -1,    15,    -1,    22,     6,     7,     8,    26,    15,
      16,    -1,    25,    26,     6,     7,     8,     9,    10,    20,
      -1,     9,    10,    24,    12,    -1,     9,    10,    20,    12,
       9,    10,    24,    12,     9,    10
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     3,     4,     5,    30,    28,     6,     8,    26,
      31,    33,    38,    33,    15,    25,    32,     0,     1,    28,
       7,     9,    10,    20,    24,    34,    35,    36,    37,    38,
      40,    41,     9,    10,    11,    12,    13,    15,    22,    25,
      33,    39,    42,    28,    20,    40,    23,    27,    15,    16,
      20,    21,     7,    20,    24,    38,    42,    15,    25,    33,
      33,    33,    39,     9,    11,    33,    39,    40,    40,    15,
      40,    16,    40,    40,    21,    20,    40,    33,    33,    42,
      42,    27,    40,    12,    42,    12,    42,    12,    42,    12,
      42,    12
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    29,    30,    30,    30,    30,    30,    30,    31,    31,
      31,    31,    31,    31,    31,    31,    31,    31,    31,    31,
      32,    32,    33,    33,    33,    34,    34,    34,    34,    34,
      34,    34,    35,    35,    35,    35,    35,    35,    35,    35,
      35,    36,    36,    36,    37,    37,    37,    37,    37,    37,
      38,    38,    39,    39,    39,    39,    39,    39,    39,    40,
      40,    41,    41,    42,    42
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     2,     2,     3,     2,     2,     1,     2,
       3,     3,     3,     2,     3,     3,     3,     4,     4,     2,
       1,     1,     1,     5,     3,     1,     2,     3,     3,     2,
       3,     3,     1,     2,     2,     3,     3,     4,     1,     2,
       3,     1,     2,     3,     1,     2,     1,     2,     2,     3,
       1,     1,     1,     3,     4,     5,     6,     7,     8,     1,
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
      yyerror (input, molList, lastAtom, lastBond, branchPoints, scanner, start_token, YY_("syntax error: cannot back up")); \
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
                  Type, Value, input, molList, lastAtom, lastBond, branchPoints, scanner, start_token); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (branchPoints);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, input, molList, lastAtom, lastBond, branchPoints, scanner, start_token);
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
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
                                              , input, molList, lastAtom, lastBond, branchPoints, scanner, start_token);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, molList, lastAtom, lastBond, branchPoints, scanner, start_token); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (branchPoints);
  YYUSE (scanner);
  YYUSE (start_token);
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
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
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
#line 103 "smiles.yy" /* yacc.c:1646  */
    {
// the molList has already been updated, no need to do anything
}
#line 1399 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 3:
#line 106 "smiles.yy" /* yacc.c:1646  */
    {
  lastAtom = (yyvsp[0].atom);
}
#line 1407 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 4:
#line 109 "smiles.yy" /* yacc.c:1646  */
    {
  lastBond = (yyvsp[0].bond);
}
#line 1415 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 5:
#line 112 "smiles.yy" /* yacc.c:1646  */
    {
  yyclearin;
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1426 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 118 "smiles.yy" /* yacc.c:1646  */
    {
  YYACCEPT;
}
#line 1434 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 121 "smiles.yy" /* yacc.c:1646  */
    {
  yyclearin;
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1445 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 8:
#line 132 "smiles.yy" /* yacc.c:1646  */
    {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  RDKit::RWMol *curMol = (*molList)[ sz ];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1);
  curMol->addAtom((yyvsp[0].atom),true,true);
  //delete (yyvsp[0].atom);
  (yyval.moli) = sz;
}
#line 1460 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 9:
#line 143 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  //delete $2;
}
#line 1474 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 10:
#line 153 "smiles.yy" /* yacc.c:1646  */
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
  mp->addBond((yyvsp[-1].bond),true);
  //delete $3;
}
#line 1498 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 11:
#line 173 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  //delete $3;
}
#line 1510 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 12:
#line 181 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1520 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 13:
#line 187 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  mp->setAtomBookmark(atom,(yyvsp[0].ival));

  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::UNSPECIFIED);
  mp->setBondBookmark(newB,(yyvsp[0].ival));
  newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);

}
#line 1543 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 14:
#line 206 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     (yyvsp[-1].bond)->getBondType());
  newB->setBondDir((yyvsp[-1].bond)->getBondDir());
  mp->setAtomBookmark(atom,(yyvsp[0].ival));
  mp->setBondBookmark(newB,(yyvsp[0].ival));

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
  delete (yyvsp[-1].bond);
}
#line 1565 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 15:
#line 224 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::SINGLE);
  mp->setAtomBookmark(atom,(yyvsp[0].ival));
  mp->setBondBookmark(newB,(yyvsp[0].ival));

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
}
#line 1585 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 16:
#line 240 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  //delete $3;
  branchPoints->push_back(atomIdx1);
}
#line 1600 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 17:
#line 250 "smiles.yy" /* yacc.c:1646  */
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
  mp->addBond((yyvsp[-1].bond),true);
  //delete $4;
  branchPoints->push_back(atomIdx1);
}
#line 1625 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 18:
#line 270 "smiles.yy" /* yacc.c:1646  */
    {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  //delete $4;
  branchPoints->push_back(atomIdx1);
}
#line 1638 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 19:
#line 278 "smiles.yy" /* yacc.c:1646  */
    {
  if(branchPoints->empty()) yyerror(input,molList,branchPoints,scanner,start_token,"extra close parentheses");
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->setActiveAtom(branchPoints->back());
  branchPoints->pop_back();
}
#line 1649 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 21:
#line 288 "smiles.yy" /* yacc.c:1646  */
    {
          (yyval.bond) = new Bond(Bond::SINGLE);
          }
#line 1657 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 23:
#line 298 "smiles.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setNoImplicit(true);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1667 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 24:
#line 305 "smiles.yy" /* yacc.c:1646  */
    {
  (yyval.atom) = (yyvsp[-1].atom);
  (yyvsp[-1].atom)->setNoImplicit(true);
}
#line 1676 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 26:
#line 313 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-1].atom)->setFormalCharge(1); }
#line 1682 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 27:
#line 314 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-2].atom)->setFormalCharge(2); }
#line 1688 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 28:
#line 315 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-2].atom)->setFormalCharge((yyvsp[0].ival)); }
#line 1694 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 29:
#line 316 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-1].atom)->setFormalCharge(-1); }
#line 1700 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 30:
#line 317 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-2].atom)->setFormalCharge(-2); }
#line 1706 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 31:
#line 318 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-2].atom)->setFormalCharge(-(yyvsp[0].ival)); }
#line 1712 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 32:
#line 322 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom(1); }
#line 1718 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 33:
#line 323 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-1].ival)); }
#line 1724 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 34:
#line 324 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs(1); }
#line 1730 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 35:
#line 325 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-2].ival)); (yyval.atom)->setNumExplicitHs(1);}
#line 1736 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 36:
#line 326 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs((yyvsp[0].ival)); }
#line 1742 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 37:
#line 327 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-3].ival)); (yyval.atom)->setNumExplicitHs((yyvsp[0].ival));}
#line 1748 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 39:
#line 329 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = (yyvsp[-1].atom); (yyvsp[-1].atom)->setNumExplicitHs(1);}
#line 1754 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 330 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = (yyvsp[-2].atom); (yyvsp[-2].atom)->setNumExplicitHs((yyvsp[0].ival));}
#line 1760 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 42:
#line 335 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-1].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
#line 1766 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 43:
#line 336 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[-2].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
#line 1772 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 45:
#line 341 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[0].atom)->setIsotope( (yyvsp[-1].ival) ); (yyval.atom) = (yyvsp[0].atom); }
#line 1778 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 47:
#line 343 "smiles.yy" /* yacc.c:1646  */
    { (yyvsp[0].atom)->setIsotope( (yyvsp[-1].ival) ); (yyval.atom) = (yyvsp[0].atom); }
#line 1784 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 48:
#line 344 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom((yyvsp[0].ival)); }
#line 1790 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 49:
#line 345 "smiles.yy" /* yacc.c:1646  */
    { (yyval.atom) = new Atom((yyvsp[0].ival)); (yyval.atom)->setIsotope((yyvsp[-2].ival)); }
#line 1796 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 53:
#line 355 "smiles.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 1802 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 54:
#line 356 "smiles.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 1808 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 55:
#line 357 "smiles.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 1814 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 56:
#line 358 "smiles.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 1820 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 57:
#line 359 "smiles.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 1826 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 58:
#line 360 "smiles.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 1832 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;

  case 62:
#line 370 "smiles.yy" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
#line 1838 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
    break;


#line 1842 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp" /* yacc.c:1646  */
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
      yyerror (input, molList, lastAtom, lastBond, branchPoints, scanner, start_token, YY_("syntax error"));
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
        yyerror (input, molList, lastAtom, lastBond, branchPoints, scanner, start_token, yymsgp);
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
                      yytoken, &yylval, input, molList, lastAtom, lastBond, branchPoints, scanner, start_token);
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
                  yystos[yystate], yyvsp, input, molList, lastAtom, lastBond, branchPoints, scanner, start_token);
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
  yyerror (input, molList, lastAtom, lastBond, branchPoints, scanner, start_token, YY_("memory exhausted"));
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
                  yytoken, &yylval, input, molList, lastAtom, lastBond, branchPoints, scanner, start_token);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, input, molList, lastAtom, lastBond, branchPoints, scanner, start_token);
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
#line 384 "smiles.yy" /* yacc.c:1906  */

