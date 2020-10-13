/* A Bison parser, made by GNU Bison 3.5.1.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2020 Free Software Foundation,
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

/* Undocumented macros, especially those whose name start with YY_,
   are private implementation details.  Do not rely on them.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.5.1"

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

/* First part of user prologue.  */
#line 1 "smarts.yy"


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
     SmilesParseOps::CleanupAfterParseError(*iter);
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

#line 127 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
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
#line 68 "smarts.yy"

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  int                      ival;

#line 232 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token);
/* "%code provides" blocks.  */
#line 63 "smarts.yy"

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)

#line 249 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

#endif /* !YY_YYSMARTS_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMARTS_TAB_HPP_INCLUDED  */



#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))

/* Stored state numbers (used for stacks). */
typedef yytype_uint8 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

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

#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
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
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                            \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
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

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

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
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
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
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
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
#define YYLAST   631

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  46
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  21
/* YYNRULES -- Number of rules.  */
#define YYNRULES  118
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  175

#define YYUNDEFTOK  2
#define YYMAXUTOK   300


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
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
static const yytype_int16 yyrline[] =
{
       0,   112,   112,   115,   119,   122,   125,   129,   133,   136,
     141,   144,   152,   153,   154,   155,   163,   171,   196,   216,
     221,   252,   272,   285,   286,   297,   298,   299,   303,   326,
     330,   335,   341,   350,   356,   364,   372,   385,   391,   398,
     404,   427,   430,   436,   437,   441,   458,   482,   483,   488,
     489,   494,   495,   500,   501,   502,   503,   504,   505,   506,
     509,   512,   515,   518,   521,   524,   530,   536,   543,   551,
     559,   565,   571,   577,   583,   589,   590,   597,   598,   601,
     604,   607,   610,   613,   618,   626,   638,   643,   648,   652,
     656,   660,   663,   664,   671,   672,   678,   684,   690,   695,
     702,   703,   704,   705,   706,   707,   711,   712,   713,   714,
     715,   716,   717,   722,   723,   727,   728,   737,   738
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
static const yytype_int16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300
};
# endif

#define YYPACT_NINF (-57)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-85)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     296,   -42,    20,   314,   168,    23,   -57,   -57,   -57,   -57,
     530,   217,   -57,   -57,   -57,   -57,   106,   146,   237,   248,
     -57,   250,   258,   -57,   -57,    57,   239,    99,    57,    -7,
     351,   388,   565,    20,   388,   -57,   -12,   425,   -57,   -57,
     -57,     8,     0,   -57,    31,    76,   -57,   -57,   -57,   168,
     -57,   -57,   143,   168,   -57,   -57,    10,   -57,   585,   277,
     -57,   126,   -57,   -57,   238,    20,   142,   -57,   -57,   580,
     -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,
     -57,   -57,   -57,   -57,   -57,   388,   -57,   277,   -57,   -57,
     134,   -57,   -57,   565,   565,   565,   -57,    67,   -57,    57,
      57,   -57,   -57,   -57,   168,   168,   168,   -57,   -57,   -57,
     201,    64,   -57,    57,     1,   -57,    57,   589,   165,    35,
     -57,    76,    76,   -57,   -57,    36,   565,   495,   460,    57,
      25,   -57,   -57,   -57,    44,   158,    12,   -57,    57,    61,
     -57,    57,    22,   -57,   196,   -57,   316,    85,    89,    41,
     -57,    81,   -57,    95,   -57,    57,   -57,   -57,   334,    76,
     -57,   -57,   112,   -57,   -57,   103,   -57,   353,   -57,   -57,
     -57,   371,   -57,   121,   -57
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     0,     5,     8,     0,    11,    86,    85,    87,
       0,     2,    16,    26,    25,    49,    53,    56,    57,    58,
      75,    54,    55,   113,   115,     0,   105,   102,    71,    74,
       0,     0,     0,     0,     0,     4,     0,    15,    41,    43,
      44,     0,    47,    72,    76,   114,    96,    95,    98,     0,
      97,    94,     7,    91,    92,     1,     0,    10,    71,     0,
      47,    76,   118,   117,     0,     0,     0,    22,    17,     0,
      20,   106,    59,    62,    63,    64,    60,    61,    51,   103,
     104,   100,   101,    70,    73,     0,    12,    15,    13,    42,
       0,    14,     3,     0,     0,     0,    40,     0,    50,     0,
      68,    48,   116,    99,     0,     0,     0,     6,    93,     9,
     105,   102,    29,     0,     0,    27,     0,    68,     0,     0,
      19,     0,     0,    18,    21,    45,    37,    38,    39,     0,
       0,    52,    69,    88,    89,    90,     0,    33,     0,     0,
      31,     0,     0,    23,     0,   107,     0,     0,     0,     0,
      30,     0,    28,     0,    35,     0,    24,   108,     0,    46,
      65,    66,     0,    34,    32,     0,   109,     0,    67,    36,
     110,     0,   111,     0,   112
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -57,   -57,    18,   -14,   -57,    19,   -57,     4,    21,   -57,
     -57,   -57,     2,     6,   -57,   -21,   -56,    75,   -10,     3,
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
      61,   102,   114,     6,    14,    42,    72,    73,    74,    75,
      52,    76,    77,    14,    59,    78,    80,    82,    83,    90,
      61,    35,    36,    55,    56,    84,     7,     8,   103,     9,
      68,    97,   108,    92,    87,    14,   137,     7,     8,    98,
       9,     7,     8,   138,     9,   -25,   101,   150,    83,    88,
     118,   149,    91,    89,    10,   109,    99,   154,    96,    23,
      24,   142,   100,   101,   155,   161,    14,    14,    57,    10,
     119,    14,   104,   105,   106,    23,    24,   145,   146,   147,
      96,   104,    23,    24,   120,    23,    24,   130,   123,   131,
     132,    81,    14,   129,    62,    63,   152,   126,   127,   128,
      80,    82,   158,   136,    24,   144,   139,   132,    96,    68,
     133,   134,   135,   160,   167,   102,   163,    23,    24,   148,
      14,    14,   -83,   171,    23,    24,    81,   173,   151,   -77,
     164,   153,     7,     8,    98,     9,   168,    68,   169,   162,
       7,     8,   174,     9,   124,   165,    14,    96,    96,    96,
     159,    99,    62,    63,    64,     0,    65,   117,     0,    46,
      47,   121,   122,    68,    23,    24,    48,    66,    10,   -80,
      49,     7,     8,     0,     9,   125,    50,     0,    51,     0,
     104,   105,   106,    62,    63,    64,   143,    65,   107,     0,
      46,    47,     0,    46,    47,   104,   105,    48,    66,    10,
      48,    49,     7,     8,    49,     9,     0,    50,     0,    51,
      50,     0,    51,     0,    62,    63,    64,   156,    65,    23,
      24,    46,    47,     7,     8,     0,     9,    79,    48,    66,
      10,     0,    49,     0,     0,    62,    63,    64,    50,    65,
      51,     0,    46,    47,     7,     8,     0,     9,     0,    48,
      66,    10,     0,    49,     0,    23,    24,    23,    24,    50,
     -81,    51,   -84,    46,    47,    79,    23,    24,    23,    24,
      48,   -82,    10,   -78,    49,     0,    23,    24,     0,     0,
      50,   -79,    51,     7,     8,    15,     9,    16,    17,    18,
      19,    20,    21,     0,    22,    23,    24,     1,     0,     2,
       3,     4,    25,    26,    27,     0,     0,     0,    28,    29,
       0,     0,   115,    32,    93,    94,    95,    33,     0,   116,
       7,     8,    15,     9,    16,    17,    18,    19,    20,    21,
       0,    22,    23,    24,    62,    63,     0,   157,     0,    25,
      26,    27,     0,     0,     0,    28,    29,     0,    30,    31,
      32,     0,    62,    63,    33,   166,    34,     7,     8,    15,
       9,    16,    17,    18,    19,    20,    21,     0,    22,    23,
      24,    62,    63,     0,   170,     0,    25,    26,    27,     0,
       0,     0,    58,    29,     0,    85,    31,    32,     0,    62,
      63,    33,   172,    34,     7,     8,    15,     9,    16,    17,
      18,    19,    20,    21,     0,    22,    23,    24,     0,     0,
       0,     0,     0,    25,    26,    27,     0,     0,     0,    28,
      29,     0,    85,    31,    32,     0,     0,     0,    33,     0,
      34,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,     0,    22,    23,    24,     0,     0,     0,     0,     0,
      25,    26,    27,     0,     0,     0,    28,    29,     0,     0,
       0,    32,    93,    94,    95,    33,     7,     8,    15,     9,
      16,    17,    18,    19,    20,    21,     0,    22,    23,    24,
       0,     0,     0,     0,     0,    25,    26,    27,     0,     0,
       0,    28,    29,     0,     0,     0,    32,    93,    94,     0,
      33,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,     0,    22,    23,    24,     0,     0,     0,     0,     0,
      25,    26,    27,     0,     0,     0,    28,    29,     0,     0,
       0,    32,    93,     0,     0,    33,     7,     8,    15,     9,
      16,    17,    18,    19,    20,    21,     0,    22,    23,    24,
       0,     0,     0,     0,     0,    25,    26,    27,     0,     0,
       0,    58,    29,     0,     0,     0,    32,     0,     0,     0,
      33,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,     0,    22,    23,    24,     0,     7,     8,     0,     9,
      25,    26,    27,     0,     0,     0,    28,    29,    62,    63,
       0,    32,     0,    23,    24,    33,     0,    23,    24,     0,
       0,   110,   111,    66,    10,   110,   111,   104,   105,   106,
     112,     0,     0,     0,   140,     0,     0,   113,     0,     0,
       0,   141
};

static const yytype_int16 yycheck[] =
{
      10,    45,    58,    45,     2,     3,    16,    17,    18,    19,
       4,    21,    22,    11,    10,    25,    26,    27,    28,    33,
      30,     3,     3,     0,     1,    32,     6,     7,    49,     9,
      11,    23,    53,    45,    30,    33,    35,     6,     7,     8,
       9,     6,     7,    42,     9,    45,    44,    35,    58,    31,
      64,    26,    34,    32,    34,    45,    25,    35,    37,    18,
      19,   117,    31,    61,    42,    24,    64,    65,    45,    34,
      64,    69,    37,    38,    39,    18,    19,   121,   122,    43,
      59,    37,    18,    19,    65,    18,    19,    97,    69,    99,
     100,    27,    90,    26,    18,    19,    35,    93,    94,    95,
     110,   111,   146,   113,    19,   119,   116,   117,    87,    90,
     104,   105,   106,    24,   158,   159,    35,    18,    19,   129,
     118,   119,    23,   167,    18,    19,    27,   171,   138,    23,
      35,   141,     6,     7,     8,     9,    24,   118,    35,   149,
       6,     7,    21,     9,    69,   155,   144,   126,   127,   128,
     147,    25,    18,    19,    20,    -1,    22,    31,    -1,    25,
      26,    19,    20,   144,    18,    19,    32,    33,    34,    23,
      36,     6,     7,    -1,     9,    41,    42,    -1,    44,    -1,
      37,    38,    39,    18,    19,    20,    21,    22,    45,    -1,
      25,    26,    -1,    25,    26,    37,    38,    32,    33,    34,
      32,    36,     6,     7,    36,     9,    -1,    42,    -1,    44,
      42,    -1,    44,    -1,    18,    19,    20,    21,    22,    18,
      19,    25,    26,     6,     7,    -1,     9,    26,    32,    33,
      34,    -1,    36,    -1,    -1,    18,    19,    20,    42,    22,
      44,    -1,    25,    26,     6,     7,    -1,     9,    -1,    32,
      33,    34,    -1,    36,    -1,    18,    19,    18,    19,    42,
      23,    44,    23,    25,    26,    26,    18,    19,    18,    19,
      32,    23,    34,    23,    36,    -1,    18,    19,    -1,    -1,
      42,    23,    44,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    -1,    17,    18,    19,     1,    -1,     3,
       4,     5,    25,    26,    27,    -1,    -1,    -1,    31,    32,
      -1,    -1,    35,    36,    37,    38,    39,    40,    -1,    42,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      -1,    17,    18,    19,    18,    19,    -1,    21,    -1,    25,
      26,    27,    -1,    -1,    -1,    31,    32,    -1,    34,    35,
      36,    -1,    18,    19,    40,    21,    42,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    -1,    17,    18,
      19,    18,    19,    -1,    21,    -1,    25,    26,    27,    -1,
      -1,    -1,    31,    32,    -1,    34,    35,    36,    -1,    18,
      19,    40,    21,    42,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    -1,    17,    18,    19,    -1,    -1,
      -1,    -1,    -1,    25,    26,    27,    -1,    -1,    -1,    31,
      32,    -1,    34,    35,    36,    -1,    -1,    -1,    40,    -1,
      42,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    -1,    17,    18,    19,    -1,    -1,    -1,    -1,    -1,
      25,    26,    27,    -1,    -1,    -1,    31,    32,    -1,    -1,
      -1,    36,    37,    38,    39,    40,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    -1,    17,    18,    19,
      -1,    -1,    -1,    -1,    -1,    25,    26,    27,    -1,    -1,
      -1,    31,    32,    -1,    -1,    -1,    36,    37,    38,    -1,
      40,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    -1,    17,    18,    19,    -1,    -1,    -1,    -1,    -1,
      25,    26,    27,    -1,    -1,    -1,    31,    32,    -1,    -1,
      -1,    36,    37,    -1,    -1,    40,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    -1,    17,    18,    19,
      -1,    -1,    -1,    -1,    -1,    25,    26,    27,    -1,    -1,
      -1,    31,    32,    -1,    -1,    -1,    36,    -1,    -1,    -1,
      40,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    -1,    17,    18,    19,    -1,     6,     7,    -1,     9,
      25,    26,    27,    -1,    -1,    -1,    31,    32,    18,    19,
      -1,    36,    -1,    18,    19,    40,    -1,    18,    19,    -1,
      -1,    26,    27,    33,    34,    26,    27,    37,    38,    39,
      35,    -1,    -1,    -1,    35,    -1,    -1,    42,    -1,    -1,
      -1,    42
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_int8 yystos[] =
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
      26,    27,    35,    42,    62,    35,    42,    31,    49,    59,
      51,    19,    20,    51,    63,    41,    53,    53,    53,    26,
      64,    64,    64,    59,    59,    59,    64,    35,    42,    64,
      35,    42,    62,    21,    49,    66,    66,    43,    64,    26,
      35,    64,    35,    64,    35,    42,    21,    21,    66,    65,
      24,    24,    64,    35,    35,    64,    21,    66,    24,    35,
      21,    66,    21,    66,    21
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_int8 yyr1[] =
{
       0,    46,    47,    47,    47,    47,    47,    47,    47,    47,
      47,    47,    48,    48,    48,    48,    49,    49,    49,    49,
      49,    49,    49,    50,    50,    51,    51,    51,    51,    52,
      52,    52,    52,    52,    52,    52,    52,    53,    53,    53,
      53,    53,    54,    54,    54,    55,    55,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    57,    57,    57,
      57,    57,    57,    57,    57,    58,    58,    58,    59,    59,
      59,    59,    60,    60,    61,    61,    61,    61,    61,    61,
      62,    62,    62,    62,    62,    62,    63,    63,    63,    63,
      63,    63,    63,    64,    64,    65,    65,    66,    66
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     3,     2,     1,     3,     2,     1,     3,
       2,     2,     2,     2,     2,     1,     1,     2,     3,     3,
       2,     3,     2,     3,     4,     1,     1,     3,     5,     3,
       5,     4,     6,     4,     6,     5,     7,     3,     3,     3,
       2,     1,     2,     1,     1,     3,     5,     1,     2,     1,
       2,     2,     3,     1,     1,     1,     1,     1,     1,     2,
       2,     2,     2,     2,     2,     5,     5,     6,     2,     3,
       2,     1,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     3,
       3,     1,     1,     2,     1,     1,     1,     1,     1,     2,
       2,     2,     1,     2,     2,     1,     1,     3,     4,     5,
       6,     7,     8,     1,     1,     1,     2,     1,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
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


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
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
    YYPRINT (yyo, yytoknum[yytype], *yyvaluep);
# endif
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyo, yytype, yyvaluep, input, molList, lastAtom, lastBond, scanner, start_token);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
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
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, void *scanner, int& start_token)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[+yyssp[yyi + 1 - yynrhs]],
                       &yyvsp[(yyi + 1) - (yynrhs)]
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
#   define yystrlen(S) (YY_CAST (YYPTRDIFF_T, strlen (S)))
#  else
/* Return the length of YYSTR.  */
static YYPTRDIFF_T
yystrlen (const char *yystr)
{
  YYPTRDIFF_T yylen;
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
static YYPTRDIFF_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYPTRDIFF_T yyn = 0;
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
            else
              goto append;

          append:
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

  if (yyres)
    return yystpcpy (yyres, yystr) - yyres;
  else
    return yystrlen (yystr);
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
yysyntax_error (YYPTRDIFF_T *yymsg_alloc, char **yymsg,
                yy_state_t *yyssp, int yytoken)
{
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat: reported tokens (one for the "unexpected",
     one per "expected"). */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Actual size of YYARG. */
  int yycount = 0;
  /* Cumulated lengths of YYARG.  */
  YYPTRDIFF_T yysize = 0;

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
      int yyn = yypact[+*yyssp];
      YYPTRDIFF_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
      yysize = yysize0;
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
                  YYPTRDIFF_T yysize1
                    = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
                    yysize = yysize1;
                  else
                    return 2;
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
    /* Don't count the "%s"s in the final size, but reserve room for
       the terminator.  */
    YYPTRDIFF_T yysize1 = yysize + (yystrlen (yyformat) - 2 * yycount) + 1;
    if (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)
      yysize = yysize1;
    else
      return 2;
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
          ++yyp;
          ++yyformat;
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
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1391 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 9: /* SIMPLE_ATOM_QUERY_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1397 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 10: /* COMPLEX_ATOM_QUERY_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1403 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 11: /* RINGSIZE_ATOM_QUERY_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1409 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 12: /* RINGBOND_ATOM_QUERY_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1415 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 13: /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1421 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 14: /* HYB_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1427 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 15: /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1433 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 16: /* ALIPHATIC  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1439 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 17: /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1445 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 44: /* BOND_TOKEN  */
#line 104 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1451 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 51: /* atomd  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1457 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 52: /* hydrogen_atom  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1463 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 53: /* atom_expr  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1469 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 54: /* point_query  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1475 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 55: /* recursive_query  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1481 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 56: /* atom_query  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1487 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 57: /* possible_range_query  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1493 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 58: /* simple_atom  */
#line 103 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1499 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 59: /* bond_expr  */
#line 104 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1505 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 60: /* bond_query  */
#line 104 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1511 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 61: /* bondd  */
#line 104 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1517 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
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

    yy_state_fast_t yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss;
    yy_state_t *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYPTRDIFF_T yystacksize;

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
  YYPTRDIFF_T yymsg_alloc = sizeof yymsgbuf;
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
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
# undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

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
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
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
| yyreduce -- do a reduction.  |
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
#line 112 "smarts.yy"
              {
// the molList has already been updated, no need to do anything
}
#line 1793 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 3:
#line 115 "smarts.yy"
                             {
  lastAtom = (yyvsp[-1].atom);
  YYACCEPT;
}
#line 1802 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 4:
#line 119 "smarts.yy"
                          {
  YYABORT;
}
#line 1810 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 5:
#line 122 "smarts.yy"
             {
  YYABORT;
}
#line 1818 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 6:
#line 125 "smarts.yy"
                                 {
  lastBond = (yyvsp[-1].bond);
  YYACCEPT;
}
#line 1827 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 7:
#line 129 "smarts.yy"
                       {
  delete (yyvsp[0].bond);
  YYABORT;
}
#line 1836 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 8:
#line 133 "smarts.yy"
             {
  YYABORT;
}
#line 1844 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 9:
#line 136 "smarts.yy"
                            {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1854 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 10:
#line 141 "smarts.yy"
                       {
  YYACCEPT;
}
#line 1862 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 11:
#line 144 "smarts.yy"
                  {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1872 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 15:
#line 155 "smarts.yy"
            {
  delete (yyvsp[0].atom);
  YYABORT;
}
#line 1881 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 16:
#line 163 "smarts.yy"
           {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (*molList)[ sz ]->addAtom((yyvsp[0].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
}
#line 1894 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 17:
#line 171 "smarts.yy"
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
#line 1923 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 18:
#line 196 "smarts.yy"
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
#line 1947 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 19:
#line 216 "smarts.yy"
                            {
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1956 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 20:
#line 221 "smarts.yy"
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
#line 1991 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 21:
#line 252 "smarts.yy"
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
#line 2015 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 22:
#line 272 "smarts.yy"
             {
  RWMol *m1_p = (*molList)[(yyval.moli)],*m2_p=(*molList)[(yyvsp[0].moli)];
  // FIX: handle generic bonds here
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED,Bond::NONE);
  delete m2_p;
  int sz = molList->size();
  if ( sz==(yyvsp[0].moli)+1) {
    molList->resize( sz-1 );
  }
}
#line 2030 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 23:
#line 285 "smarts.yy"
                                               { (yyval.moli) = (yyvsp[-1].moli); }
#line 2036 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 24:
#line 286 "smarts.yy"
                                                   {
  // FIX: this needs to handle arbitrary bond_exprs
  (yyval.moli) = (yyvsp[-1].moli);
  int sz     = molList->size();
  (yyvsp[-2].bond)->setOwningMol((*molList)[ sz-1 ]);
  (yyvsp[-2].bond)->setBeginAtomIdx(0);
  (*molList)[ sz-1 ]->setBondBookmark((yyvsp[-2].bond),ci_LEADING_BOND);
}
#line 2049 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 27:
#line 300 "smarts.yy"
{
  (yyval.atom) = (yyvsp[-1].atom);
}
#line 2057 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 28:
#line 304 "smarts.yy"
{
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 2066 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 29:
#line 327 "smarts.yy"
{
  (yyval.atom) = new QueryAtom(1);
}
#line 2074 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 30:
#line 331 "smarts.yy"
{
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 2083 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 31:
#line 335 "smarts.yy"
                                                  {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2094 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 32:
#line 341 "smarts.yy"
                                                                     {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-4].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-4].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2107 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 33:
#line 350 "smarts.yy"
                                                       {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2118 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 34:
#line 356 "smarts.yy"
                                                                          {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2131 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 35:
#line 364 "smarts.yy"
                                                              {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-3].ival));
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2144 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 36:
#line 372 "smarts.yy"
                                                                                 {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-5].ival));
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-5].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2159 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 37:
#line 385 "smarts.yy"
                                         {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 2170 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 38:
#line 391 "smarts.yy"
                               {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  (yyvsp[-2].atom)->setAtomicNum(0);
  delete (yyvsp[0].atom);
}
#line 2182 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 39:
#line 398 "smarts.yy"
                                 {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 2193 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 40:
#line 404 "smarts.yy"
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
#line 2221 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 42:
#line 430 "smarts.yy"
                                   {
  (yyvsp[0].atom)->getQuery()->setNegation(!((yyvsp[0].atom)->getQuery()->getNegation()));
  (yyvsp[0].atom)->setAtomicNum(0);
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[0].atom));
  (yyval.atom) = (yyvsp[0].atom);
}
#line 2232 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 45:
#line 441 "smarts.yy"
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
#line 2254 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 46:
#line 458 "smarts.yy"
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
#line 2280 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 48:
#line 483 "smarts.yy"
                     {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2290 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 50:
#line 489 "smarts.yy"
                    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2300 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 51:
#line 494 "smarts.yy"
                    { (yyval.atom) = new QueryAtom((yyvsp[0].ival)); }
#line 2306 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 52:
#line 495 "smarts.yy"
                           {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsotope((yyvsp[-2].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
}
#line 2316 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 59:
#line 506 "smarts.yy"
                                  {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
}
#line 2324 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 60:
#line 509 "smarts.yy"
                                         {
  (yyvsp[-1].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2332 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 61:
#line 512 "smarts.yy"
                                                  {
  (yyvsp[-1].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2340 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 62:
#line 515 "smarts.yy"
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
}
#line 2348 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 63:
#line 518 "smarts.yy"
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
}
#line 2356 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 64:
#line 521 "smarts.yy"
                                     {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
}
#line 2364 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 65:
#line 524 "smarts.yy"
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>((yyvsp[-1].ival),oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2375 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 66:
#line 530 "smarts.yy"
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>((yyvsp[-2].ival),oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2386 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 67:
#line 536 "smarts.yy"
                                                                                    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
}
#line 2398 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 68:
#line 543 "smarts.yy"
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)));
  newQ->setIsotope((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2411 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 69:
#line 551 "smarts.yy"
                        {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)));
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomHCountQuery((yyvsp[0].ival)),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2424 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 70:
#line 559 "smarts.yy"
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2435 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 71:
#line 565 "smarts.yy"
          {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2446 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 72:
#line 571 "smarts.yy"
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  newQ->setFormalCharge((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2457 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 73:
#line 577 "smarts.yy"
                    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
#line 2468 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 74:
#line 583 "smarts.yy"
           {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
#line 2479 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 76:
#line 590 "smarts.yy"
         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 2489 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 78:
#line 598 "smarts.yy"
                                  {
  (yyvsp[0].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
}
#line 2497 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 79:
#line 601 "smarts.yy"
                                           {
  (yyvsp[0].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
}
#line 2505 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 80:
#line 604 "smarts.yy"
                            {
  (yyvsp[0].atom)->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
}
#line 2513 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 81:
#line 607 "smarts.yy"
                            {
  (yyvsp[0].atom)->setQuery(makeAtomRingBondCountQuery(0));
}
#line 2521 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 82:
#line 610 "smarts.yy"
                              {
  (yyvsp[0].atom)->setQuery(makeAtomImplicitHCountQuery(0));
}
#line 2529 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 83:
#line 613 "smarts.yy"
             {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
#line 2539 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 84:
#line 618 "smarts.yy"
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNegativeFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
#line 2549 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 85:
#line 626 "smarts.yy"
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
#line 2566 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 86:
#line 638 "smarts.yy"
                      {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->setQuery(makeAtomTypeQuery((yyvsp[0].ival),true));
}
#line 2576 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 88:
#line 648 "smarts.yy"
                                        {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2585 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 89:
#line 652 "smarts.yy"
                               {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
}
#line 2594 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 90:
#line 656 "smarts.yy"
                                 {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2603 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 93:
#line 664 "smarts.yy"
                   {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2612 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 95:
#line 672 "smarts.yy"
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
#line 2623 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 96:
#line 678 "smarts.yy"
             {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
#line 2634 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 97:
#line 684 "smarts.yy"
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
#line 2645 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 98:
#line 690 "smarts.yy"
           {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
#line 2655 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 99:
#line 695 "smarts.yy"
                  {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
#line 2664 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 100:
#line 702 "smarts.yy"
                                   { (yyval.ival)=2; }
#line 2670 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 101:
#line 703 "smarts.yy"
                    { (yyval.ival)=(yyvsp[0].ival); }
#line 2676 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 102:
#line 704 "smarts.yy"
             { (yyval.ival)=1; }
#line 2682 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 103:
#line 705 "smarts.yy"
                          { (yyval.ival)=-2; }
#line 2688 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 104:
#line 706 "smarts.yy"
                     { (yyval.ival)=-(yyvsp[0].ival); }
#line 2694 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 105:
#line 707 "smarts.yy"
              { (yyval.ival)=-1; }
#line 2700 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 107:
#line 712 "smarts.yy"
                                          { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 2706 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 108:
#line 713 "smarts.yy"
                                                         { (yyval.ival) = (yyvsp[-1].ival); }
#line 2712 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 109:
#line 714 "smarts.yy"
                                                               { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2718 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 110:
#line 715 "smarts.yy"
                                                                     { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2724 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 111:
#line 716 "smarts.yy"
                                                                           { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2730 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 112:
#line 717 "smarts.yy"
                                                                                 { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2736 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 116:
#line 728 "smarts.yy"
                       { 
    if((yyvsp[-1].ival) >= std::numeric_limits<std::int32_t>::max()/10 || 
     (yyvsp[-1].ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yyvsp[0].ival) ){
     yysmarts_error(input,molList,lastAtom,lastBond,scanner,start_token,"number too large");
     YYABORT;
  }
  (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
#line 2748 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;


#line 2752 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

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
            yymsg = YY_CAST (char *, YYSTACK_ALLOC (YY_CAST (YYSIZE_T, yymsg_alloc)));
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
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

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


/*-----------------------------------------------------.
| yyreturn -- parsing is finished, return the result.  |
`-----------------------------------------------------*/
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
                  yystos[+*yyssp], yyvsp, input, molList, lastAtom, lastBond, scanner, start_token);
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
#line 741 "smarts.yy"

