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


  //
  //  Copyright (C) 2003-2022 Greg Landrum and other RDKit contributors
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
                RDKit::Atom* &,
                RDKit::Bond* &,
                unsigned int &,unsigned int &,
		void *,int , const char *msg  )
{
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error: " << msg << " while parsing: " << input << std::endl;
}

#line 122 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
    CHI_CLASS_TOKEN = 297,
    EOS_TOKEN = 298
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 65 "smarts.yy"

  int                      moli;
  RDKit::QueryAtom * atom;
  RDKit::QueryBond * bond;
  RDKit::Atom::ChiralType chiraltype;
  int                      ival;

#line 226 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmarts_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, void *scanner, int& start_token);
/* "%code provides" blocks.  */
#line 60 "smarts.yy"

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)

#line 243 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
#define YYFINAL  56
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   647

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  44
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  21
/* YYNRULES -- Number of rules.  */
#define YYNRULES  120
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  177

#define YYUNDEFTOK  2
#define YYMAXUTOK   298


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
      35,    36,    37,    38,    39,    40,    41,    42,    43
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   110,   110,   113,   117,   120,   123,   127,   131,   134,
     139,   142,   150,   151,   152,   153,   161,   169,   195,   216,
     221,   255,   276,   295,   296,   307,   308,   309,   313,   336,
     340,   345,   351,   360,   366,   374,   382,   395,   401,   408,
     414,   437,   440,   446,   447,   451,   468,   492,   493,   498,
     499,   504,   505,   510,   511,   512,   513,   514,   515,   516,
     519,   522,   525,   528,   531,   534,   540,   546,   553,   561,
     569,   575,   581,   587,   593,   599,   606,   613,   614,   621,
     622,   625,   628,   631,   634,   637,   642,   650,   662,   667,
     672,   676,   680,   684,   687,   688,   695,   696,   702,   708,
     714,   719,   726,   727,   728,   729,   730,   731,   735,   736,
     737,   738,   739,   740,   741,   746,   747,   751,   752,   761,
     762
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
  "MINUS_TOKEN", "PLUS_TOKEN", "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN",
  "ATOM_OPEN_TOKEN", "ATOM_CLOSE_TOKEN", "NOT_TOKEN", "AND_TOKEN",
  "OR_TOKEN", "SEMI_TOKEN", "BEGIN_RECURSE", "END_RECURSE", "COLON_TOKEN",
  "UNDERSCORE_TOKEN", "BOND_TOKEN", "CHI_CLASS_TOKEN", "EOS_TOKEN",
  "$accept", "meta_start", "bad_atom_def", "mol", "branch", "atomd",
  "hydrogen_atom", "atom_expr", "point_query", "recursive_query",
  "atom_query", "possible_range_query", "simple_atom", "bond_expr",
  "bond_query", "bondd", "charge_spec", "ring_number", "number",
  "nonzero_number", "digit", YY_NULLPTR
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
     295,   296,   297,   298
};
# endif

#define YYPACT_NINF (-57)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-87)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     164,   -39,    28,   182,   573,    27,   -57,   -57,   -57,   -57,
     404,   550,   -57,   -57,   -57,   -57,   199,   236,   273,   343,
     -57,   379,   416,   -57,   -57,    69,    22,    20,    69,    13,
     219,   256,   441,    28,   256,    69,   -57,   -22,   293,   -57,
     -57,   -57,    29,    17,   -57,   325,    96,   -57,   -57,   -57,
     573,   -57,   -57,    61,   573,   -57,   -57,    37,   -57,   606,
     145,   -57,   362,   -57,   -57,   578,    28,   123,   -57,   -57,
      48,   -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,
     -57,   -57,   -57,   -57,   -57,   -57,   256,   -57,   145,   -57,
     -57,   466,   -57,   -57,   -57,   441,   441,   441,   -57,    55,
     -57,    69,    69,   -57,   -57,   -57,   573,   573,   573,   -57,
     -57,   -57,    80,    58,   -57,    69,    -8,   -57,    69,   608,
     494,   587,   -57,    96,    96,   -57,   -57,    31,   441,   367,
     330,    69,    30,   -57,   -57,   -57,    60,   141,    78,   -57,
      69,    91,   -57,    69,    92,   -57,   522,   -57,   184,    43,
     105,   116,   -57,   100,   -57,   104,   -57,    69,   -57,   -57,
     221,    96,   -57,   -57,   115,   -57,   -57,   112,   -57,   258,
     -57,   -57,   -57,   295,   -57,   128,   -57
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     0,     5,     8,     0,    11,    88,    87,    89,
       0,     2,    16,    26,    25,    49,    53,    56,    57,    58,
      77,    54,    55,   115,   117,     0,   107,   104,    71,    74,
       0,     0,     0,     0,     0,    75,     4,     0,    15,    41,
      43,    44,     0,    47,    72,    78,   116,    98,    97,   100,
       0,    99,    96,     7,    93,    94,     1,     0,    10,    71,
       0,    47,    78,   120,   119,     0,     0,     0,    22,    17,
       0,    20,   108,    59,    62,    63,    64,    60,    61,    51,
     105,   106,   102,   103,    70,    73,     0,    12,    15,    13,
      42,     0,    14,    76,     3,     0,     0,     0,    40,     0,
      50,     0,    68,    48,   118,   101,     0,     0,     0,     6,
      95,     9,   107,   104,    29,     0,     0,    27,     0,    68,
       0,     0,    19,     0,     0,    18,    21,    45,    37,    38,
      39,     0,     0,    52,    69,    90,    91,    92,     0,    33,
       0,     0,    31,     0,     0,    23,     0,   109,     0,     0,
       0,     0,    30,     0,    28,     0,    35,     0,    24,   110,
       0,    46,    65,    66,     0,    34,    32,     0,   111,     0,
      67,    36,   112,     0,   113,     0,   114
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -57,   -57,    19,   -14,   -57,     2,   -57,    16,    -2,   -57,
     -57,   -57,    -1,    10,   -57,   -21,   -56,   113,   -10,    12,
     -23
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     5,    87,    11,    68,    12,    13,    38,    39,    40,
      41,    42,    61,    70,    54,    55,    44,    71,    45,    46,
      72
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      62,    14,    43,   116,     6,    37,    73,    74,    75,    76,
      14,    77,    78,    69,    53,    79,    81,    83,    84,    91,
      62,    94,    36,   104,   139,    93,    60,    56,    57,   105,
      90,   140,    14,   110,     7,     8,    98,     9,    23,    24,
      23,    24,    85,   -85,   103,   -86,    88,    82,    80,    84,
      89,   120,    99,    92,     7,     8,   151,     9,    98,    10,
     -25,   103,    24,   144,    14,    14,    63,    64,   122,    14,
      58,   149,   125,    23,    24,   121,    23,    24,    67,    10,
     111,   131,   106,   107,   108,    82,    98,    23,    24,   132,
      14,   133,   134,    69,   106,   106,   107,   108,    23,    24,
     147,   148,    81,    83,   109,   138,    80,   146,   141,   134,
     152,   128,   129,   130,    63,    64,   135,   136,   137,    14,
      14,   150,    69,   154,   156,   160,    98,    98,    98,   162,
     153,   157,   165,   155,    23,    24,   166,   169,   104,   170,
     163,   164,   123,   124,   171,    14,   173,   167,    69,   176,
     175,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,   161,    22,    23,    24,     1,     0,     2,     3,     4,
      25,    26,    27,    28,    29,   106,   107,   117,    32,    95,
      96,    97,    33,   126,   118,     0,     0,    35,     7,     8,
      15,     9,    16,    17,    18,    19,    20,    21,     0,    22,
      23,    24,    63,    64,     0,   159,     0,    25,    26,    27,
      28,    29,     0,    30,    31,    32,     0,    23,    24,    33,
       0,    34,   -79,     0,    35,     7,     8,    15,     9,    16,
      17,    18,    19,    20,    21,     0,    22,    23,    24,    63,
      64,     0,   168,     0,    25,    26,    27,    59,    29,     0,
      86,    31,    32,     0,    23,    24,    33,     0,    34,   -82,
       0,    35,     7,     8,    15,     9,    16,    17,    18,    19,
      20,    21,     0,    22,    23,    24,    63,    64,     0,   172,
       0,    25,    26,    27,    28,    29,     0,    86,    31,    32,
       0,    23,    24,    33,     0,    34,   -83,     0,    35,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,     0,
      22,    23,    24,    63,    64,     0,   174,     0,    25,    26,
      27,    28,    29,     0,     0,     0,    32,    95,    96,    97,
      33,     7,     8,   100,     9,    35,     7,     8,    15,     9,
      16,    17,    18,    19,    20,    21,     0,    22,    23,    24,
     101,     0,     0,   102,     0,    25,    26,    27,    28,    29,
       0,    23,    24,    32,    95,    96,   -84,    33,     7,     8,
     100,     9,    35,     7,     8,    15,     9,    16,    17,    18,
      19,    20,    21,     0,    22,    23,    24,   101,     0,     0,
     119,     0,    25,    26,    27,    28,    29,    23,    24,     0,
      32,    95,   -80,     0,    33,     0,     0,     0,     0,    35,
       7,     8,    15,     9,    16,    17,    18,    19,    20,    21,
       0,    22,    23,    24,     0,     0,     0,     0,     0,    25,
      26,    27,    59,    29,    23,    24,     0,    32,     0,   -81,
       0,    33,     0,     0,     0,     0,    35,     7,     8,    15,
       9,    16,    17,    18,    19,    20,    21,     0,    22,    23,
      24,     0,     0,     0,     0,     0,    25,    26,    27,    28,
      29,     0,     7,     8,    32,     9,     0,     0,    33,     0,
       0,     0,     0,    35,    63,    64,    65,     0,    66,     0,
       0,    47,    48,     0,     0,    49,    67,    10,     0,    50,
       7,     8,     0,     9,   127,    51,     0,    52,     0,     0,
       0,     0,    63,    64,    65,   145,    66,     0,     0,    47,
      48,     0,     0,    49,    67,    10,     0,    50,     7,     8,
       0,     9,     0,    51,     0,    52,     0,     0,     0,     0,
      63,    64,    65,   158,    66,     0,     0,    47,    48,     0,
       0,    49,    67,    10,     0,    50,     7,     8,     0,     9,
       0,    51,     0,    52,     0,     0,     0,     0,    63,    64,
      65,     0,    66,     0,     0,    47,    48,     0,     0,    49,
      67,    10,     0,    50,     7,     8,     0,     9,     0,    51,
       0,    52,     0,     7,     8,     0,     9,     0,    47,    48,
       0,     0,    49,    47,    48,     0,    50,    49,     0,    10,
       0,    50,    51,     0,    52,     0,     0,    51,    10,    52,
       0,   106,   107,   108,    23,    24,    23,    24,     0,     0,
       0,     0,   112,   113,   112,   113,     0,     0,   114,     0,
     142,     0,     0,     0,     0,   115,     0,   143
};

static const yytype_int16 yycheck[] =
{
      10,     2,     3,    59,    43,     3,    16,    17,    18,    19,
      11,    21,    22,    11,     4,    25,    26,    27,    28,    33,
      30,    43,     3,    46,    32,    35,    10,     0,     1,    50,
      32,    39,    33,    54,     6,     7,    38,     9,    18,    19,
      18,    19,    29,    23,    45,    23,    30,    27,    26,    59,
      31,    65,    23,    34,     6,     7,    26,     9,    60,    31,
      43,    62,    19,   119,    65,    66,    18,    19,    66,    70,
      43,    40,    70,    18,    19,    65,    18,    19,    30,    31,
      43,    26,    34,    35,    36,    27,    88,    18,    19,    99,
      91,   101,   102,    91,    34,    34,    35,    36,    18,    19,
     123,   124,   112,   113,    43,   115,    26,   121,   118,   119,
      32,    95,    96,    97,    18,    19,   106,   107,   108,   120,
     121,   131,   120,    32,    32,   148,   128,   129,   130,    24,
     140,    39,    32,   143,    18,    19,    32,   160,   161,    24,
      24,   151,    19,    20,    32,   146,   169,   157,   146,    21,
     173,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,   149,    17,    18,    19,     1,    -1,     3,     4,     5,
      25,    26,    27,    28,    29,    34,    35,    32,    33,    34,
      35,    36,    37,    70,    39,    -1,    -1,    42,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    -1,    17,
      18,    19,    18,    19,    -1,    21,    -1,    25,    26,    27,
      28,    29,    -1,    31,    32,    33,    -1,    18,    19,    37,
      -1,    39,    23,    -1,    42,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    -1,    17,    18,    19,    18,
      19,    -1,    21,    -1,    25,    26,    27,    28,    29,    -1,
      31,    32,    33,    -1,    18,    19,    37,    -1,    39,    23,
      -1,    42,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    -1,    17,    18,    19,    18,    19,    -1,    21,
      -1,    25,    26,    27,    28,    29,    -1,    31,    32,    33,
      -1,    18,    19,    37,    -1,    39,    23,    -1,    42,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    -1,
      17,    18,    19,    18,    19,    -1,    21,    -1,    25,    26,
      27,    28,    29,    -1,    -1,    -1,    33,    34,    35,    36,
      37,     6,     7,     8,     9,    42,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    -1,    17,    18,    19,
      25,    -1,    -1,    28,    -1,    25,    26,    27,    28,    29,
      -1,    18,    19,    33,    34,    35,    23,    37,     6,     7,
       8,     9,    42,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    -1,    17,    18,    19,    25,    -1,    -1,
      28,    -1,    25,    26,    27,    28,    29,    18,    19,    -1,
      33,    34,    23,    -1,    37,    -1,    -1,    -1,    -1,    42,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      -1,    17,    18,    19,    -1,    -1,    -1,    -1,    -1,    25,
      26,    27,    28,    29,    18,    19,    -1,    33,    -1,    23,
      -1,    37,    -1,    -1,    -1,    -1,    42,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    -1,    17,    18,
      19,    -1,    -1,    -1,    -1,    -1,    25,    26,    27,    28,
      29,    -1,     6,     7,    33,     9,    -1,    -1,    37,    -1,
      -1,    -1,    -1,    42,    18,    19,    20,    -1,    22,    -1,
      -1,    25,    26,    -1,    -1,    29,    30,    31,    -1,    33,
       6,     7,    -1,     9,    38,    39,    -1,    41,    -1,    -1,
      -1,    -1,    18,    19,    20,    21,    22,    -1,    -1,    25,
      26,    -1,    -1,    29,    30,    31,    -1,    33,     6,     7,
      -1,     9,    -1,    39,    -1,    41,    -1,    -1,    -1,    -1,
      18,    19,    20,    21,    22,    -1,    -1,    25,    26,    -1,
      -1,    29,    30,    31,    -1,    33,     6,     7,    -1,     9,
      -1,    39,    -1,    41,    -1,    -1,    -1,    -1,    18,    19,
      20,    -1,    22,    -1,    -1,    25,    26,    -1,    -1,    29,
      30,    31,    -1,    33,     6,     7,    -1,     9,    -1,    39,
      -1,    41,    -1,     6,     7,    -1,     9,    -1,    25,    26,
      -1,    -1,    29,    25,    26,    -1,    33,    29,    -1,    31,
      -1,    33,    39,    -1,    41,    -1,    -1,    39,    31,    41,
      -1,    34,    35,    36,    18,    19,    18,    19,    -1,    -1,
      -1,    -1,    26,    27,    26,    27,    -1,    -1,    32,    -1,
      32,    -1,    -1,    -1,    -1,    39,    -1,    39
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,    45,    43,     6,     7,     9,
      31,    47,    49,    50,    56,     8,    10,    11,    12,    13,
      14,    15,    17,    18,    19,    25,    26,    27,    28,    29,
      31,    32,    33,    37,    39,    42,    46,    49,    51,    52,
      53,    54,    55,    56,    60,    62,    63,    25,    26,    29,
      33,    39,    41,    57,    58,    59,     0,     1,    43,    28,
      51,    56,    62,    18,    19,    20,    22,    30,    48,    49,
      57,    61,    64,    62,    62,    62,    62,    62,    62,    62,
      26,    62,    27,    62,    62,    29,    31,    46,    51,    46,
      52,    47,    46,    62,    43,    34,    35,    36,    52,    23,
       8,    25,    28,    56,    64,    59,    34,    35,    36,    43,
      59,    43,    26,    27,    32,    39,    60,    32,    39,    28,
      47,    57,    49,    19,    20,    49,    61,    38,    51,    51,
      51,    26,    62,    62,    62,    57,    57,    57,    62,    32,
      39,    62,    32,    39,    60,    21,    47,    64,    64,    40,
      62,    26,    32,    62,    32,    62,    32,    39,    21,    21,
      64,    63,    24,    24,    62,    32,    32,    62,    21,    64,
      24,    32,    21,    64,    21,    64,    21
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_int8 yyr1[] =
{
       0,    44,    45,    45,    45,    45,    45,    45,    45,    45,
      45,    45,    46,    46,    46,    46,    47,    47,    47,    47,
      47,    47,    47,    48,    48,    49,    49,    49,    49,    50,
      50,    50,    50,    50,    50,    50,    50,    51,    51,    51,
      51,    51,    52,    52,    52,    53,    53,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    55,
      55,    55,    55,    55,    55,    55,    55,    56,    56,    56,
      57,    57,    57,    57,    58,    58,    59,    59,    59,    59,
      59,    59,    60,    60,    60,    60,    60,    60,    61,    61,
      61,    61,    61,    61,    61,    62,    62,    63,    63,    64,
      64
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
       2,     1,     1,     2,     1,     1,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       3,     3,     3,     1,     1,     2,     1,     1,     1,     1,
       1,     2,     2,     2,     1,     2,     2,     1,     1,     3,
       4,     5,     6,     7,     8,     1,     1,     1,     2,     1,
       1
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
        yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token, YY_("syntax error: cannot back up")); \
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
                  Type, Value, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, void *scanner, int& start_token)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (numAtomsParsed);
  YYUSE (numBondsParsed);
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
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, void *scanner, int& start_token)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyo, yytype, yyvaluep, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token);
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
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, void *scanner, int& start_token)
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
                                              , input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, void *scanner, int& start_token)
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (numAtomsParsed);
  YYUSE (numBondsParsed);
  YYUSE (scanner);
  YYUSE (start_token);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yytype)
    {
    case 8: /* ATOM_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1393 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 9: /* SIMPLE_ATOM_QUERY_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1399 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 10: /* COMPLEX_ATOM_QUERY_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1405 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 11: /* RINGSIZE_ATOM_QUERY_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1411 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 12: /* RINGBOND_ATOM_QUERY_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1417 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 13: /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1423 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 14: /* HYB_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1429 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 15: /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1435 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 16: /* ALIPHATIC  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1441 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 17: /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1447 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 41: /* BOND_TOKEN  */
#line 102 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1453 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 49: /* atomd  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1459 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 50: /* hydrogen_atom  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1465 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 51: /* atom_expr  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1471 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 52: /* point_query  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1477 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 53: /* recursive_query  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1483 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 54: /* atom_query  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1489 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 55: /* possible_range_query  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1495 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 56: /* simple_atom  */
#line 101 "smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1501 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 57: /* bond_expr  */
#line 102 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1507 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 58: /* bond_query  */
#line 102 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1513 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case 59: /* bondd  */
#line 102 "smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1519 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
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
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, void *scanner, int& start_token)
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
#line 110 "smarts.yy"
              {
// the molList has already been updated, no need to do anything
}
#line 1795 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 3:
#line 113 "smarts.yy"
                             {
  lastAtom = (yyvsp[-1].atom);
  YYACCEPT;
}
#line 1804 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 4:
#line 117 "smarts.yy"
                          {
  YYABORT;
}
#line 1812 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 5:
#line 120 "smarts.yy"
             {
  YYABORT;
}
#line 1820 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 6:
#line 123 "smarts.yy"
                                 {
  lastBond = (yyvsp[-1].bond);
  YYACCEPT;
}
#line 1829 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 7:
#line 127 "smarts.yy"
                       {
  delete (yyvsp[0].bond);
  YYABORT;
}
#line 1838 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 8:
#line 131 "smarts.yy"
             {
  YYABORT;
}
#line 1846 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 9:
#line 134 "smarts.yy"
                            {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1856 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 10:
#line 139 "smarts.yy"
                       {
  YYACCEPT;
}
#line 1864 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 11:
#line 142 "smarts.yy"
                  {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1874 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 15:
#line 153 "smarts.yy"
            {
  delete (yyvsp[0].atom);
  YYABORT;
}
#line 1883 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 16:
#line 161 "smarts.yy"
           {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (*molList)[ sz ]->addAtom((yyvsp[0].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
}
#line 1896 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 17:
#line 169 "smarts.yy"
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
  newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond(newB);
  delete newB;
  //delete $2;
}
#line 1926 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 18:
#line 195 "smarts.yy"
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
  (yyvsp[-1].bond)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond((yyvsp[-1].bond));
  delete (yyvsp[-1].bond);
}
#line 1951 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 19:
#line 216 "smarts.yy"
                            {
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1960 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
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
  if(!(mp->getAllBondsWithBookmark((yyvsp[0].ival)).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }
  mp->setAtomBookmark(atom,(yyvsp[0].ival));

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  if(atom->hasProp(RDKit::common_properties::_RingClosures)){
    atom->getProp(RDKit::common_properties::_RingClosures,tmp);
  }
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);

}
#line 1998 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 21:
#line 255 "smarts.yy"
                            {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();

  mp->setBondBookmark((yyvsp[-1].bond),(yyvsp[0].ival));
  (yyvsp[-1].bond)->setOwningMol(mp);
  (yyvsp[-1].bond)->setBeginAtomIdx(atom->getIdx());
  (yyvsp[-1].bond)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->setAtomBookmark(atom,(yyvsp[0].ival));

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  if(atom->hasProp(RDKit::common_properties::_RingClosures)){
    atom->getProp(RDKit::common_properties::_RingClosures,tmp);
  }
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);

}
#line 2023 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 22:
#line 276 "smarts.yy"
             {
  RWMol *m1_p = (*molList)[(yyval.moli)],*m2_p=(*molList)[(yyvsp[0].moli)];
  unsigned int origNumAts = m1_p->getNumAtoms();
  Atom *active = m1_p->getActiveAtom();
  // FIX: handle generic bonds here
  SmilesParseOps::AddFragToMol(m1_p,m2_p,Bond::UNSPECIFIED,Bond::NONE);
  delete m2_p;
  Bond *bond = m1_p->getBondBetweenAtoms(active->getIdx(),origNumAts);
  if(bond){
    bond->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }
  int sz = molList->size();
  if ( sz==(yyvsp[0].moli)+1) {
    molList->resize( sz-1 );
  }
}
#line 2044 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 23:
#line 295 "smarts.yy"
                                               { (yyval.moli) = (yyvsp[-1].moli); }
#line 2050 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 24:
#line 296 "smarts.yy"
                                                   {
  // FIX: this needs to handle arbitrary bond_exprs
  (yyval.moli) = (yyvsp[-1].moli);
  int sz     = molList->size();
  (yyvsp[-2].bond)->setOwningMol((*molList)[ sz-1 ]);
  (yyvsp[-2].bond)->setBeginAtomIdx(0);
  (*molList)[ sz-1 ]->setBondBookmark((yyvsp[-2].bond),ci_LEADING_BOND);
}
#line 2063 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 27:
#line 310 "smarts.yy"
{
  (yyval.atom) = (yyvsp[-1].atom);
}
#line 2071 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 28:
#line 314 "smarts.yy"
{
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 2080 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 29:
#line 337 "smarts.yy"
{
  (yyval.atom) = new QueryAtom(1);
}
#line 2088 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 30:
#line 341 "smarts.yy"
{
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 2097 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 31:
#line 345 "smarts.yy"
                                                  {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2108 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 32:
#line 351 "smarts.yy"
                                                                     {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-4].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-4].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2121 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 33:
#line 360 "smarts.yy"
                                                       {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2132 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 34:
#line 366 "smarts.yy"
                                                                          {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2145 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 35:
#line 374 "smarts.yy"
                                                              {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-3].ival));
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 2158 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 36:
#line 382 "smarts.yy"
                                                                                 {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-5].ival));
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-5].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 2173 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 37:
#line 395 "smarts.yy"
                                         {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 2184 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 38:
#line 401 "smarts.yy"
                               {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  (yyvsp[-2].atom)->setAtomicNum(0);
  delete (yyvsp[0].atom);
}
#line 2196 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 39:
#line 408 "smarts.yy"
                                 {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 2207 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 40:
#line 414 "smarts.yy"
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
#line 2235 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 42:
#line 440 "smarts.yy"
                                   {
  (yyvsp[0].atom)->getQuery()->setNegation(!((yyvsp[0].atom)->getQuery()->getNegation()));
  (yyvsp[0].atom)->setAtomicNum(0);
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[0].atom));
  (yyval.atom) = (yyvsp[0].atom);
}
#line 2246 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 45:
#line 451 "smarts.yy"
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
#line 2268 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 46:
#line 468 "smarts.yy"
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
#line 2294 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 48:
#line 493 "smarts.yy"
                     {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2304 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 50:
#line 499 "smarts.yy"
                    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2314 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 51:
#line 504 "smarts.yy"
                    { (yyval.atom) = new QueryAtom((yyvsp[0].ival)); }
#line 2320 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 52:
#line 505 "smarts.yy"
                           {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsotope((yyvsp[-2].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
}
#line 2330 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 59:
#line 516 "smarts.yy"
                                  {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
}
#line 2338 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 60:
#line 519 "smarts.yy"
                                         {
  (yyvsp[-1].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2346 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 61:
#line 522 "smarts.yy"
                                                  {
  (yyvsp[-1].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2354 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 62:
#line 525 "smarts.yy"
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
}
#line 2362 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 63:
#line 528 "smarts.yy"
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
}
#line 2370 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 64:
#line 531 "smarts.yy"
                                     {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
}
#line 2378 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 65:
#line 534 "smarts.yy"
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>((yyvsp[-1].ival),oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2389 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 66:
#line 540 "smarts.yy"
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>((yyvsp[-2].ival),oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2400 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 67:
#line 546 "smarts.yy"
                                                                                    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
}
#line 2412 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 68:
#line 553 "smarts.yy"
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)));
  newQ->setIsotope((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2425 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 69:
#line 561 "smarts.yy"
                        {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)));
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomHCountQuery((yyvsp[0].ival)),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2438 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 70:
#line 569 "smarts.yy"
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2449 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 71:
#line 575 "smarts.yy"
          {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2460 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 72:
#line 581 "smarts.yy"
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  newQ->setFormalCharge((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2471 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 73:
#line 587 "smarts.yy"
                    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
#line 2482 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 74:
#line 593 "smarts.yy"
           {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
#line 2493 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 75:
#line 599 "smarts.yy"
                  {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[0].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,0);
  (yyval.atom)=newQ;
}
#line 2505 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 76:
#line 606 "smarts.yy"
                         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[-1].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,(yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2517 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 78:
#line 614 "smarts.yy"
         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 2527 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 80:
#line 622 "smarts.yy"
                                  {
  (yyvsp[0].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
}
#line 2535 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 81:
#line 625 "smarts.yy"
                                           {
  (yyvsp[0].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
}
#line 2543 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 82:
#line 628 "smarts.yy"
                            {
  (yyvsp[0].atom)->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
}
#line 2551 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 83:
#line 631 "smarts.yy"
                            {
  (yyvsp[0].atom)->setQuery(makeAtomRingBondCountQuery(0));
}
#line 2559 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 84:
#line 634 "smarts.yy"
                              {
  (yyvsp[0].atom)->setQuery(makeAtomImplicitHCountQuery(0));
}
#line 2567 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 85:
#line 637 "smarts.yy"
             {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
#line 2577 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 86:
#line 642 "smarts.yy"
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNegativeFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
#line 2587 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 87:
#line 650 "smarts.yy"
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
#line 2604 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 88:
#line 662 "smarts.yy"
                      {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->setQuery(makeAtomTypeQuery((yyvsp[0].ival),true));
}
#line 2614 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 90:
#line 672 "smarts.yy"
                                        {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2623 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 91:
#line 676 "smarts.yy"
                               {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
}
#line 2632 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 92:
#line 680 "smarts.yy"
                                 {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2641 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 95:
#line 688 "smarts.yy"
                   {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2650 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 97:
#line 696 "smarts.yy"
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
#line 2661 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 98:
#line 702 "smarts.yy"
             {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
#line 2672 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 99:
#line 708 "smarts.yy"
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
#line 2683 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 100:
#line 714 "smarts.yy"
           {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
#line 2693 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 101:
#line 719 "smarts.yy"
                  {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
#line 2702 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 102:
#line 726 "smarts.yy"
                                   { (yyval.ival)=2; }
#line 2708 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 103:
#line 727 "smarts.yy"
                    { (yyval.ival)=(yyvsp[0].ival); }
#line 2714 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 104:
#line 728 "smarts.yy"
             { (yyval.ival)=1; }
#line 2720 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 105:
#line 729 "smarts.yy"
                          { (yyval.ival)=-2; }
#line 2726 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 106:
#line 730 "smarts.yy"
                     { (yyval.ival)=-(yyvsp[0].ival); }
#line 2732 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 107:
#line 731 "smarts.yy"
              { (yyval.ival)=-1; }
#line 2738 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 109:
#line 736 "smarts.yy"
                                          { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 2744 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 110:
#line 737 "smarts.yy"
                                                         { (yyval.ival) = (yyvsp[-1].ival); }
#line 2750 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 111:
#line 738 "smarts.yy"
                                                               { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2756 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 112:
#line 739 "smarts.yy"
                                                                     { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2762 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 113:
#line 740 "smarts.yy"
                                                                           { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2768 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 114:
#line 741 "smarts.yy"
                                                                                 { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2774 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 118:
#line 752 "smarts.yy"
                       { 
    if((yyvsp[-1].ival) >= std::numeric_limits<std::int32_t>::max()/10 || 
     (yyvsp[-1].ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yyvsp[0].ival) ){
     yysmarts_error(input,molList,lastAtom,lastBond,numAtomsParsed,numBondsParsed,scanner,start_token,"number too large");
     YYABORT;
  }
  (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
#line 2786 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;


#line 2790 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
      yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token, YY_("syntax error"));
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
        yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token, yymsgp);
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
                      yytoken, &yylval, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token);
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
                  yystos[yystate], yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token);
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
  yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token, YY_("memory exhausted"));
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
                  yytoken, &yylval, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[+*yyssp], yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, scanner, start_token);
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
#line 765 "smarts.yy"

