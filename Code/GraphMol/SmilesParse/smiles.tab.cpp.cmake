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
#define yyparse         yysmiles_parse
#define yylex           yysmiles_lex
#define yyerror         yysmiles_error
#define yydebug         yysmiles_debug
#define yynerrs         yysmiles_nerrs

/* First part of user prologue.  */
#line 1 "smiles.yy"


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
#include <limits>

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
      SmilesParseOps::CleanupAfterParseError(*iter);
      delete *iter;
  }
  molList->clear();
  molList->resize(0);
 }
}
void
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                RDKit::Atom* &,
                RDKit::Bond* &,
                unsigned int &,unsigned int &,
                std::list<unsigned int> *,
		void *,int, const char * msg )
{
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMILES Parse Error: " << msg << " while parsing: " << input << std::endl;
}

void
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::list<unsigned int> *,
		void *,int, const char * msg )
{
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMILES Parse Error: " << msg << " while parsing: " << input << std::endl;
}



#line 138 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"

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
#line 82 "smiles.yy"

  int                      moli;
  RDKit::Atom * atom;
  RDKit::Bond * bond;
  int                      ival;

#line 226 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int yysmiles_parse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token);
/* "%code provides" blocks.  */
#line 77 "smiles.yy"

#define YY_DECL int yylex \
               (YYSTYPE * yylval_param , yyscan_t yyscanner, int& start_token)

#line 243 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"

#endif /* !YY_YYSMILES_SCRATCH_RDKIT_GIT_CODE_GRAPHMOL_SMILESPARSE_SMILES_TAB_HPP_INCLUDED  */



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
typedef yytype_int8 yy_state_t;

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
#define YYFINAL  33
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   149

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  29
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  15
/* YYNRULES -- Number of rules.  */
#define YYNRULES  71
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  105

#define YYUNDEFTOK  2
#define YYMAXUTOK   283


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
      25,    26,    27,    28
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   112,   112,   115,   119,   122,   126,   130,   133,   138,
     141,   149,   150,   151,   152,   160,   171,   182,   203,   212,
     218,   239,   263,   283,   294,   316,   325,   338,   339,   345,
     346,   352,   360,   361,   362,   363,   364,   365,   366,   370,
     371,   372,   373,   374,   375,   376,   377,   378,   382,   383,
     384,   388,   389,   390,   391,   392,   393,   397,   398,   402,
     403,   404,   405,   406,   407,   408,   412,   413,   417,   418,
     429,   430
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
  "ATOM_CLOSE_TOKEN", "EOS_TOKEN", "$accept", "meta_start", "bad_atom_def",
  "mol", "bondd", "atomd", "charge_element", "h_element", "chiral_element",
  "element", "simple_atom", "ring_number", "number", "nonzero_number",
  "digit", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_int16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283
};
# endif

#define YYPACT_NINF (-30)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-30)

#define yytable_value_is_error(Yyn) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
      62,   -10,    11,    51,     7,     1,   -30,   -30,   -30,   106,
      82,   -30,   -30,   -30,   -30,   -30,    -6,    76,    -4,    76,
      76,   -30,    -1,   -30,    26,    31,    41,    48,   113,    96,
     -30,   -30,    70,   -30,    73,   -30,   -20,   -30,   -30,   -30,
     103,   -30,    11,     2,    27,     2,   -30,   -30,   -30,    -4,
      76,   -30,   -30,   -30,   -20,   -30,   -30,   125,    30,    -4,
      58,   -30,    49,    -4,   -30,   -30,   -30,   -30,    -4,   -30,
      11,    11,   -30,   -30,   -30,   -30,    96,    96,   -30,   -30,
     -30,   -30,   -30,   -30,   -30,   -30,   -30,    -4,   -30,    83,
     -30,   -30,   -30,   115,   -30,   -30,   -30,   129,   -30,   133,
     -30,   137,   -30,   105,   -30
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     0,     0,     7,     0,    10,    58,    57,     0,
       2,    15,    29,    53,    68,    66,    39,     0,     0,     0,
       0,     4,     0,    14,    32,    45,    48,    51,     0,    67,
      28,    27,     6,     1,     0,     9,     0,    51,    70,    71,
       0,    26,     0,     0,     0,     0,    16,    20,    59,    41,
       0,    13,    55,    11,    14,    12,     3,    36,    33,    46,
      49,    54,    40,     0,    52,    69,     5,     8,     0,    31,
       0,     0,    23,    19,    18,    22,     0,     0,    17,    21,
      43,    37,    38,    34,    35,    47,    50,    42,    56,     0,
      25,    24,    60,     0,    44,    30,    61,     0,    62,     0,
      63,     0,    64,     0,    65
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -30,   -30,    13,   -30,   -30,    10,    12,   -30,   -30,   -30,
       6,    44,   -14,   -30,   -29
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     5,    53,    10,    32,    11,    23,    24,    25,    26,
      12,    47,    28,    29,    48
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      65,    33,    34,    68,    52,    14,    15,    69,     7,    27,
       8,    38,    39,    22,    49,    37,    21,     7,     6,     8,
      46,    36,    30,    37,    44,    37,    37,    56,     9,    35,
      51,    54,    31,    55,    64,    80,    76,     9,    77,    14,
      15,    57,    58,    82,    84,    85,    83,    92,    93,    88,
      72,    59,    73,    74,    89,    78,    37,     7,    13,     8,
      14,    15,    60,     1,    97,     2,     3,     4,    99,    87,
     101,    16,   103,    94,    17,    18,   -29,    19,    20,    86,
      90,    91,     7,    13,     8,    14,    15,    75,     7,    79,
       8,    38,    39,    40,    41,    42,    16,    43,    66,    17,
      18,    67,    50,    20,    44,    38,    39,    45,     9,     7,
      95,     8,     7,    13,     8,    14,    15,   104,    70,     7,
      61,     8,     0,     0,    38,    39,    16,    96,    71,     9,
      18,     0,     0,    62,    14,    15,     0,    63,    38,    39,
      81,    98,    38,    39,     0,   100,    38,    39,     0,   102
};

static const yytype_int8 yycheck[] =
{
      29,     0,     1,    23,    18,     9,    10,    27,     6,     3,
       8,     9,    10,     3,    20,     9,     3,     6,    28,     8,
      10,     9,    15,    17,    22,    19,    20,    28,    26,    28,
      17,    19,    25,    20,    28,    49,     9,    26,    11,     9,
      10,    15,    16,    57,    58,    59,    16,    76,    77,    63,
      40,    20,    42,    43,    68,    45,    50,     6,     7,     8,
       9,    10,    21,     1,    93,     3,     4,     5,    97,    20,
      99,    20,   101,    87,    23,    24,    28,    26,    27,    21,
      70,    71,     6,     7,     8,     9,    10,    43,     6,    45,
       8,     9,    10,    11,    12,    13,    20,    15,    28,    23,
      24,    28,    26,    27,    22,     9,    10,    25,    26,     6,
      27,     8,     6,     7,     8,     9,    10,    12,    15,     6,
       7,     8,    -1,    -1,     9,    10,    20,    12,    25,    26,
      24,    -1,    -1,    20,     9,    10,    -1,    24,     9,    10,
      15,    12,     9,    10,    -1,    12,     9,    10,    -1,    12
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,    30,    28,     6,     8,    26,
      32,    34,    39,     7,     9,    10,    20,    23,    24,    26,
      27,    31,    34,    35,    36,    37,    38,    39,    41,    42,
      15,    25,    33,     0,     1,    28,    35,    39,     9,    10,
      11,    12,    13,    15,    22,    25,    34,    40,    43,    20,
      26,    31,    41,    31,    35,    31,    28,    15,    16,    20,
      21,     7,    20,    24,    39,    43,    28,    28,    23,    27,
      15,    25,    34,    34,    34,    40,     9,    11,    34,    40,
      41,    15,    41,    16,    41,    41,    21,    20,    41,    41,
      34,    34,    43,    43,    41,    27,    12,    43,    12,    43,
      12,    43,    12,    43,    12
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_int8 yyr1[] =
{
       0,    29,    30,    30,    30,    30,    30,    30,    30,    30,
      30,    31,    31,    31,    31,    32,    32,    32,    32,    32,
      32,    32,    32,    32,    32,    32,    32,    33,    33,    34,
      34,    34,    35,    35,    35,    35,    35,    35,    35,    36,
      36,    36,    36,    36,    36,    36,    36,    36,    37,    37,
      37,    38,    38,    38,    38,    38,    38,    39,    39,    40,
      40,    40,    40,    40,    40,    40,    41,    41,    42,    42,
      43,    43
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     3,     2,     3,     2,     1,     3,     2,
       2,     2,     2,     2,     1,     1,     2,     3,     3,     3,
       2,     3,     3,     3,     4,     4,     2,     1,     1,     1,
       5,     3,     1,     2,     3,     3,     2,     3,     3,     1,
       2,     2,     3,     3,     4,     1,     2,     3,     1,     2,
       3,     1,     2,     1,     2,     2,     3,     1,     1,     1,
       3,     4,     5,     6,     7,     8,     1,     1,     1,     2,
       1,     1
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
        yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, YY_("syntax error: cannot back up")); \
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
                  Type, Value, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (numAtomsParsed);
  YYUSE (numBondsParsed);
  YYUSE (branchPoints);
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
yy_symbol_print (FILE *yyo, int yytype, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
{
  YYFPRINTF (yyo, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyo, yytype, yyvaluep, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
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
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp, int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
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
                                              , input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
{
  YYUSE (yyvaluep);
  YYUSE (input);
  YYUSE (molList);
  YYUSE (lastAtom);
  YYUSE (lastBond);
  YYUSE (numAtomsParsed);
  YYUSE (numBondsParsed);
  YYUSE (branchPoints);
  YYUSE (scanner);
  YYUSE (start_token);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yytype)
    {
    case 6: /* AROMATIC_ATOM_TOKEN  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1247 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 7: /* ATOM_TOKEN  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1253 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 8: /* ORGANIC_ATOM_TOKEN  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1259 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 25: /* BOND_TOKEN  */
#line 104 "smiles.yy"
            { delete ((*yyvaluep).bond); }
#line 1265 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 33: /* bondd  */
#line 104 "smiles.yy"
            { delete ((*yyvaluep).bond); }
#line 1271 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 34: /* atomd  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1277 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 35: /* charge_element  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1283 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 36: /* h_element  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1289 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 37: /* chiral_element  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1295 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 38: /* element  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1301 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
        break;

    case 39: /* simple_atom  */
#line 103 "smiles.yy"
            { delete ((*yyvaluep).atom); }
#line 1307 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
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
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
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
#line 112 "smiles.yy"
              {
// the molList has already been updated, no need to do anything
}
#line 1583 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 3:
#line 115 "smiles.yy"
                             {
  lastAtom = (yyvsp[-1].atom);
  YYACCEPT;
}
#line 1592 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 4:
#line 119 "smiles.yy"
                          {
  YYABORT;
}
#line 1600 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 5:
#line 122 "smiles.yy"
                             {
  lastBond = (yyvsp[-1].bond);
  YYACCEPT;
}
#line 1609 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 6:
#line 126 "smiles.yy"
                   {
  delete (yyvsp[0].bond);
  YYABORT;
}
#line 1618 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 7:
#line 130 "smiles.yy"
             {
  YYABORT;
}
#line 1626 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 8:
#line 133 "smiles.yy"
                            {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1636 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 9:
#line 138 "smiles.yy"
                       {
  YYACCEPT;
}
#line 1644 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 10:
#line 141 "smiles.yy"
                  {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1654 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 14:
#line 152 "smiles.yy"
                 {
  delete (yyvsp[0].atom);
  YYABORT;
}
#line 1663 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 15:
#line 160 "smiles.yy"
           {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  RDKit::RWMol *curMol = (*molList)[ sz ];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1);
  curMol->addAtom((yyvsp[0].atom), true, true);
  //delete $1;
  (yyval.moli) = sz;
}
#line 1678 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 16:
#line 171 "smiles.yy"
                  {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $2;
}
#line 1693 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 17:
#line 182 "smiles.yy"
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
  mp->addBond((yyvsp[-1].bond),true);
  //delete $3;
}
#line 1718 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 18:
#line 203 "smiles.yy"
                        {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $3;
}
#line 1731 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 19:
#line 212 "smiles.yy"
                            {
  RWMol *mp = (*molList)[(yyval.moli)];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1741 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 20:
#line 218 "smiles.yy"
                  {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  mp->setAtomBookmark(atom,(yyvsp[0].ival));

  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::UNSPECIFIED);
  mp->setBondBookmark(newB,(yyvsp[0].ival));
  newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  if(!(mp->getAllBondsWithBookmark((yyvsp[0].ival)).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
}
#line 1766 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 21:
#line 239 "smiles.yy"
                             {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     (yyvsp[-1].bond)->getBondType());
  if((yyvsp[-1].bond)->hasProp(RDKit::common_properties::_unspecifiedOrder)){
    newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  }
  newB->setBondDir((yyvsp[-1].bond)->getBondDir());
  mp->setAtomBookmark(atom,(yyvsp[0].ival));
  mp->setBondBookmark(newB,(yyvsp[0].ival));
  if(!(mp->getAllBondsWithBookmark((yyvsp[0].ival)).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
  delete (yyvsp[-1].bond);
}
#line 1794 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 22:
#line 263 "smiles.yy"
                              {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::SINGLE);
  newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->setAtomBookmark(atom,(yyvsp[0].ival));
  mp->setBondBookmark(newB,(yyvsp[0].ival));
  if(!(mp->getAllBondsWithBookmark((yyvsp[0].ival)).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yyvsp[0].ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
}
#line 1818 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 23:
#line 283 "smiles.yy"
                             {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $3;
  branchPoints->push_back(atomIdx1);
}
#line 1834 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 24:
#line 294 "smiles.yy"
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
  mp->addBond((yyvsp[-1].bond),true);

  //delete $4;
  branchPoints->push_back(atomIdx1);
}
#line 1861 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 25:
#line 316 "smiles.yy"
                                         {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $4;
  branchPoints->push_back(atomIdx1);
}
#line 1875 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 26:
#line 325 "smiles.yy"
                        {
  if(branchPoints->empty()){
     yyerror(input,molList,branchPoints,scanner,start_token,"extra close parentheses");
     yyErrorCleanup(molList);
     YYABORT;
  }
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->setActiveAtom(branchPoints->back());
  branchPoints->pop_back();
}
#line 1890 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 28:
#line 339 "smiles.yy"
                        {
          (yyval.bond) = new Bond(Bond::SINGLE);
          }
#line 1898 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 30:
#line 347 "smiles.yy"
{
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setNoImplicit(true);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1908 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 31:
#line 353 "smiles.yy"
{
  (yyval.atom) = (yyvsp[-1].atom);
  (yyvsp[-1].atom)->setNoImplicit(true);
}
#line 1917 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 33:
#line 361 "smiles.yy"
                       { (yyvsp[-1].atom)->setFormalCharge(1); }
#line 1923 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 34:
#line 362 "smiles.yy"
                                  { (yyvsp[-2].atom)->setFormalCharge(2); }
#line 1929 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 35:
#line 363 "smiles.yy"
                              { (yyvsp[-2].atom)->setFormalCharge((yyvsp[0].ival)); }
#line 1935 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 36:
#line 364 "smiles.yy"
                        { (yyvsp[-1].atom)->setFormalCharge(-1); }
#line 1941 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 37:
#line 365 "smiles.yy"
                                    { (yyvsp[-2].atom)->setFormalCharge(-2); }
#line 1947 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 38:
#line 366 "smiles.yy"
                               { (yyvsp[-2].atom)->setFormalCharge(-(yyvsp[0].ival)); }
#line 1953 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 39:
#line 370 "smiles.yy"
                        { (yyval.atom) = new Atom(1); }
#line 1959 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 40:
#line 371 "smiles.yy"
                                 { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-1].ival)); }
#line 1965 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 41:
#line 372 "smiles.yy"
                                  { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs(1); }
#line 1971 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 42:
#line 373 "smiles.yy"
                                         { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-2].ival)); (yyval.atom)->setNumExplicitHs(1);}
#line 1977 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 43:
#line 374 "smiles.yy"
                                         { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs((yyvsp[0].ival)); }
#line 1983 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 44:
#line 375 "smiles.yy"
                                                { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-3].ival)); (yyval.atom)->setNumExplicitHs((yyvsp[0].ival));}
#line 1989 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 46:
#line 377 "smiles.yy"
                                                        { (yyval.atom) = (yyvsp[-1].atom); (yyvsp[-1].atom)->setNumExplicitHs(1);}
#line 1995 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 47:
#line 378 "smiles.yy"
                                                { (yyval.atom) = (yyvsp[-2].atom); (yyvsp[-2].atom)->setNumExplicitHs((yyvsp[0].ival));}
#line 2001 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 49:
#line 383 "smiles.yy"
                   { (yyvsp[-1].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
#line 2007 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 50:
#line 384 "smiles.yy"
                            { (yyvsp[-2].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
#line 2013 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 52:
#line 389 "smiles.yy"
                                           { (yyvsp[0].atom)->setIsotope( (yyvsp[-1].ival) ); (yyval.atom) = (yyvsp[0].atom); }
#line 2019 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 54:
#line 391 "smiles.yy"
                                                   { (yyvsp[0].atom)->setIsotope( (yyvsp[-1].ival) ); (yyval.atom) = (yyvsp[0].atom); }
#line 2025 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 55:
#line 392 "smiles.yy"
                                                 { (yyval.atom) = new Atom((yyvsp[0].ival)); }
#line 2031 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 56:
#line 393 "smiles.yy"
                                                         { (yyval.atom) = new Atom((yyvsp[0].ival)); (yyval.atom)->setIsotope((yyvsp[-2].ival)); }
#line 2037 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 60:
#line 403 "smiles.yy"
                                          { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 2043 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 61:
#line 404 "smiles.yy"
                                                         { (yyval.ival) = (yyvsp[-1].ival); }
#line 2049 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 62:
#line 405 "smiles.yy"
                                                               { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2055 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 63:
#line 406 "smiles.yy"
                                                                     { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2061 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 64:
#line 407 "smiles.yy"
                                                                           { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2067 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 65:
#line 408 "smiles.yy"
                                                                                 { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2073 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 69:
#line 418 "smiles.yy"
                       { 
  if((yyvsp[-1].ival) >= std::numeric_limits<std::int32_t>::max()/10 || 
     (yyvsp[-1].ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yyvsp[0].ival) ){
     yyerror(input,molList,branchPoints,scanner,start_token,"number too large");
     yyErrorCleanup(molList);
     YYABORT;
  }
  (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); 
  }
#line 2087 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;


#line 2091 "/scratch/RDKit_git/Code/GraphMol/SmilesParse/smiles.tab.cpp"

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
      yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, YY_("syntax error"));
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
        yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, yymsgp);
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
                      yytoken, &yylval, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
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
                  yystos[yystate], yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
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
  yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, YY_("memory exhausted"));
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
                  yytoken, &yylval, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[+*yyssp], yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
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
#line 440 "smiles.yy"

