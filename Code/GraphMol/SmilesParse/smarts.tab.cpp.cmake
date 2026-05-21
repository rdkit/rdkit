/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
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
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

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


  //
  //  Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors
  //
  //   @@ All Rights Reserved  @@
  //
#include <cstring>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <RDGeneral/RDLog.h>

#define YYDEBUG 1
#include "smarts.tab.hpp"

extern int yysmarts_lex(YYSTYPE *,void *, int &, unsigned int&);

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
  const std::uint64_t SMARTS_H_MASK = 0x1;
  const std::uint64_t SMARTS_CHARGE_MASK = 0x2;

  void atom_expr_and_point_query(QueryAtom *atom_expr, QueryAtom *point_query) {
    atom_expr->expandQuery(point_query->getQuery()->copy(), Queries::COMPOSITE_AND, true);
    if (atom_expr->getChiralTag() == Atom::CHI_UNSPECIFIED) {
      atom_expr->setChiralTag(point_query->getChiralTag());
      int perm;
      if (point_query->getPropIfPresent(common_properties::_chiralPermutation, perm)) {
        atom_expr->setProp(common_properties::_chiralPermutation, perm);
      }
    }
    if (point_query->getFlags() & SMARTS_H_MASK) {
      if (!(atom_expr->getFlags() & SMARTS_H_MASK)) {
        atom_expr->setNumExplicitHs(point_query->getNumExplicitHs());
        atom_expr->setNoImplicit(true);
        atom_expr->getFlags() |= SMARTS_H_MASK;
      } else if (atom_expr->getNumExplicitHs() != point_query->getNumExplicitHs()) {
        // conflicting queries...
        atom_expr->setNumExplicitHs(0);
        atom_expr->setNoImplicit(true);
      }
    }
    if (point_query->getFlags() & SMARTS_CHARGE_MASK) {
      if (!(atom_expr->getFlags() & SMARTS_CHARGE_MASK)) {
        atom_expr->setFormalCharge(point_query->getFormalCharge());
        atom_expr->getFlags() |= SMARTS_CHARGE_MASK;
      } else if (atom_expr->getFormalCharge() != point_query->getFormalCharge()) {
        // conflicting queries...
        atom_expr->setFormalCharge(0);
      }
    }
  }

}
void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                RDKit::Atom* &,
                RDKit::Bond* &,
                unsigned int &,
                unsigned int &,
                std::vector<std::pair<unsigned int, unsigned int>>&,
                void *,
                int,
                unsigned int bad_token_position,
                const char *msg  )
{
  yyErrorCleanup(ms);
  SmilesParseOps::detail::printSyntaxErrorMessage(input,
                                                  msg,
                                                  bad_token_position,
                                                  "SMARTS");
}

void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::vector<std::pair<unsigned int, unsigned int>>&,
                void *,
                int,
                unsigned int bad_token_position,
                const char * msg )
{
  yyErrorCleanup(ms);
  SmilesParseOps::detail::printSyntaxErrorMessage(input,
                                                  msg,
                                                  bad_token_position,
                                                  "SMARTS");
}

void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                unsigned int bad_token_position, const char * msg )
{
  yyErrorCleanup(ms);
  SmilesParseOps::detail::printSyntaxErrorMessage(input,
                                                  msg,
                                                  bad_token_position,
                                                  "SMARTS");
}




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

#include "smarts.tab.hpp"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_START_MOL = 3,                  /* START_MOL  */
  YYSYMBOL_START_ATOM = 4,                 /* START_ATOM  */
  YYSYMBOL_START_BOND = 5,                 /* START_BOND  */
  YYSYMBOL_AROMATIC_ATOM_TOKEN = 6,        /* AROMATIC_ATOM_TOKEN  */
  YYSYMBOL_ORGANIC_ATOM_TOKEN = 7,         /* ORGANIC_ATOM_TOKEN  */
  YYSYMBOL_ATOM_TOKEN = 8,                 /* ATOM_TOKEN  */
  YYSYMBOL_SIMPLE_ATOM_QUERY_TOKEN = 9,    /* SIMPLE_ATOM_QUERY_TOKEN  */
  YYSYMBOL_COMPLEX_ATOM_QUERY_TOKEN = 10,  /* COMPLEX_ATOM_QUERY_TOKEN  */
  YYSYMBOL_MIN_RINGSIZE_ATOM_QUERY_TOKEN = 11, /* MIN_RINGSIZE_ATOM_QUERY_TOKEN  */
  YYSYMBOL_RINGSIZE_ATOM_QUERY_TOKEN = 12, /* RINGSIZE_ATOM_QUERY_TOKEN  */
  YYSYMBOL_RINGBOND_ATOM_QUERY_TOKEN = 13, /* RINGBOND_ATOM_QUERY_TOKEN  */
  YYSYMBOL_IMPLICIT_H_ATOM_QUERY_TOKEN = 14, /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
  YYSYMBOL_HYB_TOKEN = 15,                 /* HYB_TOKEN  */
  YYSYMBOL_HETERONEIGHBOR_ATOM_QUERY_TOKEN = 16, /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
  YYSYMBOL_ALIPHATIC = 17,                 /* ALIPHATIC  */
  YYSYMBOL_ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN = 18, /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
  YYSYMBOL_ZERO_TOKEN = 19,                /* ZERO_TOKEN  */
  YYSYMBOL_NONZERO_DIGIT_TOKEN = 20,       /* NONZERO_DIGIT_TOKEN  */
  YYSYMBOL_GROUP_OPEN_TOKEN = 21,          /* GROUP_OPEN_TOKEN  */
  YYSYMBOL_GROUP_CLOSE_TOKEN = 22,         /* GROUP_CLOSE_TOKEN  */
  YYSYMBOL_SEPARATOR_TOKEN = 23,           /* SEPARATOR_TOKEN  */
  YYSYMBOL_RANGE_OPEN_TOKEN = 24,          /* RANGE_OPEN_TOKEN  */
  YYSYMBOL_RANGE_CLOSE_TOKEN = 25,         /* RANGE_CLOSE_TOKEN  */
  YYSYMBOL_HASH_TOKEN = 26,                /* HASH_TOKEN  */
  YYSYMBOL_MINUS_TOKEN = 27,               /* MINUS_TOKEN  */
  YYSYMBOL_PLUS_TOKEN = 28,                /* PLUS_TOKEN  */
  YYSYMBOL_H_TOKEN = 29,                   /* H_TOKEN  */
  YYSYMBOL_AT_TOKEN = 30,                  /* AT_TOKEN  */
  YYSYMBOL_PERCENT_TOKEN = 31,             /* PERCENT_TOKEN  */
  YYSYMBOL_ATOM_OPEN_TOKEN = 32,           /* ATOM_OPEN_TOKEN  */
  YYSYMBOL_ATOM_CLOSE_TOKEN = 33,          /* ATOM_CLOSE_TOKEN  */
  YYSYMBOL_NOT_TOKEN = 34,                 /* NOT_TOKEN  */
  YYSYMBOL_AND_TOKEN = 35,                 /* AND_TOKEN  */
  YYSYMBOL_OR_TOKEN = 36,                  /* OR_TOKEN  */
  YYSYMBOL_SEMI_TOKEN = 37,                /* SEMI_TOKEN  */
  YYSYMBOL_BEGIN_RECURSE = 38,             /* BEGIN_RECURSE  */
  YYSYMBOL_END_RECURSE = 39,               /* END_RECURSE  */
  YYSYMBOL_COLON_TOKEN = 40,               /* COLON_TOKEN  */
  YYSYMBOL_UNDERSCORE_TOKEN = 41,          /* UNDERSCORE_TOKEN  */
  YYSYMBOL_BOND_TOKEN = 42,                /* BOND_TOKEN  */
  YYSYMBOL_CHI_CLASS_TOKEN = 43,           /* CHI_CLASS_TOKEN  */
  YYSYMBOL_BAD_CHARACTER = 44,             /* BAD_CHARACTER  */
  YYSYMBOL_EOS_TOKEN = 45,                 /* EOS_TOKEN  */
  YYSYMBOL_YYACCEPT = 46,                  /* $accept  */
  YYSYMBOL_meta_start = 47,                /* meta_start  */
  YYSYMBOL_bad_atom_def = 48,              /* bad_atom_def  */
  YYSYMBOL_mol = 49,                       /* mol  */
  YYSYMBOL_atomd = 50,                     /* atomd  */
  YYSYMBOL_hydrogen_atom = 51,             /* hydrogen_atom  */
  YYSYMBOL_atom_expr = 52,                 /* atom_expr  */
  YYSYMBOL_point_query = 53,               /* point_query  */
  YYSYMBOL_recursive_query = 54,           /* recursive_query  */
  YYSYMBOL_atom_query = 55,                /* atom_query  */
  YYSYMBOL_possible_range_query = 56,      /* possible_range_query  */
  YYSYMBOL_simple_atom = 57,               /* simple_atom  */
  YYSYMBOL_bond_expr = 58,                 /* bond_expr  */
  YYSYMBOL_bond_query = 59,                /* bond_query  */
  YYSYMBOL_bondd = 60,                     /* bondd  */
  YYSYMBOL_charge_spec = 61,               /* charge_spec  */
  YYSYMBOL_ring_number = 62,               /* ring_number  */
  YYSYMBOL_number = 63,                    /* number  */
  YYSYMBOL_nonzero_number = 64,            /* nonzero_number  */
  YYSYMBOL_digit = 65,                     /* digit  */
  YYSYMBOL_branch_open_token = 66          /* branch_open_token  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




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

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
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
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
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

#if !defined yyoverflow

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
#endif /* !defined yyoverflow */

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
#define YYFINAL  57
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   606

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  46
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  21
/* YYNRULES -- Number of rules.  */
#define YYNRULES  128
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  189

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   300


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

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
       0,   190,   190,   193,   197,   200,   203,   207,   211,   214,
     220,   225,   228,   236,   237,   238,   239,   247,   256,   270,
     291,   297,   321,   342,   359,   382,   396,   397,   398,   402,
     425,   429,   434,   440,   449,   456,   465,   474,   488,   495,
     503,   510,   515,   520,   523,   529,   530,   534,   551,   575,
     576,   581,   582,   587,   588,   593,   594,   595,   596,   597,
     598,   599,   600,   604,   608,   612,   616,   620,   624,   628,
     635,   642,   651,   660,   669,   679,   689,   699,   708,   716,
     723,   729,   735,   742,   756,   757,   764,   765,   769,   773,
     777,   781,   785,   790,   798,   810,   815,   820,   825,   830,
     835,   838,   839,   847,   848,   854,   860,   866,   871,   878,
     879,   880,   881,   882,   883,   887,   888,   889,   890,   891,
     892,   893,   898,   899,   903,   904,   913,   914,   918
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "START_MOL",
  "START_ATOM", "START_BOND", "AROMATIC_ATOM_TOKEN", "ORGANIC_ATOM_TOKEN",
  "ATOM_TOKEN", "SIMPLE_ATOM_QUERY_TOKEN", "COMPLEX_ATOM_QUERY_TOKEN",
  "MIN_RINGSIZE_ATOM_QUERY_TOKEN", "RINGSIZE_ATOM_QUERY_TOKEN",
  "RINGBOND_ATOM_QUERY_TOKEN", "IMPLICIT_H_ATOM_QUERY_TOKEN", "HYB_TOKEN",
  "HETERONEIGHBOR_ATOM_QUERY_TOKEN", "ALIPHATIC",
  "ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN", "ZERO_TOKEN",
  "NONZERO_DIGIT_TOKEN", "GROUP_OPEN_TOKEN", "GROUP_CLOSE_TOKEN",
  "SEPARATOR_TOKEN", "RANGE_OPEN_TOKEN", "RANGE_CLOSE_TOKEN", "HASH_TOKEN",
  "MINUS_TOKEN", "PLUS_TOKEN", "H_TOKEN", "AT_TOKEN", "PERCENT_TOKEN",
  "ATOM_OPEN_TOKEN", "ATOM_CLOSE_TOKEN", "NOT_TOKEN", "AND_TOKEN",
  "OR_TOKEN", "SEMI_TOKEN", "BEGIN_RECURSE", "END_RECURSE", "COLON_TOKEN",
  "UNDERSCORE_TOKEN", "BOND_TOKEN", "CHI_CLASS_TOKEN", "BAD_CHARACTER",
  "EOS_TOKEN", "$accept", "meta_start", "bad_atom_def", "mol", "atomd",
  "hydrogen_atom", "atom_expr", "point_query", "recursive_query",
  "atom_query", "possible_range_query", "simple_atom", "bond_expr",
  "bond_query", "bondd", "charge_spec", "ring_number", "number",
  "nonzero_number", "digit", "branch_open_token", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-61)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-94)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
      55,   -31,    46,   195,   533,     3,   -61,   -61,   -61,   -61,
     423,   516,   -61,   -61,   -61,   -61,   250,   288,   360,   397,
     435,   -61,   478,   481,   -61,   -61,    80,    19,     5,    80,
       0,   233,   271,   461,    46,   271,    80,   -61,    -9,   309,
     -61,   -61,   -61,    21,    40,   -61,    64,    92,   -61,   -61,
     -61,   533,   -61,   -61,   143,   533,   -61,   -61,    42,   -61,
     -61,   551,   157,   -61,   342,   -61,   -61,   -61,   -61,    46,
     137,   -61,   545,   -61,   -61,    35,   -61,   -61,   126,   -61,
     -61,   -61,   -61,   -61,   -61,   -61,   -61,   -61,   -61,   -61,
     -61,   271,   -61,   157,   -61,   -61,   487,   -61,   -61,   -61,
     461,   461,   461,   -61,   162,   -61,    80,    80,   -61,   -61,
     -61,   533,   533,   533,   -61,   -61,   -61,   199,    60,   -61,
      80,   -13,   -61,    80,   566,   -61,    92,    92,   -61,   -61,
     -61,    96,    80,    23,    33,   461,   -61,   385,   347,    80,
      88,   -61,   -61,   -61,    56,   163,    93,   -61,    80,    97,
     -61,    80,    43,   -61,   117,   -61,   109,   129,   101,   127,
     211,   -61,   122,   -61,   141,   -61,    80,   -61,   235,   -61,
     -61,   171,    92,   -61,   -61,   187,   -61,   -61,   183,   -61,
     273,   -61,   -61,   -61,   -61,   311,   -61,   198,   -61
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     5,     8,     0,    12,    95,    94,    96,
       0,     2,    17,    27,    26,    51,    55,    58,    59,    60,
      61,    84,    56,    57,   122,   124,     0,   114,   111,    78,
      81,     0,     0,     0,     0,     0,    82,     4,     0,    16,
      43,    45,    46,     0,    49,    79,    85,   123,   105,   104,
     107,     0,   106,   103,     7,   100,   101,     1,     0,     9,
      11,    78,     0,    49,    85,   127,   126,   128,    25,     0,
       0,    18,     0,    21,   115,     0,    62,    65,     0,    66,
      67,    68,    63,    64,    53,   112,   113,   109,   110,    77,
      80,     0,    13,    16,    14,    44,     0,    15,    83,     3,
       0,     0,     0,    41,     0,    52,     0,    75,    50,   125,
     108,     0,     0,     0,     6,   102,    10,   114,   111,    30,
       0,     0,    28,     0,    75,    20,     0,     0,    19,    22,
      23,     0,     0,     0,    47,    38,    42,    39,    40,     0,
       0,    54,    76,    97,    98,    99,     0,    34,     0,     0,
      32,     0,     0,   116,     0,    24,     0,     0,     0,     0,
       0,    31,     0,    29,     0,    36,     0,   117,     0,    72,
      73,     0,    48,    69,    70,     0,    35,    33,     0,   118,
       0,    74,    71,    37,   119,     0,   120,     0,   121
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -61,   -61,     2,   200,    29,   -61,    18,    24,   -61,   -61,
     -61,    20,    31,   -61,   -40,   -60,   145,   -10,    74,   -45,
     -61
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     5,    92,    11,    12,    13,    39,    40,    41,    42,
      43,    63,    72,    55,    56,    45,    73,    46,    47,    74,
      75
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      64,   121,   109,    57,    58,    37,    76,    77,    79,    80,
      81,   110,    82,    83,     6,   115,    84,    86,    88,    89,
     147,    64,    14,    44,    24,    25,    98,   148,    62,   -92,
      90,    14,    38,    87,    94,    54,    99,    97,    24,    25,
      71,     7,     8,   -93,     9,   104,    85,    59,    60,    93,
     157,    89,     7,     8,    14,     9,     1,    95,     2,     3,
       4,    48,    49,   103,   152,    50,   108,    10,   133,    51,
       7,     8,   105,     9,   158,    52,   165,    53,    10,    24,
      25,   153,   154,   166,   108,   -26,   103,   116,    87,    14,
     106,   111,    14,   107,   140,    14,   141,   142,   125,    24,
      25,   128,     7,     8,   130,     9,   131,    86,    88,   168,
     146,    65,    66,   149,   142,   160,    14,   103,   135,   137,
     138,    25,   156,   180,   136,    71,   161,   109,    10,   159,
     163,   111,   112,   113,   169,   185,    65,    66,   162,   167,
     187,   164,   143,   144,   145,    24,    25,   171,    24,    25,
     175,    14,   173,   132,   170,   176,   178,   126,   127,   103,
     155,   103,   103,     7,     8,    15,     9,    16,    17,    18,
      19,    20,    21,    22,   177,    23,    24,    25,   111,   112,
     113,    24,    25,    26,    27,    28,    29,    30,   114,   139,
     122,    33,   100,   101,   102,    34,   181,   123,   111,   112,
      36,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,    22,   182,    23,    24,    25,   183,   129,    24,    25,
     188,    26,    27,    28,    29,    30,    85,    31,    32,    33,
      24,    25,   172,    34,    96,    35,   174,     0,    36,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,    22,
       0,    23,    24,    25,    65,    66,     0,   179,     0,    26,
      27,    28,    61,    30,     0,    91,    32,    33,     0,    24,
      25,    34,     0,    35,   -86,     0,    36,     7,     8,    15,
       9,    16,    17,    18,    19,    20,    21,    22,     0,    23,
      24,    25,    65,    66,     0,   184,     0,    26,    27,    28,
      29,    30,     0,    91,    32,    33,     0,    24,    25,    34,
       0,    35,   -89,     0,    36,     7,     8,    15,     9,    16,
      17,    18,    19,    20,    21,    22,     0,    23,    24,    25,
      65,    66,     0,   186,     0,    26,    27,    28,    29,    30,
       0,     0,     0,    33,   100,   101,   102,    34,     7,     8,
     105,     9,    36,     7,     8,    15,     9,    16,    17,    18,
      19,    20,    21,    22,     0,    23,    24,    25,   106,     0,
       0,   124,     0,    26,    27,    28,    29,    30,     0,    24,
      25,    33,   100,   101,    78,    34,     0,     0,     0,     0,
      36,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,    22,     0,    23,    24,    25,     0,     0,     0,     0,
       0,    26,    27,    28,    29,    30,    24,    25,     0,    33,
     100,   -90,     0,    34,     0,     0,     0,     0,    36,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,    22,
       0,    23,    24,    25,     0,     0,     0,     0,     0,    26,
      27,    28,    61,    30,    24,    25,     0,    33,     0,   -91,
       0,    34,     0,     0,     0,     0,    36,     7,     8,    15,
       9,    16,    17,    18,    19,    20,    21,    22,     0,    23,
      24,    25,     0,     0,     0,     0,     0,    26,    27,    28,
      29,    30,     0,     7,     8,    33,     9,    24,    25,    34,
      24,    25,   -87,     0,    36,   -88,    65,    66,    67,    68,
      69,     0,     0,    48,    49,     0,     0,    50,    70,    10,
       0,    51,     7,     8,     0,     9,   134,    52,     0,    53,
       0,     0,     0,     0,     0,    65,    66,    67,    68,    69,
       0,     0,    48,    49,     0,     0,    50,    70,    10,     0,
      51,     7,     8,     0,     9,     0,    52,     0,    53,    48,
      49,     0,     0,    50,    65,    66,     0,    51,     0,     0,
      24,    25,     0,    52,     0,    53,    70,    10,   117,   118,
     111,   112,   113,     0,   119,    24,    25,     0,     0,     0,
       0,   120,     0,   117,   118,     0,     0,     0,     0,   150,
       0,     0,     0,     0,     0,     0,   151
};

static const yytype_int16 yycheck[] =
{
      10,    61,    47,     0,     1,     3,    16,    17,    18,    19,
      20,    51,    22,    23,    45,    55,    26,    27,    28,    29,
      33,    31,     2,     3,    19,    20,    36,    40,    10,    24,
      30,    11,     3,    28,    32,     4,    45,    35,    19,    20,
      11,     6,     7,    24,     9,    24,    27,    44,    45,    31,
      27,    61,     6,     7,    34,     9,     1,    33,     3,     4,
       5,    26,    27,    39,   124,    30,    46,    32,    78,    34,
       6,     7,     8,     9,    41,    40,    33,    42,    32,    19,
      20,   126,   127,    40,    64,    45,    62,    45,    28,    69,
      26,    35,    72,    29,   104,    75,   106,   107,    69,    19,
      20,    72,     6,     7,    75,     9,    75,   117,   118,   154,
     120,    19,    20,   123,   124,    27,    96,    93,   100,   101,
     102,    20,   132,   168,   100,    96,    33,   172,    32,   139,
      33,    35,    36,    37,    25,   180,    19,    20,   148,    22,
     185,   151,   111,   112,   113,    19,    20,   157,    19,    20,
     160,   131,    25,    27,    25,    33,   166,    20,    21,   135,
     131,   137,   138,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    33,    18,    19,    20,    35,    36,
      37,    19,    20,    26,    27,    28,    29,    30,    45,    27,
      33,    34,    35,    36,    37,    38,    25,    40,    35,    36,
      43,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    25,    18,    19,    20,    33,    72,    19,    20,
      22,    26,    27,    28,    29,    30,    27,    32,    33,    34,
      19,    20,   158,    38,    34,    40,    25,    -1,    43,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      -1,    18,    19,    20,    19,    20,    -1,    22,    -1,    26,
      27,    28,    29,    30,    -1,    32,    33,    34,    -1,    19,
      20,    38,    -1,    40,    24,    -1,    43,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    -1,    18,
      19,    20,    19,    20,    -1,    22,    -1,    26,    27,    28,
      29,    30,    -1,    32,    33,    34,    -1,    19,    20,    38,
      -1,    40,    24,    -1,    43,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    -1,    18,    19,    20,
      19,    20,    -1,    22,    -1,    26,    27,    28,    29,    30,
      -1,    -1,    -1,    34,    35,    36,    37,    38,     6,     7,
       8,     9,    43,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    -1,    18,    19,    20,    26,    -1,
      -1,    29,    -1,    26,    27,    28,    29,    30,    -1,    19,
      20,    34,    35,    36,    24,    38,    -1,    -1,    -1,    -1,
      43,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    -1,    18,    19,    20,    -1,    -1,    -1,    -1,
      -1,    26,    27,    28,    29,    30,    19,    20,    -1,    34,
      35,    24,    -1,    38,    -1,    -1,    -1,    -1,    43,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      -1,    18,    19,    20,    -1,    -1,    -1,    -1,    -1,    26,
      27,    28,    29,    30,    19,    20,    -1,    34,    -1,    24,
      -1,    38,    -1,    -1,    -1,    -1,    43,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    -1,    18,
      19,    20,    -1,    -1,    -1,    -1,    -1,    26,    27,    28,
      29,    30,    -1,     6,     7,    34,     9,    19,    20,    38,
      19,    20,    24,    -1,    43,    24,    19,    20,    21,    22,
      23,    -1,    -1,    26,    27,    -1,    -1,    30,    31,    32,
      -1,    34,     6,     7,    -1,     9,    39,    40,    -1,    42,
      -1,    -1,    -1,    -1,    -1,    19,    20,    21,    22,    23,
      -1,    -1,    26,    27,    -1,    -1,    30,    31,    32,    -1,
      34,     6,     7,    -1,     9,    -1,    40,    -1,    42,    26,
      27,    -1,    -1,    30,    19,    20,    -1,    34,    -1,    -1,
      19,    20,    -1,    40,    -1,    42,    31,    32,    27,    28,
      35,    36,    37,    -1,    33,    19,    20,    -1,    -1,    -1,
      -1,    40,    -1,    27,    28,    -1,    -1,    -1,    -1,    33,
      -1,    -1,    -1,    -1,    -1,    -1,    40
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,    47,    45,     6,     7,     9,
      32,    49,    50,    51,    57,     8,    10,    11,    12,    13,
      14,    15,    16,    18,    19,    20,    26,    27,    28,    29,
      30,    32,    33,    34,    38,    40,    43,    48,    50,    52,
      53,    54,    55,    56,    57,    61,    63,    64,    26,    27,
      30,    34,    40,    42,    58,    59,    60,     0,     1,    44,
      45,    29,    52,    57,    63,    19,    20,    21,    22,    23,
      31,    50,    58,    62,    65,    66,    63,    63,    24,    63,
      63,    63,    63,    63,    63,    27,    63,    28,    63,    63,
      30,    32,    48,    52,    48,    53,    49,    48,    63,    45,
      35,    36,    37,    53,    24,     8,    26,    29,    57,    65,
      60,    35,    36,    37,    45,    60,    45,    27,    28,    33,
      40,    61,    33,    40,    29,    50,    20,    21,    50,    62,
      50,    58,    27,    63,    39,    52,    53,    52,    52,    27,
      63,    63,    63,    58,    58,    58,    63,    33,    40,    63,
      33,    40,    61,    65,    65,    50,    63,    27,    41,    63,
      27,    33,    63,    33,    63,    33,    40,    22,    65,    25,
      25,    63,    64,    25,    25,    63,    33,    33,    63,    22,
      65,    25,    25,    33,    22,    65,    22,    65,    22
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    46,    47,    47,    47,    47,    47,    47,    47,    47,
      47,    47,    47,    48,    48,    48,    48,    49,    49,    49,
      49,    49,    49,    49,    49,    49,    50,    50,    50,    50,
      51,    51,    51,    51,    51,    51,    51,    51,    52,    52,
      52,    52,    52,    52,    53,    53,    53,    54,    54,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    56,    56,    56,    56,
      56,    56,    56,    56,    57,    57,    57,    58,    58,    58,
      58,    59,    59,    60,    60,    60,    60,    60,    60,    61,
      61,    61,    61,    61,    61,    62,    62,    62,    62,    62,
      62,    62,    63,    63,    64,    64,    65,    65,    66
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     3,     2,     1,     3,     2,     1,     2,
       3,     2,     2,     2,     2,     2,     1,     1,     2,     3,
       3,     2,     3,     3,     4,     2,     1,     1,     3,     5,
       3,     5,     4,     6,     4,     6,     5,     7,     3,     3,
       3,     2,     3,     1,     2,     1,     1,     3,     5,     1,
       2,     1,     2,     2,     3,     1,     1,     1,     1,     1,
       1,     1,     2,     2,     2,     2,     2,     2,     2,     5,
       5,     6,     5,     5,     6,     2,     3,     2,     1,     1,
       2,     1,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     3,     3,     3,
       1,     1,     2,     1,     1,     1,     1,     1,     2,     2,
       2,     1,     2,     2,     1,     1,     3,     4,     5,     6,
       7,     8,     1,     1,     1,     2,     1,     1,     1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


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
        yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position, YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


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




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::vector<std::pair<unsigned int, unsigned int>>& branchPoints, void *scanner, int& start_token, unsigned int& current_token_position)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  YY_USE (input);
  YY_USE (molList);
  YY_USE (lastAtom);
  YY_USE (lastBond);
  YY_USE (numAtomsParsed);
  YY_USE (numBondsParsed);
  YY_USE (branchPoints);
  YY_USE (scanner);
  YY_USE (start_token);
  YY_USE (current_token_position);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::vector<std::pair<unsigned int, unsigned int>>& branchPoints, void *scanner, int& start_token, unsigned int& current_token_position)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position);
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
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::vector<std::pair<unsigned int, unsigned int>>& branchPoints, void *scanner, int& start_token, unsigned int& current_token_position)
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
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)], input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
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






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::vector<std::pair<unsigned int, unsigned int>>& branchPoints, void *scanner, int& start_token, unsigned int& current_token_position)
{
  YY_USE (yyvaluep);
  YY_USE (input);
  YY_USE (molList);
  YY_USE (lastAtom);
  YY_USE (lastBond);
  YY_USE (numAtomsParsed);
  YY_USE (numBondsParsed);
  YY_USE (branchPoints);
  YY_USE (scanner);
  YY_USE (start_token);
  YY_USE (current_token_position);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yykind)
    {
    case YYSYMBOL_ATOM_TOKEN: /* ATOM_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_SIMPLE_ATOM_QUERY_TOKEN: /* SIMPLE_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_COMPLEX_ATOM_QUERY_TOKEN: /* COMPLEX_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_MIN_RINGSIZE_ATOM_QUERY_TOKEN: /* MIN_RINGSIZE_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_RINGSIZE_ATOM_QUERY_TOKEN: /* RINGSIZE_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_RINGBOND_ATOM_QUERY_TOKEN: /* RINGBOND_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_IMPLICIT_H_ATOM_QUERY_TOKEN: /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_HYB_TOKEN: /* HYB_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_HETERONEIGHBOR_ATOM_QUERY_TOKEN: /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_ALIPHATIC: /* ALIPHATIC  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN: /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_BOND_TOKEN: /* BOND_TOKEN  */
            { delete ((*yyvaluep).bond); }
        break;

    case YYSYMBOL_atomd: /* atomd  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_hydrogen_atom: /* hydrogen_atom  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_atom_expr: /* atom_expr  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_point_query: /* point_query  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_recursive_query: /* recursive_query  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_atom_query: /* atom_query  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_possible_range_query: /* possible_range_query  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_simple_atom: /* simple_atom  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_bond_expr: /* bond_expr  */
            { delete ((*yyvaluep).bond); }
        break;

    case YYSYMBOL_bond_query: /* bond_query  */
            { delete ((*yyvaluep).bond); }
        break;

    case YYSYMBOL_bondd: /* bondd  */
            { delete ((*yyvaluep).bond); }
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
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::vector<std::pair<unsigned int, unsigned int>>& branchPoints, void *scanner, int& start_token, unsigned int& current_token_position)
{
/* Lookahead token kind.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

    /* Number of syntax errors so far.  */
    int yynerrs = 0;

    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

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
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
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
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
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

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex (&yylval, scanner, start_token, current_token_position);
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
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
  case 2: /* meta_start: START_MOL mol  */
              {
// the molList has already been updated, no need to do anything
}
    break;

  case 3: /* meta_start: START_ATOM atomd EOS_TOKEN  */
                             {
  lastAtom = (yyvsp[-1].atom);
  YYACCEPT;
}
    break;

  case 4: /* meta_start: START_ATOM bad_atom_def  */
                          {
  YYABORT;
}
    break;

  case 5: /* meta_start: START_ATOM  */
             {
  YYABORT;
}
    break;

  case 6: /* meta_start: START_BOND bond_expr EOS_TOKEN  */
                                 {
  lastBond = (yyvsp[-1].bond);
  YYACCEPT;
}
    break;

  case 7: /* meta_start: START_BOND bond_expr  */
                       {
  delete (yyvsp[0].bond);
  YYABORT;
}
    break;

  case 8: /* meta_start: START_BOND  */
             {
  YYABORT;
}
    break;

  case 9: /* meta_start: meta_start BAD_CHARACTER  */
                           {
  yyerrok;
  yyErrorCleanup(molList);
  yyerror(input, molList, current_token_position, "syntax error");
  YYABORT;
}
    break;

  case 10: /* meta_start: meta_start error EOS_TOKEN  */
                            {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
    break;

  case 11: /* meta_start: meta_start EOS_TOKEN  */
                       {
  YYACCEPT;
}
    break;

  case 12: /* meta_start: error EOS_TOKEN  */
                  {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
    break;

  case 16: /* bad_atom_def: atom_expr  */
            {
  delete (yyvsp[0].atom);
  YYABORT;
}
    break;

  case 17: /* mol: atomd  */
           {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1);
  (*molList)[ sz ]->addAtom((yyvsp[0].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
}
    break;

  case 18: /* mol: mol atomd  */
                  {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);

  QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(a1,mp->getAtomWithIdx(atomIdx2));
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond(newB,true);
}
    break;

  case 19: /* mol: mol bond_expr atomd  */
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
    break;

  case 20: /* mol: mol SEPARATOR_TOKEN atomd  */
                            {
  RWMol *mp = (*molList)[(yyval.moli)];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom((yyvsp[0].atom),true,true);
}
    break;

  case 21: /* mol: mol ring_number  */
                  {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();

  QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(atom, nullptr);
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
    break;

  case 22: /* mol: mol bond_expr ring_number  */
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
    break;

  case 23: /* mol: mol branch_open_token atomd  */
                              {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);

  QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(a1,mp->getAtomWithIdx(atomIdx2));
  newB->setOwningMol(mp);
  newB->setBeginAtomIdx(atomIdx1);
  newB->setEndAtomIdx(atomIdx2);
  newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond(newB);
  delete newB;

  branchPoints.push_back({atomIdx1, (yyvsp[-1].ival)});
}
    break;

  case 24: /* mol: mol branch_open_token bond_expr atomd  */
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
  branchPoints.push_back({atomIdx1, (yyvsp[-2].ival)});

}
    break;

  case 25: /* mol: mol GROUP_CLOSE_TOKEN  */
                        {
  if(branchPoints.empty()){
     yyerror(input,molList,branchPoints,scanner,start_token, current_token_position, "extra close parentheses");
     yyErrorCleanup(molList);
     YYABORT;
  }
  RWMol *mp = (*molList)[(yyval.moli)];
  mp->setActiveAtom(branchPoints.back().first);
  branchPoints.pop_back();
}
    break;

  case 28: /* atomd: ATOM_OPEN_TOKEN atom_expr ATOM_CLOSE_TOKEN  */
{
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 29: /* atomd: ATOM_OPEN_TOKEN atom_expr COLON_TOKEN number ATOM_CLOSE_TOKEN  */
{
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
    break;

  case 30: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN ATOM_CLOSE_TOKEN  */
{
  (yyval.atom) = new QueryAtom(1);
}
    break;

  case 31: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN  */
{
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
    break;

  case 32: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN ATOM_CLOSE_TOKEN  */
                                                  {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
    break;

  case 33: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN  */
                                                                     {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-4].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-4].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
    break;

  case 34: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN charge_spec ATOM_CLOSE_TOKEN  */
                                                       {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->getFlags() |= SMARTS_CHARGE_MASK;
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
    break;

  case 35: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN  */
                                                                          {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->getFlags() |= SMARTS_CHARGE_MASK;
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
    break;

  case 36: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN charge_spec ATOM_CLOSE_TOKEN  */
                                                              {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-3].ival));
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->getFlags() |= SMARTS_CHARGE_MASK;
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
    break;

  case 37: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN  */
                                                                                 {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-5].ival));
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->getFlags() |= SMARTS_CHARGE_MASK;
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-5].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
    break;

  case 38: /* atom_expr: atom_expr AND_TOKEN atom_expr  */
                                         {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if ((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) { (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag()); }
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
  (yyval.atom) = (yyvsp[-2].atom);
}
    break;

  case 39: /* atom_expr: atom_expr OR_TOKEN atom_expr  */
                               {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if ((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) { (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag()); }
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  (yyvsp[-2].atom)->setAtomicNum(0);
  delete (yyvsp[0].atom);
  (yyval.atom) = (yyvsp[-2].atom);
}
    break;

  case 40: /* atom_expr: atom_expr SEMI_TOKEN atom_expr  */
                                 {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if ((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) { (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag()); }
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
  (yyval.atom) = (yyvsp[-2].atom);
}
    break;

  case 41: /* atom_expr: atom_expr point_query  */
                        {
  atom_expr_and_point_query((yyvsp[-1].atom), (yyvsp[0].atom));
  delete (yyvsp[0].atom);
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 42: /* atom_expr: atom_expr AND_TOKEN point_query  */
                                  {
  atom_expr_and_point_query((yyvsp[-2].atom), (yyvsp[0].atom));
  delete (yyvsp[0].atom);
  (yyval.atom) = (yyvsp[-2].atom);
}
    break;

  case 44: /* point_query: NOT_TOKEN point_query  */
                                   {
  (yyvsp[0].atom)->getQuery()->setNegation(!((yyvsp[0].atom)->getQuery()->getNegation()));
  (yyvsp[0].atom)->setAtomicNum(0);
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[0].atom));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 47: /* recursive_query: BEGIN_RECURSE mol END_RECURSE  */
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
    break;

  case 48: /* recursive_query: BEGIN_RECURSE mol END_RECURSE UNDERSCORE_TOKEN nonzero_number  */
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
    break;

  case 50: /* atom_query: number simple_atom  */
                     {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
    break;

  case 52: /* atom_query: number ATOM_TOKEN  */
                    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
    break;

  case 53: /* atom_query: HASH_TOKEN number  */
                    { (yyval.atom) = new QueryAtom((yyvsp[0].ival)); }
    break;

  case 54: /* atom_query: number HASH_TOKEN number  */
                           {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsotope((yyvsp[-2].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
}
    break;

  case 62: /* atom_query: COMPLEX_ATOM_QUERY_TOKEN number  */
                                  {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 63: /* atom_query: HETERONEIGHBOR_ATOM_QUERY_TOKEN number  */
                                         {
  (yyvsp[-1].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 64: /* atom_query: ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN number  */
                                                  {
  (yyvsp[-1].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 65: /* atom_query: MIN_RINGSIZE_ATOM_QUERY_TOKEN number  */
                                       {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 66: /* atom_query: RINGSIZE_ATOM_QUERY_TOKEN number  */
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomInRingOfSizeQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 67: /* atom_query: RINGBOND_ATOM_QUERY_TOKEN number  */
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 68: /* atom_query: IMPLICIT_H_ATOM_QUERY_TOKEN number  */
                                     {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 69: /* atom_query: possible_range_query RANGE_OPEN_TOKEN MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>((yyvsp[-1].ival),oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-4].atom);
}
    break;

  case 70: /* atom_query: possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN RANGE_CLOSE_TOKEN  */
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>((yyvsp[-2].ival),oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-4].atom);
}
    break;

  case 71: /* atom_query: possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
                                                                                    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-5].atom);
}
    break;

  case 72: /* atom_query: RINGSIZE_ATOM_QUERY_TOKEN RANGE_OPEN_TOKEN MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
                                                                                  {
  int lv = -1;
  int uv = (yyvsp[-1].ival);
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>(uv,[lv,uv](Atom const *at) {
            return queryAtomIsInRingOfSize(at, lv, uv);
          },std::string("greater_AtomRingSize"));
  (yyvsp[-4].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-4].atom);
}
    break;

  case 73: /* atom_query: RINGSIZE_ATOM_QUERY_TOKEN RANGE_OPEN_TOKEN number MINUS_TOKEN RANGE_CLOSE_TOKEN  */
                                                                                  {
  int lv = (yyvsp[-2].ival);
  int uv = -1;
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(lv,[lv,uv](Atom const *at) {
            return queryAtomIsInRingOfSize(at, lv, uv);
          },std::string("less_AtomRingSize"));
  (yyvsp[-4].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-4].atom);
}
    break;

  case 74: /* atom_query: RINGSIZE_ATOM_QUERY_TOKEN RANGE_OPEN_TOKEN number MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
                                                                                         {
  int lv = (yyvsp[-3].ival);
  int uv = (yyvsp[-1].ival);
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery(lv,uv,false,false,[lv,uv](Atom const *at) {
            return queryAtomIsInRingOfSize(at, lv, uv);
          },std::string("range_AtomRingSize"));
  (yyvsp[-5].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-5].atom);
}
    break;

  case 75: /* atom_query: number H_TOKEN  */
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)));
  newQ->setIsotope((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs(1);
  newQ->setNoImplicit(true);
  newQ->getFlags() |= SMARTS_H_MASK;
  (yyval.atom)=newQ;
}
    break;

  case 76: /* atom_query: number H_TOKEN number  */
                        {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)));
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomHCountQuery((yyvsp[0].ival)),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs((yyvsp[0].ival));
  newQ->setNoImplicit(true);
  newQ->getFlags() |= SMARTS_H_MASK;
  (yyval.atom)=newQ;
}
    break;

  case 77: /* atom_query: H_TOKEN number  */
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  newQ->setNumExplicitHs((yyvsp[0].ival));
  newQ->setNoImplicit(true);
  newQ->getFlags() |= SMARTS_H_MASK;
  (yyval.atom)=newQ;
  
}
    break;

  case 78: /* atom_query: H_TOKEN  */
          {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  newQ->setNoImplicit(true);
  newQ->getFlags() |= SMARTS_H_MASK;
  (yyval.atom)=newQ;
}
    break;

  case 79: /* atom_query: charge_spec  */
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  newQ->setFormalCharge((yyvsp[0].ival));
  newQ->getFlags() |= SMARTS_CHARGE_MASK;
  (yyval.atom)=newQ;
}
    break;

  case 80: /* atom_query: AT_TOKEN AT_TOKEN  */
                    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
    break;

  case 81: /* atom_query: AT_TOKEN  */
           {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
    break;

  case 82: /* atom_query: CHI_CLASS_TOKEN  */
                  {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[0].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,0);
  (yyval.atom)=newQ;
}
    break;

  case 83: /* atom_query: CHI_CLASS_TOKEN number  */
                         {
  if((yyvsp[0].ival)==0){
    yyerror(input,molList,branchPoints,scanner,start_token, current_token_position,
            "chiral permutation cannot be zero");
    yyErrorCleanup(molList);
    YYABORT;
  }

  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[-1].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,(yyvsp[0].ival));
  (yyval.atom)=newQ;
}
    break;

  case 85: /* atom_query: number  */
         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
    break;

  case 87: /* possible_range_query: HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
                                  {
  (yyvsp[0].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 88: /* possible_range_query: ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
                                           {
  (yyvsp[0].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 89: /* possible_range_query: MIN_RINGSIZE_ATOM_QUERY_TOKEN  */
                                {
  (yyvsp[0].atom)->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 90: /* possible_range_query: RINGBOND_ATOM_QUERY_TOKEN  */
                            {
  (yyvsp[0].atom)->setQuery(makeAtomRingBondCountQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 91: /* possible_range_query: IMPLICIT_H_ATOM_QUERY_TOKEN  */
                              {
  (yyvsp[0].atom)->setQuery(makeAtomImplicitHCountQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 92: /* possible_range_query: PLUS_TOKEN  */
             {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
    break;

  case 93: /* possible_range_query: MINUS_TOKEN  */
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNegativeFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
    break;

  case 94: /* simple_atom: ORGANIC_ATOM_TOKEN  */
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
    break;

  case 95: /* simple_atom: AROMATIC_ATOM_TOKEN  */
                      {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->setQuery(makeAtomTypeQuery((yyvsp[0].ival),true));
}
    break;

  case 97: /* bond_expr: bond_expr AND_TOKEN bond_expr  */
                                        {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-2].bond);
}
    break;

  case 98: /* bond_expr: bond_expr OR_TOKEN bond_expr  */
                               {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-2].bond);
}
    break;

  case 99: /* bond_expr: bond_expr SEMI_TOKEN bond_expr  */
                                 {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-2].bond);
}
    break;

  case 102: /* bond_query: bond_query bondd  */
                   {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-1].bond);
}
    break;

  case 104: /* bondd: MINUS_TOKEN  */
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
    break;

  case 105: /* bondd: HASH_TOKEN  */
             {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
    break;

  case 106: /* bondd: COLON_TOKEN  */
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
    break;

  case 107: /* bondd: AT_TOKEN  */
           {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
    break;

  case 108: /* bondd: NOT_TOKEN bondd  */
                  {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
    break;

  case 109: /* charge_spec: PLUS_TOKEN PLUS_TOKEN  */
                                   { (yyval.ival)=2; }
    break;

  case 110: /* charge_spec: PLUS_TOKEN number  */
                    { (yyval.ival)=(yyvsp[0].ival); }
    break;

  case 111: /* charge_spec: PLUS_TOKEN  */
             { (yyval.ival)=1; }
    break;

  case 112: /* charge_spec: MINUS_TOKEN MINUS_TOKEN  */
                          { (yyval.ival)=-2; }
    break;

  case 113: /* charge_spec: MINUS_TOKEN number  */
                     { (yyval.ival)=-(yyvsp[0].ival); }
    break;

  case 114: /* charge_spec: MINUS_TOKEN  */
              { (yyval.ival)=-1; }
    break;

  case 116: /* ring_number: PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit  */
                                          { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
    break;

  case 117: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN  */
                                                         { (yyval.ival) = (yyvsp[-1].ival); }
    break;

  case 118: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN  */
                                                               { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 119: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN  */
                                                                     { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 120: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN  */
                                                                           { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 121: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN  */
                                                                                 { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 125: /* nonzero_number: nonzero_number digit  */
                       {
    if((yyvsp[-1].ival) >= std::numeric_limits<std::int32_t>::max()/10 ||
     (yyvsp[-1].ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yyvsp[0].ival) ){
     yysmarts_error(input,molList,lastAtom,lastBond,numAtomsParsed,numBondsParsed,branchPoints,scanner,start_token, current_token_position, "number too large");
     YYABORT;
  }
  (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
    break;

  case 128: /* branch_open_token: GROUP_OPEN_TOKEN  */
                                    { (yyval.ival) = current_token_position; }
    break;



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
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

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
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position, YY_("syntax error"));
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
                      yytoken, &yylval, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position);
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
  ++yynerrs;

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

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
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
                  YY_ACCESSING_SYMBOL (yystate), yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position, YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, current_token_position);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}


