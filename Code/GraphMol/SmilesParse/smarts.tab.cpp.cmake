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
  YYSYMBOL_RINGSIZE_ATOM_QUERY_TOKEN = 11, /* RINGSIZE_ATOM_QUERY_TOKEN  */
  YYSYMBOL_RINGBOND_ATOM_QUERY_TOKEN = 12, /* RINGBOND_ATOM_QUERY_TOKEN  */
  YYSYMBOL_IMPLICIT_H_ATOM_QUERY_TOKEN = 13, /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
  YYSYMBOL_HYB_TOKEN = 14,                 /* HYB_TOKEN  */
  YYSYMBOL_HETERONEIGHBOR_ATOM_QUERY_TOKEN = 15, /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
  YYSYMBOL_ALIPHATIC = 16,                 /* ALIPHATIC  */
  YYSYMBOL_ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN = 17, /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
  YYSYMBOL_ZERO_TOKEN = 18,                /* ZERO_TOKEN  */
  YYSYMBOL_NONZERO_DIGIT_TOKEN = 19,       /* NONZERO_DIGIT_TOKEN  */
  YYSYMBOL_GROUP_OPEN_TOKEN = 20,          /* GROUP_OPEN_TOKEN  */
  YYSYMBOL_GROUP_CLOSE_TOKEN = 21,         /* GROUP_CLOSE_TOKEN  */
  YYSYMBOL_SEPARATOR_TOKEN = 22,           /* SEPARATOR_TOKEN  */
  YYSYMBOL_RANGE_OPEN_TOKEN = 23,          /* RANGE_OPEN_TOKEN  */
  YYSYMBOL_RANGE_CLOSE_TOKEN = 24,         /* RANGE_CLOSE_TOKEN  */
  YYSYMBOL_HASH_TOKEN = 25,                /* HASH_TOKEN  */
  YYSYMBOL_MINUS_TOKEN = 26,               /* MINUS_TOKEN  */
  YYSYMBOL_PLUS_TOKEN = 27,                /* PLUS_TOKEN  */
  YYSYMBOL_H_TOKEN = 28,                   /* H_TOKEN  */
  YYSYMBOL_AT_TOKEN = 29,                  /* AT_TOKEN  */
  YYSYMBOL_PERCENT_TOKEN = 30,             /* PERCENT_TOKEN  */
  YYSYMBOL_ATOM_OPEN_TOKEN = 31,           /* ATOM_OPEN_TOKEN  */
  YYSYMBOL_ATOM_CLOSE_TOKEN = 32,          /* ATOM_CLOSE_TOKEN  */
  YYSYMBOL_NOT_TOKEN = 33,                 /* NOT_TOKEN  */
  YYSYMBOL_AND_TOKEN = 34,                 /* AND_TOKEN  */
  YYSYMBOL_OR_TOKEN = 35,                  /* OR_TOKEN  */
  YYSYMBOL_SEMI_TOKEN = 36,                /* SEMI_TOKEN  */
  YYSYMBOL_BEGIN_RECURSE = 37,             /* BEGIN_RECURSE  */
  YYSYMBOL_END_RECURSE = 38,               /* END_RECURSE  */
  YYSYMBOL_COLON_TOKEN = 39,               /* COLON_TOKEN  */
  YYSYMBOL_UNDERSCORE_TOKEN = 40,          /* UNDERSCORE_TOKEN  */
  YYSYMBOL_BOND_TOKEN = 41,                /* BOND_TOKEN  */
  YYSYMBOL_CHI_CLASS_TOKEN = 42,           /* CHI_CLASS_TOKEN  */
  YYSYMBOL_BAD_CHARACTER = 43,             /* BAD_CHARACTER  */
  YYSYMBOL_EOS_TOKEN = 44,                 /* EOS_TOKEN  */
  YYSYMBOL_YYACCEPT = 45,                  /* $accept  */
  YYSYMBOL_meta_start = 46,                /* meta_start  */
  YYSYMBOL_bad_atom_def = 47,              /* bad_atom_def  */
  YYSYMBOL_mol = 48,                       /* mol  */
  YYSYMBOL_atomd = 49,                     /* atomd  */
  YYSYMBOL_hydrogen_atom = 50,             /* hydrogen_atom  */
  YYSYMBOL_atom_expr = 51,                 /* atom_expr  */
  YYSYMBOL_point_query = 52,               /* point_query  */
  YYSYMBOL_recursive_query = 53,           /* recursive_query  */
  YYSYMBOL_atom_query = 54,                /* atom_query  */
  YYSYMBOL_possible_range_query = 55,      /* possible_range_query  */
  YYSYMBOL_simple_atom = 56,               /* simple_atom  */
  YYSYMBOL_bond_expr = 57,                 /* bond_expr  */
  YYSYMBOL_bond_query = 58,                /* bond_query  */
  YYSYMBOL_bondd = 59,                     /* bondd  */
  YYSYMBOL_charge_spec = 60,               /* charge_spec  */
  YYSYMBOL_ring_number = 61,               /* ring_number  */
  YYSYMBOL_number = 62,                    /* number  */
  YYSYMBOL_nonzero_number = 63,            /* nonzero_number  */
  YYSYMBOL_digit = 64,                     /* digit  */
  YYSYMBOL_branch_open_token = 65          /* branch_open_token  */
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
#define YYFINAL  56
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   598

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  45
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  21
/* YYNRULES -- Number of rules.  */
#define YYNRULES  123
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  178

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   299


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
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   187,   187,   190,   194,   197,   200,   204,   208,   211,
     217,   222,   225,   233,   234,   235,   236,   244,   253,   267,
     288,   294,   318,   339,   356,   379,   393,   394,   395,   399,
     422,   426,   431,   437,   446,   453,   462,   471,   485,   492,
     500,   507,   512,   517,   520,   526,   527,   531,   548,   572,
     573,   578,   579,   584,   585,   590,   591,   592,   593,   594,
     595,   596,   600,   604,   608,   612,   616,   620,   627,   634,
     642,   652,   662,   671,   679,   686,   692,   698,   705,   712,
     713,   720,   721,   725,   729,   733,   737,   741,   746,   754,
     766,   771,   776,   781,   786,   791,   794,   795,   803,   804,
     810,   816,   822,   827,   834,   835,   836,   837,   838,   839,
     843,   844,   845,   846,   847,   848,   849,   854,   855,   859,
     860,   869,   870,   874
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
  "RINGSIZE_ATOM_QUERY_TOKEN", "RINGBOND_ATOM_QUERY_TOKEN",
  "IMPLICIT_H_ATOM_QUERY_TOKEN", "HYB_TOKEN",
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

#define YYPACT_NINF (-59)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-89)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     171,   -16,    30,   189,    87,    23,   -59,   -59,   -59,   -59,
     411,   501,   -59,   -59,   -59,   -59,    26,   243,   280,   350,
     -59,   386,   423,   -59,   -59,    39,   122,    41,    39,     5,
     226,   263,   448,    30,   263,    39,   -59,    -4,   300,   -59,
     -59,   -59,    25,    28,   -59,    46,    64,   -59,   -59,   -59,
      87,   -59,   -59,     7,    87,   -59,   -59,    31,   -59,   -59,
     557,   152,   -59,   332,   -59,   -59,   -59,   -59,    30,   163,
     -59,   116,   -59,   -59,   529,   -59,   -59,   -59,   -59,   -59,
     -59,   -59,   -59,   -59,   -59,   -59,   -59,   -59,   263,   -59,
     152,   -59,   -59,   473,   -59,   -59,   -59,   448,   448,   448,
     -59,    70,   -59,    39,    39,   -59,   -59,   -59,    87,    87,
      87,   -59,   -59,   -59,    82,    59,   -59,    39,   -18,   -59,
      39,   559,   -59,    64,    64,   -59,   -59,   -59,   538,     6,
     448,   -59,   374,   337,    39,    43,   -59,   -59,   -59,    45,
     158,    60,   -59,    39,    83,   -59,    39,    -1,   -59,   191,
     -59,    98,    78,   205,   -59,    97,   -59,   100,   -59,    39,
     -59,   228,    64,   -59,   -59,   118,   -59,   -59,   121,   -59,
     265,   -59,   -59,   -59,   302,   -59,   134,   -59
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     0,     5,     8,     0,    12,    90,    89,    91,
       0,     2,    17,    27,    26,    51,    55,    58,    59,    60,
      79,    56,    57,   117,   119,     0,   109,   106,    73,    76,
       0,     0,     0,     0,     0,    77,     4,     0,    16,    43,
      45,    46,     0,    49,    74,    80,   118,   100,    99,   102,
       0,   101,    98,     7,    95,    96,     1,     0,     9,    11,
      73,     0,    49,    80,   122,   121,   123,    25,     0,     0,
      18,     0,    21,   110,     0,    61,    64,    65,    66,    62,
      63,    53,   107,   108,   104,   105,    72,    75,     0,    13,
      16,    14,    44,     0,    15,    78,     3,     0,     0,     0,
      41,     0,    52,     0,    70,    50,   120,   103,     0,     0,
       0,     6,    97,    10,   109,   106,    30,     0,     0,    28,
       0,    70,    20,     0,     0,    19,    22,    23,     0,    47,
      38,    42,    39,    40,     0,     0,    54,    71,    92,    93,
      94,     0,    34,     0,     0,    32,     0,     0,   111,     0,
      24,     0,     0,     0,    31,     0,    29,     0,    36,     0,
     112,     0,    48,    67,    68,     0,    35,    33,     0,   113,
       0,    69,    37,   114,     0,   115,     0,   116
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -59,   -59,    -2,   135,    16,   -59,     0,    24,   -59,   -59,
     -59,     2,    29,   -59,   -28,   -58,   102,   -10,    54,   -43,
     -59
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     5,    89,    11,    12,    13,    38,    39,    40,    41,
      42,    62,    71,    54,    55,    44,    72,    45,    46,    73,
      74
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      63,    36,   118,   106,    14,    43,    75,    76,    77,    78,
      61,    79,    80,    14,   142,    81,    83,    85,    86,    37,
      63,   143,   107,    56,    57,    95,   112,    70,     6,    91,
      90,   158,    94,    53,    87,    14,     7,     8,   159,     9,
      96,   108,   109,   110,    23,    24,   151,   105,   101,   -81,
      86,   111,     7,     8,   102,     9,    92,    23,    24,    23,
      24,    10,   100,   147,   -87,   105,    58,    59,    84,   153,
      14,   103,   -26,    14,   104,   113,    14,    23,    24,   108,
     148,   149,    64,    65,   122,   100,    84,   125,    23,    24,
     127,   135,   154,   136,   137,    14,   134,   130,   132,   133,
      23,    24,   163,   128,    83,    85,   161,   141,    82,    70,
     144,   137,    47,    48,   100,   156,    49,    24,   170,   106,
      50,   131,     7,     8,   152,     9,    51,   174,    52,   166,
      14,   176,   167,   155,    64,    65,   157,   138,   139,   140,
      23,    24,   171,   165,   150,   -88,    69,    10,    82,   168,
     108,   109,   110,   172,   100,   177,   100,   100,     7,     8,
      15,     9,    16,    17,    18,    19,    20,    21,    93,    22,
      23,    24,     1,   126,     2,     3,     4,    25,    26,    27,
      28,    29,   123,   124,   119,    32,    97,    98,    99,    33,
       0,   120,   108,   109,    35,     7,     8,    15,     9,    16,
      17,    18,    19,    20,    21,   162,    22,    23,    24,    64,
      65,     0,   160,     0,    25,    26,    27,    28,    29,     0,
      30,    31,    32,    23,    24,     0,    33,     0,    34,   164,
       0,    35,     7,     8,    15,     9,    16,    17,    18,    19,
      20,    21,     0,    22,    23,    24,    64,    65,     0,   169,
       0,    25,    26,    27,    60,    29,     0,    88,    31,    32,
       0,    23,    24,    33,     0,    34,   -84,     0,    35,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,     0,
      22,    23,    24,    64,    65,     0,   173,     0,    25,    26,
      27,    28,    29,     0,    88,    31,    32,     0,    23,    24,
      33,     0,    34,   -85,     0,    35,     7,     8,    15,     9,
      16,    17,    18,    19,    20,    21,     0,    22,    23,    24,
      64,    65,     0,   175,     0,    25,    26,    27,    28,    29,
       0,     0,     0,    32,    97,    98,    99,    33,     7,     8,
     102,     9,    35,     7,     8,    15,     9,    16,    17,    18,
      19,    20,    21,     0,    22,    23,    24,   103,     0,     0,
     121,     0,    25,    26,    27,    28,    29,     0,    23,    24,
      32,    97,    98,   -86,    33,     0,     0,     0,     0,    35,
       7,     8,    15,     9,    16,    17,    18,    19,    20,    21,
       0,    22,    23,    24,     0,     0,     0,     0,     0,    25,
      26,    27,    28,    29,    23,    24,     0,    32,    97,   -82,
       0,    33,     0,     0,     0,     0,    35,     7,     8,    15,
       9,    16,    17,    18,    19,    20,    21,     0,    22,    23,
      24,     0,     0,     0,     0,     0,    25,    26,    27,    60,
      29,    23,    24,     0,    32,     0,   -83,     0,    33,     0,
       0,     0,     0,    35,     7,     8,    15,     9,    16,    17,
      18,    19,    20,    21,     0,    22,    23,    24,     0,     0,
       0,     0,     0,    25,    26,    27,    28,    29,     0,     7,
       8,    32,     9,     0,     0,    33,     0,     0,     0,     0,
      35,    64,    65,    66,    67,    68,     0,     0,    47,    48,
       0,     0,    49,    69,    10,     0,    50,     7,     8,     0,
       9,   129,    51,     0,    52,     0,     0,     0,     0,    64,
      65,    66,    67,    68,     0,     0,    47,    48,     0,     0,
      49,    69,    10,     0,    50,     7,     8,     0,     9,     0,
      51,     0,    52,     0,     7,     8,     0,     9,     0,     0,
       0,     0,     0,     0,    47,    48,     0,     0,    49,     0,
      10,     0,    50,     0,     0,     0,     0,     0,    51,    10,
      52,     0,   108,   109,   110,    23,    24,    23,    24,     0,
       0,     0,     0,   114,   115,   114,   115,     0,     0,   116,
       0,   145,     0,     0,     0,     0,   117,     0,   146
};

static const yytype_int16 yycheck[] =
{
      10,     3,    60,    46,     2,     3,    16,    17,    18,    19,
      10,    21,    22,    11,    32,    25,    26,    27,    28,     3,
      30,    39,    50,     0,     1,    35,    54,    11,    44,    31,
      30,    32,    34,     4,    29,    33,     6,     7,    39,     9,
      44,    34,    35,    36,    18,    19,    40,    45,    23,    23,
      60,    44,     6,     7,     8,     9,    32,    18,    19,    18,
      19,    31,    38,   121,    23,    63,    43,    44,    27,    26,
      68,    25,    44,    71,    28,    44,    74,    18,    19,    34,
     123,   124,    18,    19,    68,    61,    27,    71,    18,    19,
      74,   101,    32,   103,   104,    93,    26,    97,    98,    99,
      18,    19,    24,    74,   114,   115,   149,   117,    26,    93,
     120,   121,    25,    26,    90,    32,    29,    19,   161,   162,
      33,    97,     6,     7,   134,     9,    39,   170,    41,    32,
     128,   174,    32,   143,    18,    19,   146,   108,   109,   110,
      18,    19,    24,   153,   128,    23,    30,    31,    26,   159,
      34,    35,    36,    32,   130,    21,   132,   133,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    33,    17,
      18,    19,     1,    71,     3,     4,     5,    25,    26,    27,
      28,    29,    19,    20,    32,    33,    34,    35,    36,    37,
      -1,    39,    34,    35,    42,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,   151,    17,    18,    19,    18,
      19,    -1,    21,    -1,    25,    26,    27,    28,    29,    -1,
      31,    32,    33,    18,    19,    -1,    37,    -1,    39,    24,
      -1,    42,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    -1,    17,    18,    19,    18,    19,    -1,    21,
      -1,    25,    26,    27,    28,    29,    -1,    31,    32,    33,
      -1,    18,    19,    37,    -1,    39,    23,    -1,    42,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    -1,
      17,    18,    19,    18,    19,    -1,    21,    -1,    25,    26,
      27,    28,    29,    -1,    31,    32,    33,    -1,    18,    19,
      37,    -1,    39,    23,    -1,    42,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    -1,    17,    18,    19,
      18,    19,    -1,    21,    -1,    25,    26,    27,    28,    29,
      -1,    -1,    -1,    33,    34,    35,    36,    37,     6,     7,
       8,     9,    42,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    -1,    17,    18,    19,    25,    -1,    -1,
      28,    -1,    25,    26,    27,    28,    29,    -1,    18,    19,
      33,    34,    35,    23,    37,    -1,    -1,    -1,    -1,    42,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      -1,    17,    18,    19,    -1,    -1,    -1,    -1,    -1,    25,
      26,    27,    28,    29,    18,    19,    -1,    33,    34,    23,
      -1,    37,    -1,    -1,    -1,    -1,    42,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    -1,    17,    18,
      19,    -1,    -1,    -1,    -1,    -1,    25,    26,    27,    28,
      29,    18,    19,    -1,    33,    -1,    23,    -1,    37,    -1,
      -1,    -1,    -1,    42,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    -1,    17,    18,    19,    -1,    -1,
      -1,    -1,    -1,    25,    26,    27,    28,    29,    -1,     6,
       7,    33,     9,    -1,    -1,    37,    -1,    -1,    -1,    -1,
      42,    18,    19,    20,    21,    22,    -1,    -1,    25,    26,
      -1,    -1,    29,    30,    31,    -1,    33,     6,     7,    -1,
       9,    38,    39,    -1,    41,    -1,    -1,    -1,    -1,    18,
      19,    20,    21,    22,    -1,    -1,    25,    26,    -1,    -1,
      29,    30,    31,    -1,    33,     6,     7,    -1,     9,    -1,
      39,    -1,    41,    -1,     6,     7,    -1,     9,    -1,    -1,
      -1,    -1,    -1,    -1,    25,    26,    -1,    -1,    29,    -1,
      31,    -1,    33,    -1,    -1,    -1,    -1,    -1,    39,    31,
      41,    -1,    34,    35,    36,    18,    19,    18,    19,    -1,
      -1,    -1,    -1,    26,    27,    26,    27,    -1,    -1,    32,
      -1,    32,    -1,    -1,    -1,    -1,    39,    -1,    39
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,    46,    44,     6,     7,     9,
      31,    48,    49,    50,    56,     8,    10,    11,    12,    13,
      14,    15,    17,    18,    19,    25,    26,    27,    28,    29,
      31,    32,    33,    37,    39,    42,    47,    49,    51,    52,
      53,    54,    55,    56,    60,    62,    63,    25,    26,    29,
      33,    39,    41,    57,    58,    59,     0,     1,    43,    44,
      28,    51,    56,    62,    18,    19,    20,    21,    22,    30,
      49,    57,    61,    64,    65,    62,    62,    62,    62,    62,
      62,    62,    26,    62,    27,    62,    62,    29,    31,    47,
      51,    47,    52,    48,    47,    62,    44,    34,    35,    36,
      52,    23,     8,    25,    28,    56,    64,    59,    34,    35,
      36,    44,    59,    44,    26,    27,    32,    39,    60,    32,
      39,    28,    49,    19,    20,    49,    61,    49,    57,    38,
      51,    52,    51,    51,    26,    62,    62,    62,    57,    57,
      57,    62,    32,    39,    62,    32,    39,    60,    64,    64,
      49,    40,    62,    26,    32,    62,    32,    62,    32,    39,
      21,    64,    63,    24,    24,    62,    32,    32,    62,    21,
      64,    24,    32,    21,    64,    21,    64,    21
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    45,    46,    46,    46,    46,    46,    46,    46,    46,
      46,    46,    46,    47,    47,    47,    47,    48,    48,    48,
      48,    48,    48,    48,    48,    48,    49,    49,    49,    49,
      50,    50,    50,    50,    50,    50,    50,    50,    51,    51,
      51,    51,    51,    51,    52,    52,    52,    53,    53,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    55,    55,    55,    55,    55,    55,    55,    55,    56,
      56,    56,    57,    57,    57,    57,    58,    58,    59,    59,
      59,    59,    59,    59,    60,    60,    60,    60,    60,    60,
      61,    61,    61,    61,    61,    61,    61,    62,    62,    63,
      63,    64,    64,    65
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
       1,     2,     2,     2,     2,     2,     2,     5,     5,     6,
       2,     3,     2,     1,     1,     2,     1,     1,     2,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     3,     3,     3,     1,     1,     2,     1,     1,
       1,     1,     1,     2,     2,     2,     1,     2,     2,     1,
       1,     3,     4,     5,     6,     7,     8,     1,     1,     1,
       2,     1,     1,     1
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
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
  (yyval.atom) = (yyvsp[-2].atom);
}
    break;

  case 39: /* atom_expr: atom_expr OR_TOKEN atom_expr  */
                               {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  (yyvsp[-2].atom)->setAtomicNum(0);
  delete (yyvsp[0].atom);
  (yyval.atom) = (yyvsp[-2].atom);
}
    break;

  case 40: /* atom_expr: atom_expr SEMI_TOKEN atom_expr  */
                                 {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
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

  case 61: /* atom_query: COMPLEX_ATOM_QUERY_TOKEN number  */
                                  {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 62: /* atom_query: HETERONEIGHBOR_ATOM_QUERY_TOKEN number  */
                                         {
  (yyvsp[-1].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 63: /* atom_query: ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN number  */
                                                  {
  (yyvsp[-1].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 64: /* atom_query: RINGSIZE_ATOM_QUERY_TOKEN number  */
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 65: /* atom_query: RINGBOND_ATOM_QUERY_TOKEN number  */
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 66: /* atom_query: IMPLICIT_H_ATOM_QUERY_TOKEN number  */
                                     {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
  (yyval.atom) = (yyvsp[-1].atom);
}
    break;

  case 67: /* atom_query: possible_range_query RANGE_OPEN_TOKEN MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>((yyvsp[-1].ival),oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-4].atom);
}
    break;

  case 68: /* atom_query: possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN RANGE_CLOSE_TOKEN  */
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>((yyvsp[-2].ival),oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-4].atom);
}
    break;

  case 69: /* atom_query: possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
                                                                                    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
  (yyval.atom) = (yyvsp[-5].atom);
}
    break;

  case 70: /* atom_query: number H_TOKEN  */
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

  case 71: /* atom_query: number H_TOKEN number  */
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

  case 72: /* atom_query: H_TOKEN number  */
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  newQ->setNumExplicitHs((yyvsp[0].ival));
  newQ->setNoImplicit(true);
  newQ->getFlags() |= SMARTS_H_MASK;
  (yyval.atom)=newQ;
  
}
    break;

  case 73: /* atom_query: H_TOKEN  */
          {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  newQ->setNoImplicit(true);
  newQ->getFlags() |= SMARTS_H_MASK;
  (yyval.atom)=newQ;
}
    break;

  case 74: /* atom_query: charge_spec  */
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  newQ->setFormalCharge((yyvsp[0].ival));
  newQ->getFlags() |= SMARTS_CHARGE_MASK;
  (yyval.atom)=newQ;
}
    break;

  case 75: /* atom_query: AT_TOKEN AT_TOKEN  */
                    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
    break;

  case 76: /* atom_query: AT_TOKEN  */
           {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
    break;

  case 77: /* atom_query: CHI_CLASS_TOKEN  */
                  {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[0].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,0);
  (yyval.atom)=newQ;
}
    break;

  case 78: /* atom_query: CHI_CLASS_TOKEN number  */
                         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[-1].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,(yyvsp[0].ival));
  (yyval.atom)=newQ;
}
    break;

  case 80: /* atom_query: number  */
         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
    break;

  case 82: /* possible_range_query: HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
                                  {
  (yyvsp[0].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 83: /* possible_range_query: ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
                                           {
  (yyvsp[0].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 84: /* possible_range_query: RINGSIZE_ATOM_QUERY_TOKEN  */
                            {
  (yyvsp[0].atom)->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 85: /* possible_range_query: RINGBOND_ATOM_QUERY_TOKEN  */
                            {
  (yyvsp[0].atom)->setQuery(makeAtomRingBondCountQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 86: /* possible_range_query: IMPLICIT_H_ATOM_QUERY_TOKEN  */
                              {
  (yyvsp[0].atom)->setQuery(makeAtomImplicitHCountQuery(0));
  (yyval.atom) = (yyvsp[0].atom);
}
    break;

  case 87: /* possible_range_query: PLUS_TOKEN  */
             {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
    break;

  case 88: /* possible_range_query: MINUS_TOKEN  */
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNegativeFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
    break;

  case 89: /* simple_atom: ORGANIC_ATOM_TOKEN  */
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

  case 90: /* simple_atom: AROMATIC_ATOM_TOKEN  */
                      {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->setQuery(makeAtomTypeQuery((yyvsp[0].ival),true));
}
    break;

  case 92: /* bond_expr: bond_expr AND_TOKEN bond_expr  */
                                        {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-2].bond);
}
    break;

  case 93: /* bond_expr: bond_expr OR_TOKEN bond_expr  */
                               {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-2].bond);
}
    break;

  case 94: /* bond_expr: bond_expr SEMI_TOKEN bond_expr  */
                                 {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-2].bond);
}
    break;

  case 97: /* bond_query: bond_query bondd  */
                   {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
  (yyval.bond) = (yyvsp[-1].bond);
}
    break;

  case 99: /* bondd: MINUS_TOKEN  */
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
    break;

  case 100: /* bondd: HASH_TOKEN  */
             {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
    break;

  case 101: /* bondd: COLON_TOKEN  */
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
    break;

  case 102: /* bondd: AT_TOKEN  */
           {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
    break;

  case 103: /* bondd: NOT_TOKEN bondd  */
                  {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
    break;

  case 104: /* charge_spec: PLUS_TOKEN PLUS_TOKEN  */
                                   { (yyval.ival)=2; }
    break;

  case 105: /* charge_spec: PLUS_TOKEN number  */
                    { (yyval.ival)=(yyvsp[0].ival); }
    break;

  case 106: /* charge_spec: PLUS_TOKEN  */
             { (yyval.ival)=1; }
    break;

  case 107: /* charge_spec: MINUS_TOKEN MINUS_TOKEN  */
                          { (yyval.ival)=-2; }
    break;

  case 108: /* charge_spec: MINUS_TOKEN number  */
                     { (yyval.ival)=-(yyvsp[0].ival); }
    break;

  case 109: /* charge_spec: MINUS_TOKEN  */
              { (yyval.ival)=-1; }
    break;

  case 111: /* ring_number: PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit  */
                                          { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
    break;

  case 112: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN  */
                                                         { (yyval.ival) = (yyvsp[-1].ival); }
    break;

  case 113: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN  */
                                                               { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 114: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN  */
                                                                     { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 115: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN  */
                                                                           { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 116: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN  */
                                                                                 { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 120: /* nonzero_number: nonzero_number digit  */
                       {
    if((yyvsp[-1].ival) >= std::numeric_limits<std::int32_t>::max()/10 ||
     (yyvsp[-1].ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yyvsp[0].ival) ){
     yysmarts_error(input,molList,lastAtom,lastBond,numAtomsParsed,numBondsParsed,branchPoints,scanner,start_token, current_token_position, "number too large");
     YYABORT;
  }
  (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
    break;

  case 123: /* branch_open_token: GROUP_OPEN_TOKEN  */
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


