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
#line 1 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"


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
                std::list<unsigned int> *,
		void *,int , const char *msg  )
{
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error: " << msg << " while parsing: " << input << std::endl;
}

void
yysmarts_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::list<unsigned int> *,
		void *,int, const char * msg )
{
  yyErrorCleanup(ms);
  BOOST_LOG(rdErrorLog) << "SMARTS Parse Error: " << msg << " while parsing: " << input << std::endl;
}



#line 136 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
  YYSYMBOL_digit = 64                      /* digit  */
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
#define YYLAST   610

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  45
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  20
/* YYNRULES -- Number of rules.  */
#define YYNRULES  121
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  176

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
       0,   127,   127,   130,   134,   137,   140,   144,   148,   151,
     156,   161,   164,   172,   173,   174,   175,   183,   192,   206,
     227,   233,   257,   278,   295,   318,   332,   333,   334,   338,
     361,   365,   370,   376,   385,   391,   399,   407,   420,   426,
     433,   439,   462,   465,   471,   472,   476,   493,   517,   518,
     523,   524,   529,   530,   535,   536,   537,   538,   539,   540,
     541,   544,   547,   550,   553,   556,   559,   565,   571,   578,
     586,   594,   600,   606,   612,   618,   624,   631,   638,   639,
     646,   647,   650,   653,   656,   659,   662,   667,   675,   687,
     692,   697,   701,   705,   709,   712,   713,   720,   721,   727,
     733,   739,   744,   751,   752,   753,   754,   755,   756,   760,
     761,   762,   763,   764,   765,   766,   771,   772,   776,   777,
     786,   787
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
  "nonzero_number", "digit", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-57)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-88)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     238,   -31,    45,   182,    16,    29,   -57,   -57,   -57,   -57,
     404,   494,   -57,   -57,   -57,   -57,   100,   199,   236,   273,
     -57,   379,   416,   -57,   -57,    59,   343,   117,    59,    27,
     219,   256,   441,    45,   256,    59,   -57,    22,   293,   -57,
     -57,   -57,    25,    24,   -57,   105,   130,   -57,   -57,   -57,
      16,   -57,   -57,   131,    16,   -57,   -57,    35,   -57,   -57,
     556,   145,   -57,   325,   -57,   -57,   522,   -57,    45,    18,
     -57,   550,   -57,   -57,   -57,   -57,   -57,   -57,   -57,   -57,
     -57,   -57,   -57,   -57,   -57,   -57,   -57,   256,   -57,   145,
     -57,   -57,   466,   -57,   -57,   -57,   441,   441,   441,   -57,
       9,   -57,    59,    59,   -57,   -57,   -57,    16,    16,    16,
     -57,   -57,   -57,    89,   185,   -57,    59,   -13,   -57,    59,
     571,   -57,    52,   -57,   130,   130,   -57,   -57,    23,   441,
     367,   330,    59,    58,   -57,   -57,   -57,    37,   134,    57,
     -57,    59,    63,   -57,    59,    21,   -57,   -57,     3,    79,
     108,    56,   -57,    92,   -57,   114,   -57,    59,   -57,   124,
     130,   -57,   -57,   137,   -57,   -57,   144,   -57,   258,   -57,
     -57,   -57,   295,   -57,   162,   -57
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     0,     5,     8,     0,    12,    89,    88,    90,
       0,     2,    17,    27,    26,    50,    54,    57,    58,    59,
      78,    55,    56,   116,   118,     0,   108,   105,    72,    75,
       0,     0,     0,     0,     0,    76,     4,     0,    16,    42,
      44,    45,     0,    48,    73,    79,   117,    99,    98,   101,
       0,   100,    97,     7,    94,    95,     1,     0,     9,    11,
      72,     0,    48,    79,   121,   120,     0,    25,     0,     0,
      18,     0,    21,   109,    60,    63,    64,    65,    61,    62,
      52,   106,   107,   103,   104,    71,    74,     0,    13,    16,
      14,    43,     0,    15,    77,     3,     0,     0,     0,    41,
       0,    51,     0,    69,    49,   119,   102,     0,     0,     0,
       6,    96,    10,   108,   105,    30,     0,     0,    28,     0,
      69,    23,     0,    20,     0,     0,    19,    22,    46,    38,
      39,    40,     0,     0,    53,    70,    91,    92,    93,     0,
      34,     0,     0,    32,     0,     0,    24,   110,     0,     0,
       0,     0,    31,     0,    29,     0,    36,     0,   111,     0,
      47,    66,    67,     0,    35,    33,     0,   112,     0,    68,
      37,   113,     0,   114,     0,   115
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -57,   -57,     2,   152,    28,   -57,     4,     8,   -57,   -57,
     -57,    -1,    19,   -57,    -7,   -56,   115,   -10,    49,   -43
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     5,    88,    11,    12,    13,    38,    39,    40,    41,
      42,    62,    71,    54,    55,    44,    72,    45,    46,    73
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      63,    14,    43,   105,   117,    36,    74,    75,    76,    77,
      14,    78,    79,     6,    61,    80,    82,    84,    85,   140,
      63,    64,    65,    53,   158,    94,   141,    23,    24,    56,
      57,    37,    14,    90,    89,   132,    93,   124,   125,    70,
      91,    47,    48,   106,   104,    49,    99,   111,   100,    50,
      85,     7,     8,   156,     9,    51,    86,    52,     7,     8,
     157,     9,   104,   149,   145,    14,    95,    14,   -26,    99,
      14,   107,    58,    59,    23,    24,    10,    23,    24,   112,
     162,   147,   148,    10,   151,   122,   107,   108,   109,   152,
     133,    14,   134,   135,   121,   154,   123,    99,    24,   126,
     129,   130,   131,    82,    84,   159,   139,    23,    24,   142,
     135,     7,     8,   101,     9,    81,   168,   105,    23,    24,
      70,    14,   150,   -80,   164,   172,   136,   137,   138,   174,
     102,   153,   161,   103,   155,    23,    24,    99,    99,    99,
     -86,   163,    64,    65,    83,   167,   165,   166,    64,    65,
     146,     7,     8,    15,     9,    16,    17,    18,    19,    20,
      21,   169,    22,    23,    24,   107,   108,   109,   107,   108,
      25,    26,    27,    28,    29,   110,   170,   118,    32,    96,
      97,    98,    33,   175,   119,    92,   127,    35,     7,     8,
      15,     9,    16,    17,    18,    19,    20,    21,   160,    22,
      23,    24,     0,    23,    24,     0,     0,    25,    26,    27,
      28,    29,    83,    30,    31,    32,     0,    23,    24,    33,
       0,    34,   -83,     0,    35,     7,     8,    15,     9,    16,
      17,    18,    19,    20,    21,     0,    22,    23,    24,     1,
       0,     2,     3,     4,    25,    26,    27,    60,    29,     0,
      87,    31,    32,     0,    23,    24,    33,     0,    34,   -84,
       0,    35,     7,     8,    15,     9,    16,    17,    18,    19,
      20,    21,     0,    22,    23,    24,    64,    65,     0,   171,
       0,    25,    26,    27,    28,    29,     0,    87,    31,    32,
       0,    23,    24,    33,     0,    34,   -85,     0,    35,     7,
       8,    15,     9,    16,    17,    18,    19,    20,    21,     0,
      22,    23,    24,    64,    65,     0,   173,     0,    25,    26,
      27,    28,    29,     0,     0,     0,    32,    96,    97,    98,
      33,     7,     8,   101,     9,    35,     7,     8,    15,     9,
      16,    17,    18,    19,    20,    21,     0,    22,    23,    24,
     102,     0,     0,   120,     0,    25,    26,    27,    28,    29,
       0,    23,    24,    32,    96,    97,   -87,    33,     0,    81,
       0,     0,    35,     7,     8,    15,     9,    16,    17,    18,
      19,    20,    21,     0,    22,    23,    24,     0,     0,     0,
       0,     0,    25,    26,    27,    28,    29,    23,    24,     0,
      32,    96,   -81,     0,    33,     0,     0,     0,     0,    35,
       7,     8,    15,     9,    16,    17,    18,    19,    20,    21,
       0,    22,    23,    24,     0,     0,     0,     0,     0,    25,
      26,    27,    60,    29,    23,    24,     0,    32,     0,   -82,
       0,    33,     0,     0,     0,     0,    35,     7,     8,    15,
       9,    16,    17,    18,    19,    20,    21,     0,    22,    23,
      24,     0,     0,     0,     0,     0,    25,    26,    27,    28,
      29,     0,     7,     8,    32,     9,     0,     0,    33,     0,
       0,     0,     0,    35,    64,    65,    66,    67,    68,     0,
       0,    47,    48,     0,     0,    49,    69,    10,     0,    50,
       7,     8,     0,     9,   128,    51,     0,    52,     0,     0,
       0,     0,    64,    65,    66,    67,    68,     0,     0,    47,
      48,     0,     0,    49,    69,    10,     0,    50,     7,     8,
       0,     9,     0,    51,     0,    52,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    47,    48,     0,
       0,    49,     0,    10,     0,    50,     7,     8,     0,     9,
       0,    51,     0,    52,     0,     0,     0,     0,    64,    65,
       0,     0,     0,     0,    23,    24,     0,     0,     0,     0,
      69,    10,   113,   114,   107,   108,   109,     0,   115,    23,
      24,     0,     0,     0,     0,   116,     0,   113,   114,     0,
       0,     0,     0,   143,     0,     0,     0,     0,     0,     0,
     144
};

static const yytype_int16 yycheck[] =
{
      10,     2,     3,    46,    60,     3,    16,    17,    18,    19,
      11,    21,    22,    44,    10,    25,    26,    27,    28,    32,
      30,    18,    19,     4,    21,    35,    39,    18,    19,     0,
       1,     3,    33,    31,    30,    26,    34,    19,    20,    11,
      32,    25,    26,    50,    45,    29,    38,    54,    23,    33,
      60,     6,     7,    32,     9,    39,    29,    41,     6,     7,
      39,     9,    63,    40,   120,    66,    44,    68,    44,    61,
      71,    34,    43,    44,    18,    19,    31,    18,    19,    44,
      24,   124,   125,    31,    26,    66,    34,    35,    36,    32,
     100,    92,   102,   103,    66,    32,    68,    89,    19,    71,
      96,    97,    98,   113,   114,   148,   116,    18,    19,   119,
     120,     6,     7,     8,     9,    26,   159,   160,    18,    19,
      92,   122,   132,    23,    32,   168,   107,   108,   109,   172,
      25,   141,    24,    28,   144,    18,    19,   129,   130,   131,
      23,   151,    18,    19,    27,    21,    32,   157,    18,    19,
     122,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    24,    17,    18,    19,    34,    35,    36,    34,    35,
      25,    26,    27,    28,    29,    44,    32,    32,    33,    34,
      35,    36,    37,    21,    39,    33,    71,    42,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,   149,    17,
      18,    19,    -1,    18,    19,    -1,    -1,    25,    26,    27,
      28,    29,    27,    31,    32,    33,    -1,    18,    19,    37,
      -1,    39,    23,    -1,    42,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    -1,    17,    18,    19,     1,
      -1,     3,     4,     5,    25,    26,    27,    28,    29,    -1,
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
      -1,    18,    19,    33,    34,    35,    23,    37,    -1,    26,
      -1,    -1,    42,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    -1,    17,    18,    19,    -1,    -1,    -1,
      -1,    -1,    25,    26,    27,    28,    29,    18,    19,    -1,
      33,    34,    23,    -1,    37,    -1,    -1,    -1,    -1,    42,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      -1,    17,    18,    19,    -1,    -1,    -1,    -1,    -1,    25,
      26,    27,    28,    29,    18,    19,    -1,    33,    -1,    23,
      -1,    37,    -1,    -1,    -1,    -1,    42,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    -1,    17,    18,
      19,    -1,    -1,    -1,    -1,    -1,    25,    26,    27,    28,
      29,    -1,     6,     7,    33,     9,    -1,    -1,    37,    -1,
      -1,    -1,    -1,    42,    18,    19,    20,    21,    22,    -1,
      -1,    25,    26,    -1,    -1,    29,    30,    31,    -1,    33,
       6,     7,    -1,     9,    38,    39,    -1,    41,    -1,    -1,
      -1,    -1,    18,    19,    20,    21,    22,    -1,    -1,    25,
      26,    -1,    -1,    29,    30,    31,    -1,    33,     6,     7,
      -1,     9,    -1,    39,    -1,    41,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    25,    26,    -1,
      -1,    29,    -1,    31,    -1,    33,     6,     7,    -1,     9,
      -1,    39,    -1,    41,    -1,    -1,    -1,    -1,    18,    19,
      -1,    -1,    -1,    -1,    18,    19,    -1,    -1,    -1,    -1,
      30,    31,    26,    27,    34,    35,    36,    -1,    32,    18,
      19,    -1,    -1,    -1,    -1,    39,    -1,    26,    27,    -1,
      -1,    -1,    -1,    32,    -1,    -1,    -1,    -1,    -1,    -1,
      39
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
      49,    57,    61,    64,    62,    62,    62,    62,    62,    62,
      62,    26,    62,    27,    62,    62,    29,    31,    47,    51,
      47,    52,    48,    47,    62,    44,    34,    35,    36,    52,
      23,     8,    25,    28,    56,    64,    59,    34,    35,    36,
      44,    59,    44,    26,    27,    32,    39,    60,    32,    39,
      28,    49,    57,    49,    19,    20,    49,    61,    38,    51,
      51,    51,    26,    62,    62,    62,    57,    57,    57,    62,
      32,    39,    62,    32,    39,    60,    49,    64,    64,    40,
      62,    26,    32,    62,    32,    62,    32,    39,    21,    64,
      63,    24,    24,    62,    32,    32,    62,    21,    64,    24,
      32,    21,    64,    21,    64,    21
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    45,    46,    46,    46,    46,    46,    46,    46,    46,
      46,    46,    46,    47,    47,    47,    47,    48,    48,    48,
      48,    48,    48,    48,    48,    48,    49,    49,    49,    49,
      50,    50,    50,    50,    50,    50,    50,    50,    51,    51,
      51,    51,    51,    52,    52,    52,    53,    53,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      55,    55,    55,    55,    55,    55,    55,    55,    56,    56,
      56,    57,    57,    57,    57,    58,    58,    59,    59,    59,
      59,    59,    59,    60,    60,    60,    60,    60,    60,    61,
      61,    61,    61,    61,    61,    61,    62,    62,    63,    63,
      64,    64
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     3,     2,     1,     3,     2,     1,     2,
       3,     2,     2,     2,     2,     2,     1,     1,     2,     3,
       3,     2,     3,     3,     4,     2,     1,     1,     3,     5,
       3,     5,     4,     6,     4,     6,     5,     7,     3,     3,
       3,     2,     1,     2,     1,     1,     3,     5,     1,     2,
       1,     2,     2,     3,     1,     1,     1,     1,     1,     1,
       2,     2,     2,     2,     2,     2,     5,     5,     6,     2,
       3,     2,     1,     1,     2,     1,     1,     2,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     3,     3,     1,     1,     2,     1,     1,     1,
       1,     1,     2,     2,     2,     1,     2,     2,     1,     1,
       3,     4,     5,     6,     7,     8,     1,     1,     1,     2,
       1,     1
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
        yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, YY_("syntax error: cannot back up")); \
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
                  Kind, Value, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
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
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
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
                 int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
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
                       &yyvsp[(yyi + 1) - (yynrhs)], input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
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
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, RDKit::Atom* &lastAtom, RDKit::Bond* &lastBond, unsigned &numAtomsParsed, unsigned &numBondsParsed, std::list<unsigned int> *branchPoints, void *scanner, int& start_token)
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
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yykind)
    {
    case YYSYMBOL_ATOM_TOKEN: /* ATOM_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1162 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_SIMPLE_ATOM_QUERY_TOKEN: /* SIMPLE_ATOM_QUERY_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1168 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_COMPLEX_ATOM_QUERY_TOKEN: /* COMPLEX_ATOM_QUERY_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1174 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_RINGSIZE_ATOM_QUERY_TOKEN: /* RINGSIZE_ATOM_QUERY_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1180 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_RINGBOND_ATOM_QUERY_TOKEN: /* RINGBOND_ATOM_QUERY_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1186 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_IMPLICIT_H_ATOM_QUERY_TOKEN: /* IMPLICIT_H_ATOM_QUERY_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1192 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_HYB_TOKEN: /* HYB_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1198 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_HETERONEIGHBOR_ATOM_QUERY_TOKEN: /* HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1204 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_ALIPHATIC: /* ALIPHATIC  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1210 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN: /* ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1216 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_BOND_TOKEN: /* BOND_TOKEN  */
#line 119 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1222 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_atomd: /* atomd  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1228 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_hydrogen_atom: /* hydrogen_atom  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1234 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_atom_expr: /* atom_expr  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1240 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_point_query: /* point_query  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1246 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_recursive_query: /* recursive_query  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1252 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_atom_query: /* atom_query  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1258 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_possible_range_query: /* possible_range_query  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1264 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_simple_atom: /* simple_atom  */
#line 118 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).atom); }
#line 1270 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_bond_expr: /* bond_expr  */
#line 119 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1276 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_bond_query: /* bond_query  */
#line 119 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1282 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
        break;

    case YYSYMBOL_bondd: /* bondd  */
#line 119 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            { delete ((*yyvaluep).bond); }
#line 1288 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
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
      yychar = yylex (&yylval, scanner, start_token);
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
#line 127 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
              {
// the molList has already been updated, no need to do anything
}
#line 1566 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 3: /* meta_start: START_ATOM atomd EOS_TOKEN  */
#line 130 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                             {
  lastAtom = (yyvsp[-1].atom);
  YYACCEPT;
}
#line 1575 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 4: /* meta_start: START_ATOM bad_atom_def  */
#line 134 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                          {
  YYABORT;
}
#line 1583 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 5: /* meta_start: START_ATOM  */
#line 137 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
             {
  YYABORT;
}
#line 1591 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 6: /* meta_start: START_BOND bond_expr EOS_TOKEN  */
#line 140 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                 {
  lastBond = (yyvsp[-1].bond);
  YYACCEPT;
}
#line 1600 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 7: /* meta_start: START_BOND bond_expr  */
#line 144 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                       {
  delete (yyvsp[0].bond);
  YYABORT;
}
#line 1609 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 8: /* meta_start: START_BOND  */
#line 148 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
             {
  YYABORT;
}
#line 1617 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 9: /* meta_start: meta_start BAD_CHARACTER  */
#line 151 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                           {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1627 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 10: /* meta_start: meta_start error EOS_TOKEN  */
#line 156 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                            {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1637 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 11: /* meta_start: meta_start EOS_TOKEN  */
#line 161 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                       {
  YYACCEPT;
}
#line 1645 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 12: /* meta_start: error EOS_TOKEN  */
#line 164 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                  {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
#line 1655 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 16: /* bad_atom_def: atom_expr  */
#line 175 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
            {
  delete (yyvsp[0].atom);
  YYABORT;
}
#line 1664 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 17: /* mol: atomd  */
#line 183 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
           {
  int sz     = molList->size();
  molList->resize( sz + 1);
  (*molList)[ sz ] = new RWMol();
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1);
  (*molList)[ sz ]->addAtom((yyvsp[0].atom),true,true);
  //delete $1;
  (yyval.moli) = sz;
}
#line 1678 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 18: /* mol: mol atomd  */
#line 192 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 1696 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 19: /* mol: mol bond_expr atomd  */
#line 206 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 1721 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 20: /* mol: mol SEPARATOR_TOKEN atomd  */
#line 227 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                            {
  RWMol *mp = (*molList)[(yyval.moli)];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom((yyvsp[0].atom),true,true);
}
#line 1731 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 21: /* mol: mol ring_number  */
#line 233 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 1759 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 22: /* mol: mol bond_expr ring_number  */
#line 257 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 1784 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 23: /* mol: mol GROUP_OPEN_TOKEN atomd  */
#line 278 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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

  branchPoints->push_back(atomIdx1);
}
#line 1805 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 24: /* mol: mol GROUP_OPEN_TOKEN bond_expr atomd  */
#line 295 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
  branchPoints->push_back(atomIdx1);

}
#line 1831 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 25: /* mol: mol GROUP_CLOSE_TOKEN  */
#line 318 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 1846 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 28: /* atomd: ATOM_OPEN_TOKEN atom_expr ATOM_CLOSE_TOKEN  */
#line 335 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
{
  (yyval.atom) = (yyvsp[-1].atom);
}
#line 1854 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 29: /* atomd: ATOM_OPEN_TOKEN atom_expr COLON_TOKEN number ATOM_CLOSE_TOKEN  */
#line 339 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
{
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1863 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 30: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN ATOM_CLOSE_TOKEN  */
#line 362 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
{
  (yyval.atom) = new QueryAtom(1);
}
#line 1871 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 31: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN  */
#line 366 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
{
  (yyval.atom) = new QueryAtom(1);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
#line 1880 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 32: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN ATOM_CLOSE_TOKEN  */
#line 370 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                  {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1891 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 33: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN  */
#line 376 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                     {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-4].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-4].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1904 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 34: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN charge_spec ATOM_CLOSE_TOKEN  */
#line 385 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                       {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1915 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 35: /* hydrogen_atom: ATOM_OPEN_TOKEN H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN  */
#line 391 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                          {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1928 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 36: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN charge_spec ATOM_CLOSE_TOKEN  */
#line 399 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                              {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-3].ival));
  newQ->setFormalCharge((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=newQ;
}
#line 1941 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 37: /* hydrogen_atom: ATOM_OPEN_TOKEN number H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN  */
#line 407 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                                 {
  QueryAtom *newQ = new QueryAtom(1);
  newQ->setIsotope((yyvsp[-5].ival));
  newQ->setFormalCharge((yyvsp[-3].ival));
  newQ->expandQuery(makeAtomIsotopeQuery((yyvsp[-5].ival)),Queries::COMPOSITE_AND,true);
  newQ->expandQuery(makeAtomFormalChargeQuery((yyvsp[-3].ival)),Queries::COMPOSITE_AND,true);
  newQ->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));

  (yyval.atom)=newQ;
}
#line 1956 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 38: /* atom_expr: atom_expr AND_TOKEN atom_expr  */
#line 420 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                         {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 1967 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 39: /* atom_expr: atom_expr OR_TOKEN atom_expr  */
#line 426 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                               {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  (yyvsp[-2].atom)->setAtomicNum(0);
  delete (yyvsp[0].atom);
}
#line 1979 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 40: /* atom_expr: atom_expr SEMI_TOKEN atom_expr  */
#line 433 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                 {
  (yyvsp[-2].atom)->expandQuery((yyvsp[0].atom)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  if((yyvsp[-2].atom)->getChiralTag()==Atom::CHI_UNSPECIFIED) (yyvsp[-2].atom)->setChiralTag((yyvsp[0].atom)->getChiralTag());
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[-2].atom));
  delete (yyvsp[0].atom);
}
#line 1990 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 41: /* atom_expr: atom_expr point_query  */
#line 439 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 2018 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 43: /* point_query: NOT_TOKEN point_query  */
#line 465 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                   {
  (yyvsp[0].atom)->getQuery()->setNegation(!((yyvsp[0].atom)->getQuery()->getNegation()));
  (yyvsp[0].atom)->setAtomicNum(0);
  SmilesParseOps::ClearAtomChemicalProps((yyvsp[0].atom));
  (yyval.atom) = (yyvsp[0].atom);
}
#line 2029 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 46: /* recursive_query: BEGIN_RECURSE mol END_RECURSE  */
#line 476 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 2051 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 47: /* recursive_query: BEGIN_RECURSE mol END_RECURSE UNDERSCORE_TOKEN nonzero_number  */
#line 493 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 2077 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 49: /* atom_query: number simple_atom  */
#line 518 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                     {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2087 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 51: /* atom_query: number ATOM_TOKEN  */
#line 524 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                    {
  (yyvsp[0].atom)->setIsotope((yyvsp[-1].ival));
  (yyvsp[0].atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)),Queries::COMPOSITE_AND,true);
  (yyval.atom)=(yyvsp[0].atom);
}
#line 2097 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 52: /* atom_query: HASH_TOKEN number  */
#line 529 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                    { (yyval.atom) = new QueryAtom((yyvsp[0].ival)); }
#line 2103 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 53: /* atom_query: number HASH_TOKEN number  */
#line 530 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                           {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsotope((yyvsp[-2].ival));
  (yyval.atom)->expandQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)),Queries::COMPOSITE_AND,true);
}
#line 2113 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 60: /* atom_query: COMPLEX_ATOM_QUERY_TOKEN number  */
#line 541 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                  {
  static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-1].atom)->getQuery())->setVal((yyvsp[0].ival));
}
#line 2121 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 61: /* atom_query: HETERONEIGHBOR_ATOM_QUERY_TOKEN number  */
#line 544 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                         {
  (yyvsp[-1].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2129 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 62: /* atom_query: ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN number  */
#line 547 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                  {
  (yyvsp[-1].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery((yyvsp[0].ival)));
}
#line 2137 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 63: /* atom_query: RINGSIZE_ATOM_QUERY_TOKEN number  */
#line 550 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomMinRingSizeQuery((yyvsp[0].ival)));
}
#line 2145 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 64: /* atom_query: RINGBOND_ATOM_QUERY_TOKEN number  */
#line 553 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                   {
  (yyvsp[-1].atom)->setQuery(makeAtomRingBondCountQuery((yyvsp[0].ival)));
}
#line 2153 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 65: /* atom_query: IMPLICIT_H_ATOM_QUERY_TOKEN number  */
#line 556 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                     {
  (yyvsp[-1].atom)->setQuery(makeAtomImplicitHCountQuery((yyvsp[0].ival)));
}
#line 2161 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 66: /* atom_query: possible_range_query RANGE_OPEN_TOKEN MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
#line 559 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>((yyvsp[-1].ival),oq->getDataFunc(),
    std::string("greater_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2172 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 67: /* atom_query: possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN RANGE_CLOSE_TOKEN  */
#line 565 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                             {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-4].atom)->getQuery());
  ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>((yyvsp[-2].ival),oq->getDataFunc(),
    std::string("less_")+oq->getDescription());
  (yyvsp[-4].atom)->setQuery(nq);
}
#line 2183 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 68: /* atom_query: possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN number RANGE_CLOSE_TOKEN  */
#line 571 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                                    {
  ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>((yyvsp[-5].atom)->getQuery());
  ATOM_RANGE_QUERY *nq = makeAtomRangeQuery((yyvsp[-3].ival),(yyvsp[-1].ival),false,false,
    oq->getDataFunc(),
    std::string("range_")+oq->getDescription());
  (yyvsp[-5].atom)->setQuery(nq);
}
#line 2195 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 69: /* atom_query: number H_TOKEN  */
#line 578 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-1].ival)));
  newQ->setIsotope((yyvsp[-1].ival));
  newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2208 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 70: /* atom_query: number H_TOKEN number  */
#line 586 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                        {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[-2].ival)));
  newQ->setIsotope((yyvsp[-2].ival));
  newQ->expandQuery(makeAtomHCountQuery((yyvsp[0].ival)),Queries::COMPOSITE_AND,true);
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2221 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 71: /* atom_query: H_TOKEN number  */
#line 594 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                 {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery((yyvsp[0].ival)));
  newQ->setNumExplicitHs((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2232 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 72: /* atom_query: H_TOKEN  */
#line 600 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
          {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomHCountQuery(1));
  newQ->setNumExplicitHs(1);
  (yyval.atom)=newQ;
}
#line 2243 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 73: /* atom_query: charge_spec  */
#line 606 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery((yyvsp[0].ival)));
  newQ->setFormalCharge((yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2254 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 74: /* atom_query: AT_TOKEN AT_TOKEN  */
#line 612 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                    {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  (yyval.atom)=newQ;
}
#line 2265 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 75: /* atom_query: AT_TOKEN  */
#line 618 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
           {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
  (yyval.atom)=newQ;
}
#line 2276 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 76: /* atom_query: CHI_CLASS_TOKEN  */
#line 624 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                  {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[0].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,0);
  (yyval.atom)=newQ;
}
#line 2288 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 77: /* atom_query: CHI_CLASS_TOKEN number  */
#line 631 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNullQuery());
  newQ->setChiralTag((yyvsp[-1].chiraltype));
  newQ->setProp(common_properties::_chiralPermutation,(yyvsp[0].ival));
  (yyval.atom)=newQ;
}
#line 2300 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 79: /* atom_query: number  */
#line 639 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
         {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomIsotopeQuery((yyvsp[0].ival)));
  (yyval.atom)=newQ;
}
#line 2310 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 81: /* possible_range_query: HETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 647 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                  {
  (yyvsp[0].atom)->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
}
#line 2318 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 82: /* possible_range_query: ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN  */
#line 650 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                           {
  (yyvsp[0].atom)->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
}
#line 2326 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 83: /* possible_range_query: RINGSIZE_ATOM_QUERY_TOKEN  */
#line 653 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                            {
  (yyvsp[0].atom)->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
}
#line 2334 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 84: /* possible_range_query: RINGBOND_ATOM_QUERY_TOKEN  */
#line 656 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                            {
  (yyvsp[0].atom)->setQuery(makeAtomRingBondCountQuery(0));
}
#line 2342 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 85: /* possible_range_query: IMPLICIT_H_ATOM_QUERY_TOKEN  */
#line 659 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                              {
  (yyvsp[0].atom)->setQuery(makeAtomImplicitHCountQuery(0));
}
#line 2350 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 86: /* possible_range_query: PLUS_TOKEN  */
#line 662 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
             {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
#line 2360 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 87: /* possible_range_query: MINUS_TOKEN  */
#line 667 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
              {
  QueryAtom *newQ = new QueryAtom();
  newQ->setQuery(makeAtomNegativeFormalChargeQuery(0));
  (yyval.atom) = newQ;
}
#line 2370 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 88: /* simple_atom: ORGANIC_ATOM_TOKEN  */
#line 675 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
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
#line 2387 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 89: /* simple_atom: AROMATIC_ATOM_TOKEN  */
#line 687 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                      {
  (yyval.atom) = new QueryAtom((yyvsp[0].ival));
  (yyval.atom)->setIsAromatic(true);
  (yyval.atom)->setQuery(makeAtomTypeQuery((yyvsp[0].ival),true));
}
#line 2397 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 91: /* bond_expr: bond_expr AND_TOKEN bond_expr  */
#line 697 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                        {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2406 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 92: /* bond_expr: bond_expr OR_TOKEN bond_expr  */
#line 701 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                               {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_OR,true);
  delete (yyvsp[0].bond);
}
#line 2415 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 93: /* bond_expr: bond_expr SEMI_TOKEN bond_expr  */
#line 705 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                 {
  (yyvsp[-2].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2424 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 96: /* bond_query: bond_query bondd  */
#line 713 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                   {
  (yyvsp[-1].bond)->expandQuery((yyvsp[0].bond)->getQuery()->copy(),Queries::COMPOSITE_AND,true);
  delete (yyvsp[0].bond);
}
#line 2433 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 98: /* bondd: MINUS_TOKEN  */
#line 721 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::SINGLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
  (yyval.bond) = newB;
}
#line 2444 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 99: /* bondd: HASH_TOKEN  */
#line 727 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
             {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::TRIPLE);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
  (yyval.bond) = newB;
}
#line 2455 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 100: /* bondd: COLON_TOKEN  */
#line 733 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
              {
  QueryBond *newB= new QueryBond();
  newB->setBondType(Bond::AROMATIC);
  newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
  (yyval.bond) = newB;
}
#line 2466 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 101: /* bondd: AT_TOKEN  */
#line 739 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
           {
  QueryBond *newB= new QueryBond();
  newB->setQuery(makeBondIsInRingQuery());
  (yyval.bond) = newB;
}
#line 2476 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 102: /* bondd: NOT_TOKEN bondd  */
#line 744 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                  {
  (yyvsp[0].bond)->getQuery()->setNegation(!((yyvsp[0].bond)->getQuery()->getNegation()));
  (yyval.bond) = (yyvsp[0].bond);
}
#line 2485 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 103: /* charge_spec: PLUS_TOKEN PLUS_TOKEN  */
#line 751 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                   { (yyval.ival)=2; }
#line 2491 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 104: /* charge_spec: PLUS_TOKEN number  */
#line 752 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                    { (yyval.ival)=(yyvsp[0].ival); }
#line 2497 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 105: /* charge_spec: PLUS_TOKEN  */
#line 753 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
             { (yyval.ival)=1; }
#line 2503 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 106: /* charge_spec: MINUS_TOKEN MINUS_TOKEN  */
#line 754 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                          { (yyval.ival)=-2; }
#line 2509 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 107: /* charge_spec: MINUS_TOKEN number  */
#line 755 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                     { (yyval.ival)=-(yyvsp[0].ival); }
#line 2515 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 108: /* charge_spec: MINUS_TOKEN  */
#line 756 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
              { (yyval.ival)=-1; }
#line 2521 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 110: /* ring_number: PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit  */
#line 761 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                          { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
#line 2527 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 111: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN  */
#line 762 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                         { (yyval.ival) = (yyvsp[-1].ival); }
#line 2533 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 112: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN  */
#line 763 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                               { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2539 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 113: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN  */
#line 764 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                     { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2545 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 114: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN  */
#line 765 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                           { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2551 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 115: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN  */
#line 766 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                                                                                 { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
#line 2557 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;

  case 119: /* nonzero_number: nonzero_number digit  */
#line 777 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"
                       { 
    if((yyvsp[-1].ival) >= std::numeric_limits<std::int32_t>::max()/10 || 
     (yyvsp[-1].ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yyvsp[0].ival) ){
     yysmarts_error(input,molList,lastAtom,lastBond,numAtomsParsed,numBondsParsed,branchPoints,scanner,start_token,"number too large");
     YYABORT;
  }
  (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival); }
#line 2569 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"
    break;


#line 2573 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.tab.cpp"

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
      yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, YY_("syntax error"));
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
                  YY_ACCESSING_SYMBOL (yystate), yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
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
  yyerror (input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token, YY_("memory exhausted"));
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
                  yytoken, &yylval, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp, input, molList, lastAtom, lastBond, numAtomsParsed, numBondsParsed, branchPoints, scanner, start_token);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 790 "/home/glandrum/RDKit_git/Code/GraphMol/SmilesParse/smarts.yy"

