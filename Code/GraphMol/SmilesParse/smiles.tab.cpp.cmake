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
#define yyparse         yysmiles_parse
#define yylex           yysmiles_lex
#define yyerror         yysmiles_error
#define yydebug         yysmiles_debug
#define yynerrs         yysmiles_nerrs

/* First part of user prologue.  */


  // $Id$
  //
  //  Copyright (C) 2001-2016 Randal Henne, Greg Landrum and Rational Discovery LLC
  //
  //   @@ All Rights Reserved  @@
  //

#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesParseOps.h>
#include <RDGeneral/RDLog.h>

#define YYDEBUG 1
#include "smiles.tab.hpp"

extern int yysmiles_lex(YYSTYPE *,void *,int &, unsigned int&);

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
                std::vector<std::pair<unsigned int, unsigned int>>&,
		void *,int, unsigned int bad_token_position, const char * msg )
{
  yyErrorCleanup(ms);
  SmilesParseOps::detail::printSyntaxErrorMessage(input, msg, bad_token_position);
}

void
yysmiles_error( const char *input,
                std::vector<RDKit::RWMol *> *ms,
                std::vector<std::pair<unsigned int, unsigned int>>&,
		void *,int, unsigned int bad_token_position, const char * msg )
{
  yyErrorCleanup(ms);
  SmilesParseOps::detail::printSyntaxErrorMessage(input, msg, bad_token_position);
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

#include "smiles.tab.hpp"
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
  YYSYMBOL_ATOM_TOKEN = 7,                 /* ATOM_TOKEN  */
  YYSYMBOL_ORGANIC_ATOM_TOKEN = 8,         /* ORGANIC_ATOM_TOKEN  */
  YYSYMBOL_NONZERO_DIGIT_TOKEN = 9,        /* NONZERO_DIGIT_TOKEN  */
  YYSYMBOL_ZERO_TOKEN = 10,                /* ZERO_TOKEN  */
  YYSYMBOL_GROUP_OPEN_TOKEN = 11,          /* GROUP_OPEN_TOKEN  */
  YYSYMBOL_GROUP_CLOSE_TOKEN = 12,         /* GROUP_CLOSE_TOKEN  */
  YYSYMBOL_SEPARATOR_TOKEN = 13,           /* SEPARATOR_TOKEN  */
  YYSYMBOL_LOOP_CONNECTOR_TOKEN = 14,      /* LOOP_CONNECTOR_TOKEN  */
  YYSYMBOL_MINUS_TOKEN = 15,               /* MINUS_TOKEN  */
  YYSYMBOL_PLUS_TOKEN = 16,                /* PLUS_TOKEN  */
  YYSYMBOL_H_TOKEN = 17,                   /* H_TOKEN  */
  YYSYMBOL_AT_TOKEN = 18,                  /* AT_TOKEN  */
  YYSYMBOL_PERCENT_TOKEN = 19,             /* PERCENT_TOKEN  */
  YYSYMBOL_COLON_TOKEN = 20,               /* COLON_TOKEN  */
  YYSYMBOL_HASH_TOKEN = 21,                /* HASH_TOKEN  */
  YYSYMBOL_BOND_TOKEN = 22,                /* BOND_TOKEN  */
  YYSYMBOL_CHI_CLASS_TOKEN = 23,           /* CHI_CLASS_TOKEN  */
  YYSYMBOL_ATOM_OPEN_TOKEN = 24,           /* ATOM_OPEN_TOKEN  */
  YYSYMBOL_ATOM_CLOSE_TOKEN = 25,          /* ATOM_CLOSE_TOKEN  */
  YYSYMBOL_EOS_TOKEN = 26,                 /* EOS_TOKEN  */
  YYSYMBOL_YYACCEPT = 27,                  /* $accept  */
  YYSYMBOL_meta_start = 28,                /* meta_start  */
  YYSYMBOL_bad_atom_def = 29,              /* bad_atom_def  */
  YYSYMBOL_mol = 30,                       /* mol  */
  YYSYMBOL_branch_open_token = 31,         /* branch_open_token  */
  YYSYMBOL_bondd = 32,                     /* bondd  */
  YYSYMBOL_atomd = 33,                     /* atomd  */
  YYSYMBOL_charge_element = 34,            /* charge_element  */
  YYSYMBOL_h_element = 35,                 /* h_element  */
  YYSYMBOL_chiral_element = 36,            /* chiral_element  */
  YYSYMBOL_element = 37,                   /* element  */
  YYSYMBOL_simple_atom = 38,               /* simple_atom  */
  YYSYMBOL_ring_number = 39,               /* ring_number  */
  YYSYMBOL_number = 40,                    /* number  */
  YYSYMBOL_nonzero_number = 41,            /* nonzero_number  */
  YYSYMBOL_digit = 42                      /* digit  */
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
#define YYFINAL  33
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   151

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  27
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  16
/* YYNRULES -- Number of rules.  */
#define YYNRULES  74
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  108

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   281


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
      25,    26
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   120,   120,   123,   127,   130,   134,   138,   141,   146,
     149,   157,   158,   159,   160,   168,   179,   190,   211,   220,
     226,   247,   271,   290,   300,   321,   329,   343,   346,   347,
     353,   354,   360,   368,   369,   370,   371,   372,   373,   374,
     378,   379,   380,   381,   382,   383,   384,   385,   386,   390,
     391,   392,   393,   394,   398,   399,   400,   401,   402,   403,
     407,   408,   412,   413,   414,   415,   416,   417,   418,   422,
     423,   427,   428,   439,   440
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
  "START_ATOM", "START_BOND", "AROMATIC_ATOM_TOKEN", "ATOM_TOKEN",
  "ORGANIC_ATOM_TOKEN", "NONZERO_DIGIT_TOKEN", "ZERO_TOKEN",
  "GROUP_OPEN_TOKEN", "GROUP_CLOSE_TOKEN", "SEPARATOR_TOKEN",
  "LOOP_CONNECTOR_TOKEN", "MINUS_TOKEN", "PLUS_TOKEN", "H_TOKEN",
  "AT_TOKEN", "PERCENT_TOKEN", "COLON_TOKEN", "HASH_TOKEN", "BOND_TOKEN",
  "CHI_CLASS_TOKEN", "ATOM_OPEN_TOKEN", "ATOM_CLOSE_TOKEN", "EOS_TOKEN",
  "$accept", "meta_start", "bad_atom_def", "mol", "branch_open_token",
  "bondd", "atomd", "charge_element", "h_element", "chiral_element",
  "element", "simple_atom", "ring_number", "number", "nonzero_number",
  "digit", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-22)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-31)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     128,    -8,    51,     5,    45,     1,   -22,   -22,   -22,    84,
      98,   -22,   -22,   -22,   -22,   -22,    19,    78,    14,    78,
      78,   -22,    15,   -22,    81,    28,    -2,    30,   117,   109,
     -22,   -22,    36,   -22,    39,   -22,    12,   -22,   -22,   -22,
     -22,   -22,    51,    25,    63,    25,   106,   -22,   -22,   -22,
      14,    78,   -22,   -22,   -22,    12,   -22,   -22,    -5,    54,
      14,    50,    14,   -22,    72,    14,   -22,   -22,   -22,   -22,
      14,   -22,   -22,   -22,   -22,   109,   109,   -22,   -22,    51,
      51,   -22,   -22,   -22,   -22,   -22,   -22,   -22,   -22,   -22,
      14,   -22,    75,   -22,   127,   -22,   -22,   -22,   -22,   -22,
     131,   -22,   135,   -22,   139,   -22,   103,   -22
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     0,     0,     7,     0,    10,    61,    60,     0,
       2,    15,    30,    56,    71,    69,    40,     0,     0,     0,
       0,     4,     0,    14,    33,    46,    49,    54,     0,    70,
      29,    28,     6,     1,     0,     9,     0,    54,    73,    74,
      27,    26,     0,     0,     0,     0,     0,    16,    20,    62,
      42,     0,    13,    58,    11,    14,    12,     3,    37,    34,
      47,    50,    52,    57,    41,     0,    55,    72,     5,     8,
       0,    32,    19,    18,    22,     0,     0,    17,    21,     0,
       0,    23,    44,    38,    39,    35,    36,    48,    51,    53,
      43,    59,     0,    63,     0,    25,    24,    45,    31,    64,
       0,    65,     0,    66,     0,    67,     0,    68
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -22,   -22,    49,   -22,   -22,   -22,    -3,    52,   -22,   -22,
     -22,     0,    37,   -12,   -22,   -21
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     5,    54,    10,    46,    32,    11,    23,    24,    25,
      26,    12,    48,    28,    29,    49
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      22,    33,    34,    27,    14,    15,    53,    47,    67,    37,
      83,     7,    13,     8,    14,    15,    61,    37,     6,    37,
      37,    62,    16,    14,    15,    17,    18,    35,    66,    19,
      20,     7,    70,     8,    38,    39,    50,    71,    82,    72,
      73,    57,    77,    81,    44,    60,    84,    86,    87,     9,
      89,    37,    21,    91,    93,    94,   -30,     7,    92,     8,
      30,    36,    68,    14,    15,    69,    52,    31,    88,    56,
      85,    55,    75,   100,    76,     9,    95,    96,    97,   102,
      74,   104,    78,   106,     7,    13,     8,    14,    15,    90,
       7,    13,     8,    14,    15,    16,    58,    59,    17,    18,
      98,    16,    51,    20,     7,    18,     8,    38,    39,    40,
      41,    42,     7,    43,     8,   107,     0,    44,    38,    39,
      45,    79,     9,     7,    63,     8,     0,     0,    80,     1,
       9,     2,     3,     4,    64,     0,    38,    39,    65,    99,
      38,    39,     0,   101,    38,    39,     0,   103,    38,    39,
       0,   105
};

static const yytype_int8 yycheck[] =
{
       3,     0,     1,     3,     9,    10,    18,    10,    29,     9,
      15,     6,     7,     8,     9,    10,    18,    17,    26,    19,
      20,    23,    17,     9,    10,    20,    21,    26,    28,    24,
      25,     6,    20,     8,     9,    10,    17,    25,    50,    42,
      43,    26,    45,    46,    19,    17,    58,    59,    60,    24,
      62,    51,     3,    65,    75,    76,    26,     6,    70,     8,
      15,     9,    26,     9,    10,    26,    17,    22,    18,    20,
      16,    19,     9,    94,    11,    24,    79,    80,    90,   100,
      43,   102,    45,   104,     6,     7,     8,     9,    10,    17,
       6,     7,     8,     9,    10,    17,    15,    16,    20,    21,
      25,    17,    24,    25,     6,    21,     8,     9,    10,    11,
      12,    13,     6,    15,     8,    12,    -1,    19,     9,    10,
      22,    15,    24,     6,     7,     8,    -1,    -1,    22,     1,
      24,     3,     4,     5,    17,    -1,     9,    10,    21,    12,
       9,    10,    -1,    12,     9,    10,    -1,    12,     9,    10,
      -1,    12
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     3,     4,     5,    28,    26,     6,     8,    24,
      30,    33,    38,     7,     9,    10,    17,    20,    21,    24,
      25,    29,    33,    34,    35,    36,    37,    38,    40,    41,
      15,    22,    32,     0,     1,    26,    34,    38,     9,    10,
      11,    12,    13,    15,    19,    22,    31,    33,    39,    42,
      17,    24,    29,    40,    29,    34,    29,    26,    15,    16,
      17,    18,    23,     7,    17,    21,    38,    42,    26,    26,
      20,    25,    33,    33,    39,     9,    11,    33,    39,    15,
      22,    33,    40,    15,    40,    16,    40,    40,    18,    40,
      17,    40,    40,    42,    42,    33,    33,    40,    25,    12,
      42,    12,    42,    12,    42,    12,    42,    12
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    27,    28,    28,    28,    28,    28,    28,    28,    28,
      28,    29,    29,    29,    29,    30,    30,    30,    30,    30,
      30,    30,    30,    30,    30,    30,    30,    31,    32,    32,
      33,    33,    33,    34,    34,    34,    34,    34,    34,    34,
      35,    35,    35,    35,    35,    35,    35,    35,    35,    36,
      36,    36,    36,    36,    37,    37,    37,    37,    37,    37,
      38,    38,    39,    39,    39,    39,    39,    39,    39,    40,
      40,    41,    41,    42,    42
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     3,     2,     3,     2,     1,     3,     2,
       2,     2,     2,     2,     1,     1,     2,     3,     3,     3,
       2,     3,     3,     3,     4,     4,     2,     1,     1,     1,
       1,     5,     3,     1,     2,     3,     3,     2,     3,     3,
       1,     2,     2,     3,     3,     4,     1,     2,     3,     1,
       2,     3,     2,     3,     1,     2,     1,     2,     2,     3,
       1,     1,     1,     3,     4,     5,     6,     7,     8,     1,
       1,     1,     2,     1,     1
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
    case YYSYMBOL_AROMATIC_ATOM_TOKEN: /* AROMATIC_ATOM_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_ATOM_TOKEN: /* ATOM_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_ORGANIC_ATOM_TOKEN: /* ORGANIC_ATOM_TOKEN  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_BOND_TOKEN: /* BOND_TOKEN  */
            { delete ((*yyvaluep).bond); }
        break;

    case YYSYMBOL_bondd: /* bondd  */
            { delete ((*yyvaluep).bond); }
        break;

    case YYSYMBOL_atomd: /* atomd  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_charge_element: /* charge_element  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_h_element: /* h_element  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_chiral_element: /* chiral_element  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_element: /* element  */
            { delete ((*yyvaluep).atom); }
        break;

    case YYSYMBOL_simple_atom: /* simple_atom  */
            { delete ((*yyvaluep).atom); }
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

  case 5: /* meta_start: START_BOND bondd EOS_TOKEN  */
                             {
  lastBond = (yyvsp[-1].bond);
  YYACCEPT;
}
    break;

  case 6: /* meta_start: START_BOND bondd  */
                   {
  delete (yyvsp[0].bond);
  YYABORT;
}
    break;

  case 7: /* meta_start: START_BOND  */
             {
  YYABORT;
}
    break;

  case 8: /* meta_start: meta_start error EOS_TOKEN  */
                            {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
    break;

  case 9: /* meta_start: meta_start EOS_TOKEN  */
                       {
  YYACCEPT;
}
    break;

  case 10: /* meta_start: error EOS_TOKEN  */
                  {
  yyerrok;
  yyErrorCleanup(molList);
  YYABORT;
}
    break;

  case 14: /* bad_atom_def: charge_element  */
                 {
  delete (yyvsp[0].atom);
  YYABORT;
}
    break;

  case 15: /* mol: atomd  */
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
    break;

  case 16: /* mol: mol atomd  */
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
    break;

  case 17: /* mol: mol BOND_TOKEN atomd  */
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
    break;

  case 18: /* mol: mol MINUS_TOKEN atomd  */
                        {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $3;
}
    break;

  case 19: /* mol: mol SEPARATOR_TOKEN atomd  */
                            {
  RWMol *mp = (*molList)[(yyval.moli)];
  (yyvsp[0].atom)->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom((yyvsp[0].atom),true,true);
}
    break;

  case 20: /* mol: mol ring_number  */
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
    break;

  case 21: /* mol: mol BOND_TOKEN ring_number  */
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
    break;

  case 22: /* mol: mol MINUS_TOKEN ring_number  */
                              {
  RWMol * mp = (*molList)[(yyval.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::SINGLE);
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
    break;

  case 23: /* mol: mol branch_open_token atomd  */
                              {
  RWMol *mp = (*molList)[(yyval.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  branchPoints.push_back({atomIdx1, (yyvsp[-1].ival)});
}
    break;

  case 24: /* mol: mol branch_open_token BOND_TOKEN atomd  */
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

  case 25: /* mol: mol branch_open_token MINUS_TOKEN atomd  */
                                          {
  RWMol *mp = (*molList)[(yyval.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yyvsp[0].atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  branchPoints.push_back({atomIdx1, (yyvsp[-2].ival)});
}
    break;

  case 26: /* mol: mol GROUP_CLOSE_TOKEN  */
                        {
  if(branchPoints.empty()){
     yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"extra close parentheses");
     yyErrorCleanup(molList);
     YYABORT;
  }
  RWMol *mp = (*molList)[(yyval.moli)];

  mp->setActiveAtom(branchPoints.back().first);
  branchPoints.pop_back();
}
    break;

  case 27: /* branch_open_token: GROUP_OPEN_TOKEN  */
                                    { (yyval.ival) = current_token_position; }
    break;

  case 29: /* bondd: MINUS_TOKEN  */
                        {
          (yyval.bond) = new Bond(Bond::SINGLE);
          }
    break;

  case 31: /* atomd: ATOM_OPEN_TOKEN charge_element COLON_TOKEN number ATOM_CLOSE_TOKEN  */
{
  (yyval.atom) = (yyvsp[-3].atom);
  (yyval.atom)->setNoImplicit(true);
  (yyval.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yyvsp[-1].ival));
}
    break;

  case 32: /* atomd: ATOM_OPEN_TOKEN charge_element ATOM_CLOSE_TOKEN  */
{
  (yyval.atom) = (yyvsp[-1].atom);
  (yyvsp[-1].atom)->setNoImplicit(true);
}
    break;

  case 34: /* charge_element: h_element PLUS_TOKEN  */
                       { (yyvsp[-1].atom)->setFormalCharge(1); }
    break;

  case 35: /* charge_element: h_element PLUS_TOKEN PLUS_TOKEN  */
                                  { (yyvsp[-2].atom)->setFormalCharge(2); }
    break;

  case 36: /* charge_element: h_element PLUS_TOKEN number  */
                              { (yyvsp[-2].atom)->setFormalCharge((yyvsp[0].ival)); }
    break;

  case 37: /* charge_element: h_element MINUS_TOKEN  */
                        { (yyvsp[-1].atom)->setFormalCharge(-1); }
    break;

  case 38: /* charge_element: h_element MINUS_TOKEN MINUS_TOKEN  */
                                    { (yyvsp[-2].atom)->setFormalCharge(-2); }
    break;

  case 39: /* charge_element: h_element MINUS_TOKEN number  */
                               { (yyvsp[-2].atom)->setFormalCharge(-(yyvsp[0].ival)); }
    break;

  case 40: /* h_element: H_TOKEN  */
                        { (yyval.atom) = new Atom(1); }
    break;

  case 41: /* h_element: number H_TOKEN  */
                                 { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-1].ival)); }
    break;

  case 42: /* h_element: H_TOKEN H_TOKEN  */
                                  { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs(1); }
    break;

  case 43: /* h_element: number H_TOKEN H_TOKEN  */
                                         { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-2].ival)); (yyval.atom)->setNumExplicitHs(1);}
    break;

  case 44: /* h_element: H_TOKEN H_TOKEN number  */
                                         { (yyval.atom) = new Atom(1); (yyval.atom)->setNumExplicitHs((yyvsp[0].ival)); }
    break;

  case 45: /* h_element: number H_TOKEN H_TOKEN number  */
                                                { (yyval.atom) = new Atom(1); (yyval.atom)->setIsotope((yyvsp[-3].ival)); (yyval.atom)->setNumExplicitHs((yyvsp[0].ival));}
    break;

  case 47: /* h_element: chiral_element H_TOKEN  */
                                                        { (yyval.atom) = (yyvsp[-1].atom); (yyvsp[-1].atom)->setNumExplicitHs(1);}
    break;

  case 48: /* h_element: chiral_element H_TOKEN number  */
                                                { (yyval.atom) = (yyvsp[-2].atom); (yyvsp[-2].atom)->setNumExplicitHs((yyvsp[0].ival));}
    break;

  case 50: /* chiral_element: element AT_TOKEN  */
                   { (yyvsp[-1].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
    break;

  case 51: /* chiral_element: element AT_TOKEN AT_TOKEN  */
                            { (yyvsp[-2].atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
    break;

  case 52: /* chiral_element: element CHI_CLASS_TOKEN  */
                          { (yyvsp[-1].atom)->setChiralTag((yyvsp[0].chiraltype)); (yyvsp[-1].atom)->setProp(common_properties::_chiralPermutation,0); }
    break;

  case 53: /* chiral_element: element CHI_CLASS_TOKEN number  */
                                 { (yyvsp[-2].atom)->setChiralTag((yyvsp[-1].chiraltype)); (yyvsp[-2].atom)->setProp(common_properties::_chiralPermutation,(yyvsp[0].ival)); }
    break;

  case 55: /* element: number simple_atom  */
                                           { (yyvsp[0].atom)->setIsotope( (yyvsp[-1].ival) ); (yyval.atom) = (yyvsp[0].atom); }
    break;

  case 57: /* element: number ATOM_TOKEN  */
                                                   { (yyvsp[0].atom)->setIsotope( (yyvsp[-1].ival) ); (yyval.atom) = (yyvsp[0].atom); }
    break;

  case 58: /* element: HASH_TOKEN number  */
                                                 { (yyval.atom) = new Atom((yyvsp[0].ival)); }
    break;

  case 59: /* element: number HASH_TOKEN number  */
                                                         { (yyval.atom) = new Atom((yyvsp[0].ival)); (yyval.atom)->setIsotope((yyvsp[-2].ival)); }
    break;

  case 63: /* ring_number: PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit  */
                                          { (yyval.ival) = (yyvsp[-1].ival)*10+(yyvsp[0].ival); }
    break;

  case 64: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN  */
                                                         { (yyval.ival) = (yyvsp[-1].ival); }
    break;

  case 65: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN  */
                                                               { (yyval.ival) = (yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 66: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN  */
                                                                     { (yyval.ival) = (yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 67: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN  */
                                                                           { (yyval.ival) = (yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 68: /* ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN  */
                                                                                 { (yyval.ival) = (yyvsp[-5].ival)*10000+(yyvsp[-4].ival)*1000+(yyvsp[-3].ival)*100+(yyvsp[-2].ival)*10+(yyvsp[-1].ival); }
    break;

  case 72: /* nonzero_number: nonzero_number digit  */
                       {
  if((yyvsp[-1].ival) >= std::numeric_limits<std::int32_t>::max()/10 ||
     (yyvsp[-1].ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yyvsp[0].ival) ){
     yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"number too large");
     yyErrorCleanup(molList);
     YYABORT;
  }
  (yyval.ival) = (yyvsp[-1].ival)*10 + (yyvsp[0].ival);
  }
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


