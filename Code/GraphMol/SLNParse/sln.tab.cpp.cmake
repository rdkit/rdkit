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
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         yysln_parse
#define yylex           yysln_lex
#define yyerror         yysln_error
#define yydebug         yysln_debug
#define yynerrs         yysln_nerrs

/* First part of user prologue.  */
#line 1 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"


  // $Id$
  //
  //  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
  //  All rights reserved.
  //
  // Redistribution and use in source and binary forms, with or without
  // modification, are permitted provided that the following conditions are
  // met:
  //
  //     * Redistributions of source code must retain the above copyright
  //       notice, this list of conditions and the following disclaimer.
  //     * Redistributions in binary form must reproduce the above
  //       copyright notice, this list of conditions and the following
  //       disclaimer in the documentation and/or other materials provided
  //       with the distribution.
  //     * Neither the name of Novartis Institutes for BioMedical Research Inc.
  //       nor the names of its contributors may be used to endorse or promote
  //       products derived from this software without specific prior
  //       written permission.
  //
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  //
  // Created by Greg Landrum, September 2006
  //

#include <cstring>
#include <cstdio>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SLNParse/SLNParseOps.h>
#include <GraphMol/SLNParse/SLNAttribs.h>
#include <GraphMol/SLNParse/SLNParse.h>
#include <RDGeneral/RDLog.h>

#define YYDEBUG 1
#include "sln.tab.hpp"

int yysln_lex(YYSTYPE *,void *);

namespace SLNParse = RDKit::SLNParse;

void
yysln_error( const char *input,
             std::vector<RDKit::RWMol *> *ms,bool doQ,
	     void *scanner,const char * msg )
{
  RDUNUSED_PARAM(ms);
  RDUNUSED_PARAM(doQ);
  RDUNUSED_PARAM(scanner);
  BOOST_LOG(rdErrorLog)<<"SLN Parse Error: "<<msg<<" while parsing: "<<input<<std::endl;

  for(auto& m : *ms) {
    SLNParse::CleanupAfterParse(m);
    delete m;
  }
  ms->clear();
  ms->resize(0);
}

#define YYPRINT(file, type, value)   yyprint (file, type, value)

static void
yyprint (FILE *file, int type, YYSTYPE value)
{
  if (type == TEXT_BLOCK)
    fprintf (file, " %s", value.text_T->c_str());
  else fprintf (file, " %d", type);
}


#line 163 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"

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

#include "sln.tab.hpp"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_TEXT_BLOCK = 3,                 /* TEXT_BLOCK  */
  YYSYMBOL_CHAR_TOKEN = 4,                 /* CHAR_TOKEN  */
  YYSYMBOL_DIGIT_TOKEN = 5,                /* DIGIT_TOKEN  */
  YYSYMBOL_H_TOKEN = 6,                    /* H_TOKEN  */
  YYSYMBOL_H_BRACKET_TOKEN = 7,            /* H_BRACKET_TOKEN  */
  YYSYMBOL_H_ASTERIX_TOKEN = 8,            /* H_ASTERIX_TOKEN  */
  YYSYMBOL_AT_TOKEN = 9,                   /* AT_TOKEN  */
  YYSYMBOL_ATOM_TOKEN = 10,                /* ATOM_TOKEN  */
  YYSYMBOL_COMPARE_TOKEN = 11,             /* COMPARE_TOKEN  */
  YYSYMBOL_OPEN_PAREN_TOKEN = 12,          /* OPEN_PAREN_TOKEN  */
  YYSYMBOL_CLOSE_PAREN_TOKEN = 13,         /* CLOSE_PAREN_TOKEN  */
  YYSYMBOL_OPEN_BRACKET_TOKEN = 14,        /* OPEN_BRACKET_TOKEN  */
  YYSYMBOL_CLOSE_BRACKET_TOKEN = 15,       /* CLOSE_BRACKET_TOKEN  */
  YYSYMBOL_OPEN_ANGLE_TOKEN = 16,          /* OPEN_ANGLE_TOKEN  */
  YYSYMBOL_CLOSE_ANGLE_TOKEN = 17,         /* CLOSE_ANGLE_TOKEN  */
  YYSYMBOL_SEPARATOR_TOKEN = 18,           /* SEPARATOR_TOKEN  */
  YYSYMBOL_ASTERIX_TOKEN = 19,             /* ASTERIX_TOKEN  */
  YYSYMBOL_EOS_TOKEN = 20,                 /* EOS_TOKEN  */
  YYSYMBOL_PLUS_TOKEN = 21,                /* PLUS_TOKEN  */
  YYSYMBOL_MINUS_TOKEN = 22,               /* MINUS_TOKEN  */
  YYSYMBOL_COLON_TOKEN = 23,               /* COLON_TOKEN  */
  YYSYMBOL_EQUALS_TOKEN = 24,              /* EQUALS_TOKEN  */
  YYSYMBOL_TILDE_TOKEN = 25,               /* TILDE_TOKEN  */
  YYSYMBOL_HASH_TOKEN = 26,                /* HASH_TOKEN  */
  YYSYMBOL_COMMA_TOKEN = 27,               /* COMMA_TOKEN  */
  YYSYMBOL_NOT_TOKEN = 28,                 /* NOT_TOKEN  */
  YYSYMBOL_AND_TOKEN = 29,                 /* AND_TOKEN  */
  YYSYMBOL_OR_TOKEN = 30,                  /* OR_TOKEN  */
  YYSYMBOL_SEMI_TOKEN = 31,                /* SEMI_TOKEN  */
  YYSYMBOL_CARET_EQUALS_TOKEN = 32,        /* CARET_EQUALS_TOKEN  */
  YYSYMBOL_COLON_EQUALS_TOKEN = 33,        /* COLON_EQUALS_TOKEN  */
  YYSYMBOL_RECURSE_TOKEN = 34,             /* RECURSE_TOKEN  */
  YYSYMBOL_NEG_RECURSE_TOKEN = 35,         /* NEG_RECURSE_TOKEN  */
  YYSYMBOL_ERROR_TOKEN = 36,               /* ERROR_TOKEN  */
  YYSYMBOL_YYACCEPT = 37,                  /* $accept  */
  YYSYMBOL_cmpd = 38,                      /* cmpd  */
  YYSYMBOL_mol = 39,                       /* mol  */
  YYSYMBOL_primmol = 40,                   /* primmol  */
  YYSYMBOL_atom = 41,                      /* atom  */
  YYSYMBOL_hatom = 42,                     /* hatom  */
  YYSYMBOL_primatom = 43,                  /* primatom  */
  YYSYMBOL_bond = 44,                      /* bond  */
  YYSYMBOL_primbond = 45,                  /* primbond  */
  YYSYMBOL_onebond = 46,                   /* onebond  */
  YYSYMBOL_attriblist = 47,                /* attriblist  */
  YYSYMBOL_ctabattriblist = 48,            /* ctabattriblist  */
  YYSYMBOL_attrib = 49,                    /* attrib  */
  YYSYMBOL_recursivequery = 50,            /* recursivequery  */
  YYSYMBOL_ctabattrib = 51,                /* ctabattrib  */
  YYSYMBOL_number = 52                     /* number  */
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
#define YYFINAL  25
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   219

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  16
/* YYNRULES -- Number of rules.  */
#define YYNRULES  72
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  116

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   291


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
      35,    36
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   147,   147,   148,   156,   161,   164,   171,   172,   177,
     187,   190,   194,   198,   202,   207,   211,   217,   221,   225,
     236,   248,   261,   273,   287,   299,   313,   314,   315,   319,
     325,   333,   341,   351,   362,   363,   366,   370,   376,   384,
     385,   392,   397,   406,   407,   420,   429,   438,   448,   460,
     464,   468,   472,   479,   484,   490,   498,   502,   513,   519,
     525,   531,   537,   543,   548,   562,   578,   592,   600,   610,
     620,   632,   633
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
  "\"end of file\"", "error", "\"invalid token\"", "TEXT_BLOCK",
  "CHAR_TOKEN", "DIGIT_TOKEN", "H_TOKEN", "H_BRACKET_TOKEN",
  "H_ASTERIX_TOKEN", "AT_TOKEN", "ATOM_TOKEN", "COMPARE_TOKEN",
  "OPEN_PAREN_TOKEN", "CLOSE_PAREN_TOKEN", "OPEN_BRACKET_TOKEN",
  "CLOSE_BRACKET_TOKEN", "OPEN_ANGLE_TOKEN", "CLOSE_ANGLE_TOKEN",
  "SEPARATOR_TOKEN", "ASTERIX_TOKEN", "EOS_TOKEN", "PLUS_TOKEN",
  "MINUS_TOKEN", "COLON_TOKEN", "EQUALS_TOKEN", "TILDE_TOKEN",
  "HASH_TOKEN", "COMMA_TOKEN", "NOT_TOKEN", "AND_TOKEN", "OR_TOKEN",
  "SEMI_TOKEN", "CARET_EQUALS_TOKEN", "COLON_EQUALS_TOKEN",
  "RECURSE_TOKEN", "NEG_RECURSE_TOKEN", "ERROR_TOKEN", "$accept", "cmpd",
  "mol", "primmol", "atom", "hatom", "primatom", "bond", "primbond",
  "onebond", "attriblist", "ctabattriblist", "attrib", "recursivequery",
  "ctabattrib", "number", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-40)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-67)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     194,   -13,   -40,    28,   -40,   -40,     5,    -6,   154,   -40,
     -40,    46,   -40,     3,   -40,   -40,    12,    15,    45,   194,
     194,     7,   -40,    17,     4,   -40,    26,    48,   -40,   209,
      50,   107,   -40,   -40,   -40,    39,   -40,   -40,   199,   170,
     -40,    54,    28,   -40,    65,   -40,   -40,   -40,    69,    75,
     -40,    45,    45,    45,   194,   -40,   -40,    45,   -40,   141,
      -2,   -40,   154,    18,    50,   127,   204,    45,   -40,    50,
     -40,    45,   -40,   -40,    57,    11,   -40,   -40,   -40,   -40,
      92,    95,    68,    71,    72,   -40,    48,   -40,    64,    86,
      50,   134,   139,    76,   152,   -40,   -40,    45,   -40,   -40,
     -40,   -40,   -40,    88,   -40,    96,    91,   -40,   -40,   -40,
     156,   -40,    97,   -40,   -40,   -40
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,     9,     0,    30,    34,     0,     2,     7,    10,
      26,    27,     6,    55,    71,    62,    58,    60,     0,     0,
       0,     0,    52,    63,     0,     1,     0,     0,     5,     0,
       0,     0,    45,    48,    46,    41,    47,    11,     0,    39,
      43,    28,     0,    35,     0,    59,    61,    56,     0,     0,
      33,     0,     0,     0,     0,    72,    31,     0,     4,    67,
       0,    53,     8,    13,     0,     0,     0,     0,    19,     0,
      12,     0,    44,    29,     0,     0,    57,    49,    50,    51,
       0,     0,     0,     0,     0,     3,     0,    20,     0,    15,
       0,     0,     0,    14,     0,    38,    36,     0,    32,    68,
      70,    69,    54,    17,    22,     0,    16,    42,    21,    40,
       0,    24,    18,    23,    37,    25
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -40,    -9,   -40,   -27,    -8,   -40,   -40,    47,   -40,    79,
     -39,   -40,   -10,   -40,    25,   -29
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     6,     7,     8,     9,    10,    11,    38,    39,    40,
      21,    60,    22,    23,    61,    24
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      37,    63,    62,    74,    65,    25,    26,    12,    47,    55,
      48,    49,    29,    75,    44,    85,    55,    45,    81,    56,
      46,    27,    50,    55,    87,    28,    96,    57,    92,    86,
      70,    13,    94,    14,    97,    88,    51,    52,    53,    91,
      93,    77,    78,    79,    54,    80,    58,    15,    13,    16,
      17,    59,    41,    67,    37,    14,    18,    37,   110,    73,
      42,   105,    19,    20,    15,    43,    16,    17,    76,    55,
      26,    99,    95,    18,   100,   101,    26,   103,    66,    19,
      20,    55,   108,    37,   -64,    27,    51,    52,    53,    28,
     -65,    27,   104,    26,   111,    28,   -64,   113,   -64,   -64,
     -64,    55,   -65,   115,   -65,   -65,   -65,   -66,    27,   112,
      98,   102,    28,     2,     3,     4,    64,     5,    72,   -66,
       0,   -66,   -66,   -66,    51,    52,    53,     0,     0,    32,
      33,    34,    35,    36,     3,     4,    30,     5,     0,    31,
      89,     3,     4,    30,     5,     0,    31,   106,     0,    32,
      33,    34,    35,    36,   107,     0,    32,    33,    34,    35,
      36,     3,     4,    30,     5,    82,    31,   109,    51,    52,
      53,   114,     0,    83,    84,     0,    32,    33,    34,    35,
      36,    51,    52,    53,    71,    51,    52,    53,     0,     0,
       0,     0,    32,    33,    34,     1,    36,     0,     0,     0,
       2,     3,     4,     0,     5,    68,     3,     4,    69,     5,
       2,     3,     4,    90,     5,     2,     3,     4,     0,     5
};

static const yytype_int8 yycheck[] =
{
       8,    30,    29,    42,    31,     0,     1,    20,    18,     5,
      19,    20,    18,    42,    11,    17,     5,     5,    57,    15,
       5,    16,    15,     5,     6,    20,    15,    23,    67,    31,
      38,     3,    71,     5,    23,    64,    29,    30,    31,    66,
      69,    51,    52,    53,    27,    54,    20,    19,     3,    21,
      22,     3,     6,    14,    62,     5,    28,    65,    97,     5,
      14,    90,    34,    35,    19,    19,    21,    22,     3,     5,
       1,     3,    15,    28,     3,     3,     1,    13,    31,    34,
      35,     5,     6,    91,    15,    16,    29,    30,    31,    20,
      15,    16,     6,     1,     6,    20,    27,     6,    29,    30,
      31,     5,    27,     6,    29,    30,    31,    15,    16,    13,
      15,    86,    20,     6,     7,     8,     9,    10,    39,    27,
      -1,    29,    30,    31,    29,    30,    31,    -1,    -1,    22,
      23,    24,    25,    26,     7,     8,     9,    10,    -1,    12,
      13,     7,     8,     9,    10,    -1,    12,    13,    -1,    22,
      23,    24,    25,    26,    15,    -1,    22,    23,    24,    25,
      26,     7,     8,     9,    10,    24,    12,    15,    29,    30,
      31,    15,    -1,    32,    33,    -1,    22,    23,    24,    25,
      26,    29,    30,    31,    14,    29,    30,    31,    -1,    -1,
      -1,    -1,    22,    23,    24,     1,    26,    -1,    -1,    -1,
       6,     7,     8,    -1,    10,     6,     7,     8,     9,    10,
       6,     7,     8,     9,    10,     6,     7,     8,    -1,    10
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     6,     7,     8,    10,    38,    39,    40,    41,
      42,    43,    20,     3,     5,    19,    21,    22,    28,    34,
      35,    47,    49,    50,    52,     0,     1,    16,    20,    18,
       9,    12,    22,    23,    24,    25,    26,    41,    44,    45,
      46,     6,    14,    19,    11,     5,     5,    49,    38,    38,
      15,    29,    30,    31,    27,     5,    15,    23,    20,     3,
      48,    51,    40,    52,     9,    40,    44,    14,     6,     9,
      41,    14,    46,     5,    47,    52,     3,    49,    49,    49,
      38,    47,    24,    32,    33,    17,    31,     6,    52,    13,
       9,    40,    47,    52,    47,    15,    15,    23,    15,     3,
       3,     3,    51,    13,     6,    52,    13,    15,     6,    15,
      47,     6,    13,     6,    15,     6
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    37,    38,    38,    38,    38,    38,    39,    39,    40,
      40,    40,    40,    40,    40,    40,    40,    40,    40,    40,
      40,    40,    40,    40,    40,    40,    41,    41,    41,    41,
      42,    42,    42,    42,    43,    43,    43,    43,    43,    44,
      44,    44,    44,    45,    45,    46,    46,    46,    46,    47,
      47,    47,    47,    48,    48,    49,    49,    49,    49,    49,
      49,    49,    49,    49,    50,    50,    50,    51,    51,    51,
      51,    52,    52
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     4,     3,     2,     2,     1,     3,     1,
       1,     2,     3,     3,     4,     4,     5,     5,     6,     3,
       4,     5,     5,     6,     6,     7,     1,     1,     2,     3,
       1,     3,     5,     3,     1,     2,     4,     6,     4,     1,
       4,     1,     4,     1,     2,     1,     1,     1,     1,     3,
       3,     3,     1,     1,     3,     1,     2,     3,     1,     2,
       1,     2,     1,     1,     2,     2,     3,     1,     3,     3,
       3,     1,     2
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
        yyerror (input, molList, doQueries, scanner, YY_("syntax error: cannot back up")); \
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
                  Kind, Value, input, molList, doQueries, scanner); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  YY_USE (input);
  YY_USE (molList);
  YY_USE (doQueries);
  YY_USE (scanner);
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
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep, input, molList, doQueries, scanner);
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
                 int yyrule, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
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
                       &yyvsp[(yyi + 1) - (yynrhs)], input, molList, doQueries, scanner);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, input, molList, doQueries, scanner); \
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
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep, const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
{
  YY_USE (yyvaluep);
  YY_USE (input);
  YY_USE (molList);
  YY_USE (doQueries);
  YY_USE (scanner);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yykind)
    {
    case YYSYMBOL_TEXT_BLOCK: /* TEXT_BLOCK  */
#line 139 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).text_T); }
#line 1049 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_ATOM_TOKEN: /* ATOM_TOKEN  */
#line 137 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).atom_T); }
#line 1055 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_COMPARE_TOKEN: /* COMPARE_TOKEN  */
#line 139 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).text_T); }
#line 1061 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_atom: /* atom  */
#line 137 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).atom_T); }
#line 1067 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_hatom: /* hatom  */
#line 137 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).atom_T); }
#line 1073 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_primatom: /* primatom  */
#line 137 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).atom_T); }
#line 1079 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_bond: /* bond  */
#line 138 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).bond_T); }
#line 1085 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_primbond: /* primbond  */
#line 138 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).bond_T); }
#line 1091 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_onebond: /* onebond  */
#line 138 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).bond_T); }
#line 1097 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_attriblist: /* attriblist  */
#line 141 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).attriblist_T); }
#line 1103 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_ctabattriblist: /* ctabattriblist  */
#line 141 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).attriblist_T); }
#line 1109 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_attrib: /* attrib  */
#line 140 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).attrib_T); }
#line 1115 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_recursivequery: /* recursivequery  */
#line 140 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).attrib_T); }
#line 1121 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
        break;

    case YYSYMBOL_ctabattrib: /* ctabattrib  */
#line 140 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
            { delete ((*yyvaluep).attrib_T); }
#line 1127 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
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
yyparse (const char *input, std::vector<RDKit::RWMol *> *molList, bool doQueries, void *scanner)
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
      yychar = yylex (&yylval, scanner);
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
  case 3: /* cmpd: cmpd OPEN_ANGLE_TOKEN ctabattriblist CLOSE_ANGLE_TOKEN  */
#line 148 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                         {
  // allowing mol<attrs><attrs> seems to be a NIBR thing, I don't
  // think it's standard SLN
  RDKit::ROMol *mol=(*molList)[(yyvsp[-3].mol_T)];
  SLNParse::parseMolAttribs(mol,*(yyvsp[-1].attriblist_T));
  delete (yyvsp[-1].attriblist_T);
  (yyval.mol_T)=(yyvsp[-3].mol_T);
}
#line 1410 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 4: /* cmpd: cmpd error EOS_TOKEN  */
#line 156 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                       {
  yyclearin;
  yyerrok;
  YYABORT;
}
#line 1420 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 5: /* cmpd: cmpd EOS_TOKEN  */
#line 161 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                 {
  YYACCEPT;
}
#line 1428 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 6: /* cmpd: error EOS_TOKEN  */
#line 164 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                  {
  yyclearin;
  yyerrok;
  YYABORT;
}
#line 1438 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 8: /* mol: mol SEPARATOR_TOKEN primmol  */
#line 172 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                              {
  (yyval.mol_T)=SLNParse::addFragToMol(*molList,(yyvsp[-2].mol_T),(yyvsp[0].mol_T));
}
#line 1446 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 9: /* primmol: H_TOKEN  */
#line 177 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                 {
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }

  (yyval.mol_T)=SLNParse::startMol(*molList,newAtom,doQueries);
}
#line 1461 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 10: /* primmol: atom  */
#line 187 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
       {
  (yyval.mol_T)=SLNParse::startMol(*molList,(yyvsp[0].atom_T),doQueries);
}
#line 1469 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 11: /* primmol: primmol atom  */
#line 190 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
               {
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),(yyvsp[0].atom_T),doQueries);
  (yyval.mol_T)=(yyvsp[-1].mol_T);
}
#line 1478 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 12: /* primmol: primmol bond atom  */
#line 194 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                    {
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),(yyvsp[0].atom_T),(yyvsp[-1].bond_T),doQueries);
  (yyval.mol_T)=(yyvsp[-2].mol_T);
}
#line 1487 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 13: /* primmol: primmol AT_TOKEN number  */
#line 198 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                          {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[0].ival_T));
  (yyval.mol_T)=(yyvsp[-2].mol_T);
}
#line 1496 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 14: /* primmol: primmol bond AT_TOKEN number  */
#line 202 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                               {
  // closeRingBond() takes ownership of the bond
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[0].ival_T),(yyvsp[-2].bond_T));
  (yyval.mol_T)=(yyvsp[-3].mol_T);
}
#line 1506 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 15: /* primmol: primmol OPEN_PAREN_TOKEN primmol CLOSE_PAREN_TOKEN  */
#line 207 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                     {
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[-1].mol_T));
  (yyval.mol_T)=(yyvsp[-3].mol_T);
}
#line 1515 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 16: /* primmol: primmol OPEN_PAREN_TOKEN bond primmol CLOSE_PAREN_TOKEN  */
#line 211 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                         {
  // addBranchToMol() takes ownership of the bond and deletes the
  // branch, so no leaks here'
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[-1].mol_T),(yyvsp[-2].bond_T));
  (yyval.mol_T)=(yyvsp[-4].mol_T);
}
#line 1526 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 17: /* primmol: primmol OPEN_PAREN_TOKEN AT_TOKEN number CLOSE_PAREN_TOKEN  */
#line 217 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                            {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[-1].ival_T));
  (yyval.mol_T)=(yyvsp[-4].mol_T);
}
#line 1535 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 18: /* primmol: primmol OPEN_PAREN_TOKEN bond AT_TOKEN number CLOSE_PAREN_TOKEN  */
#line 221 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                                 {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[-1].ival_T),(yyvsp[-3].bond_T));
  (yyval.mol_T)=(yyvsp[-5].mol_T);
}
#line 1544 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 19: /* primmol: primmol bond H_TOKEN  */
#line 225 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                       {
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }

  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,(yyvsp[-1].bond_T),doQueries);
  (yyval.mol_T)=(yyvsp[-2].mol_T);
}
#line 1560 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 20: /* primmol: primmol AT_TOKEN number H_TOKEN  */
#line 236 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                  {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[-1].ival_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);

  (yyval.mol_T)=(yyvsp[-3].mol_T);
}
#line 1577 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 21: /* primmol: primmol bond AT_TOKEN number H_TOKEN  */
#line 248 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                       {
  // closeRingBond() takes ownership of the bond
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[-1].ival_T),(yyvsp[-3].bond_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);

  (yyval.mol_T)=(yyvsp[-4].mol_T);
}
#line 1595 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 22: /* primmol: primmol OPEN_PAREN_TOKEN primmol CLOSE_PAREN_TOKEN H_TOKEN  */
#line 261 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                             {
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[-2].mol_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);

  (yyval.mol_T)=(yyvsp[-4].mol_T);
}
#line 1612 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 23: /* primmol: primmol OPEN_PAREN_TOKEN bond primmol CLOSE_PAREN_TOKEN H_TOKEN  */
#line 273 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                                  {
  // addBranchToMol() takes ownership of the bond and deletes the
  // branch, so no leaks here'
  SLNParse::addBranchToMol(*molList,(yyval.mol_T),(yyvsp[-2].mol_T),(yyvsp[-3].bond_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);

  (yyval.mol_T)=(yyvsp[-5].mol_T);
}
#line 1631 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 24: /* primmol: primmol OPEN_PAREN_TOKEN AT_TOKEN number CLOSE_PAREN_TOKEN H_TOKEN  */
#line 287 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                                     {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[-2].ival_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);

  (yyval.mol_T)=(yyvsp[-5].mol_T);
}
#line 1648 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 25: /* primmol: primmol OPEN_PAREN_TOKEN bond AT_TOKEN number CLOSE_PAREN_TOKEN H_TOKEN  */
#line 299 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                                          {
  SLNParse::closeRingBond(*molList,(yyval.mol_T),(yyvsp[-2].ival_T),(yyvsp[-4].bond_T));
  RDKit::Atom *newAtom;
  if(!doQueries){
    newAtom = new RDKit::Atom(1);
  } else {
    newAtom = new RDKit::QueryAtom(1);
  }
  SLNParse::addAtomToMol(*molList,(yyval.mol_T),newAtom,doQueries);

  (yyval.mol_T)=(yyvsp[-6].mol_T);
}
#line 1665 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 28: /* atom: primatom H_TOKEN  */
#line 315 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                   {
  (yyvsp[-1].atom_T)->setNumSpecifiedHs(1);
  (yyval.atom_T)=(yyvsp[-1].atom_T);
}
#line 1674 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 29: /* atom: primatom H_TOKEN DIGIT_TOKEN  */
#line 319 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                               {
  (yyvsp[-2].atom_T)->setNumSpecifiedHs((yyvsp[0].ival_T));
  (yyval.atom_T)=(yyvsp[-2].atom_T);
}
#line 1683 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 30: /* hatom: H_ASTERIX_TOKEN  */
#line 325 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                       {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  (yyval.atom_T)->setProp(RDKit::common_properties::_starred,1,true);
}
#line 1696 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 31: /* hatom: H_BRACKET_TOKEN number CLOSE_BRACKET_TOKEN  */
#line 333 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                              {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  (yyval.atom_T)->setProp(RDKit::common_properties::_AtomID,static_cast<unsigned int>((yyvsp[-1].ival_T)));
}
#line 1709 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 32: /* hatom: H_BRACKET_TOKEN number COLON_TOKEN attriblist CLOSE_BRACKET_TOKEN  */
#line 341 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                                    {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  (yyval.atom_T)->setProp(RDKit::common_properties::_AtomID,static_cast<unsigned int>((yyvsp[-3].ival_T)));
  SLNParse::parseAtomAttribs((yyval.atom_T),*(yyvsp[-1].attriblist_T),doQueries);
  delete (yyvsp[-1].attriblist_T);
}
#line 1724 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 33: /* hatom: H_BRACKET_TOKEN attriblist CLOSE_BRACKET_TOKEN  */
#line 351 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                 {
  if(!doQueries){
    (yyval.atom_T) = new RDKit::Atom(1);
  } else {
    (yyval.atom_T) = new RDKit::QueryAtom(1);
  }
  SLNParse::parseAtomAttribs((yyval.atom_T),*(yyvsp[-1].attriblist_T),doQueries);
  delete (yyvsp[-1].attriblist_T);
}
#line 1738 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 35: /* primatom: primatom ASTERIX_TOKEN  */
#line 363 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                        {
  (yyval.atom_T)->setProp(RDKit::common_properties::_starred,1,true);
}
#line 1746 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 36: /* primatom: primatom OPEN_BRACKET_TOKEN number CLOSE_BRACKET_TOKEN  */
#line 366 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                         {
  (yyvsp[-3].atom_T)->setProp(RDKit::common_properties::_AtomID,static_cast<unsigned int>((yyvsp[-1].ival_T)));
  (yyval.atom_T)=(yyvsp[-3].atom_T);
}
#line 1755 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 37: /* primatom: primatom OPEN_BRACKET_TOKEN number COLON_TOKEN attriblist CLOSE_BRACKET_TOKEN  */
#line 370 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                                                {
  (yyvsp[-5].atom_T)->setProp(RDKit::common_properties::_AtomID,static_cast<unsigned int>((yyvsp[-3].ival_T)));
  SLNParse::parseAtomAttribs((yyvsp[-5].atom_T),*(yyvsp[-1].attriblist_T),doQueries);
  delete (yyvsp[-1].attriblist_T);
  (yyval.atom_T)=(yyvsp[-5].atom_T);
}
#line 1766 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 38: /* primatom: primatom OPEN_BRACKET_TOKEN attriblist CLOSE_BRACKET_TOKEN  */
#line 376 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                             {
  SLNParse::parseAtomAttribs((yyvsp[-3].atom_T),*(yyvsp[-1].attriblist_T),doQueries);
  delete (yyvsp[-1].attriblist_T);
  (yyval.atom_T)=(yyvsp[-3].atom_T);
}
#line 1776 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 40: /* bond: primbond OPEN_BRACKET_TOKEN attriblist CLOSE_BRACKET_TOKEN  */
#line 385 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                             {
  SLNParse::parseBondAttribs((yyvsp[-3].bond_T),*(yyvsp[-1].attriblist_T),doQueries);
  delete (yyvsp[-1].attriblist_T);
  (yyval.bond_T) = (yyvsp[-3].bond_T);
}
#line 1786 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 41: /* bond: TILDE_TOKEN  */
#line 392 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
              {
  RDKit::Bond *bond=new RDKit::QueryBond();
  bond->setQuery(RDKit::makeBondNullQuery());
  (yyval.bond_T) = bond;
}
#line 1796 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 42: /* bond: TILDE_TOKEN OPEN_BRACKET_TOKEN attriblist CLOSE_BRACKET_TOKEN  */
#line 397 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                                                {
  RDKit::Bond *bond=new RDKit::QueryBond();
  bond->setQuery(RDKit::makeBondNullQuery());
  SLNParse::parseBondAttribs(bond,*(yyvsp[-1].attriblist_T),doQueries);
  delete (yyvsp[-1].attriblist_T);
  (yyval.bond_T) = bond;
}
#line 1808 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 44: /* primbond: primbond onebond  */
#line 407 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                   {
	if(!doQueries){
        yysln_error(input,molList,doQueries,0,"sequential bonds not allowed in non-queries");
    YYABORT;
	} else {
	  RDKit::QueryBond *b1=static_cast<RDKit::QueryBond *>((yyvsp[-1].bond_T));
	  RDKit::QueryBond *b2=static_cast<RDKit::QueryBond *>((yyvsp[0].bond_T));
	  b1->expandQuery(b2->getQuery()->copy(),Queries::COMPOSITE_OR,true);
		delete b2;
	}
}
#line 1824 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 45: /* onebond: MINUS_TOKEN  */
#line 420 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                     {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::SINGLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::SINGLE);
  }
  (yyval.bond_T) = bond;
}
#line 1838 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 46: /* onebond: EQUALS_TOKEN  */
#line 429 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
               {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::DOUBLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::DOUBLE);
  }
  (yyval.bond_T) = bond;
}
#line 1852 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 47: /* onebond: HASH_TOKEN  */
#line 438 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
             {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::TRIPLE);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::TRIPLE);
  }
  (yyval.bond_T) = bond;

}
#line 1867 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 48: /* onebond: COLON_TOKEN  */
#line 448 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
              {
  RDKit::Bond *bond;
  if(doQueries){
    bond= new RDKit::QueryBond(RDKit::Bond::AROMATIC);
  } else {
    bond= new RDKit::Bond(RDKit::Bond::AROMATIC);
  }
  (yyval.bond_T) = bond;
}
#line 1881 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 49: /* attriblist: attriblist AND_TOKEN attrib  */
#line 460 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                       {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
}
#line 1890 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 50: /* attriblist: attriblist OR_TOKEN attrib  */
#line 464 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                            {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribOr,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
}
#line 1899 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 51: /* attriblist: attriblist SEMI_TOKEN attrib  */
#line 468 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                              {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribLowPriAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
}
#line 1908 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 52: /* attriblist: attrib  */
#line 472 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
         {
  (yyval.attriblist_T) = new SLNParse::AttribListType();
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribLowPriAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
}
#line 1918 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 53: /* ctabattriblist: ctabattrib  */
#line 479 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                           {
  (yyval.attriblist_T) = new SLNParse::AttribListType();
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
}
#line 1928 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 54: /* ctabattriblist: ctabattriblist SEMI_TOKEN ctabattrib  */
#line 484 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                       {
  (yyval.attriblist_T)->push_back(std::make_pair(SLNParse::AttribAnd,
                               boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
}
#line 1937 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 55: /* attrib: TEXT_BLOCK  */
#line 490 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                   {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[0].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  (yyval.attrib_T)->op = "";
  (yyval.attrib_T)->second = "";
  delete (yyvsp[0].text_T);
}
#line 1950 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 56: /* attrib: NOT_TOKEN attrib  */
#line 498 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                   {
  (yyvsp[0].attrib_T)->negated=true;
  (yyval.attrib_T)=(yyvsp[0].attrib_T);
}
#line 1959 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 57: /* attrib: TEXT_BLOCK COMPARE_TOKEN TEXT_BLOCK  */
#line 502 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                      {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[-2].text_T);
  (yyval.attrib_T)->op = *(yyvsp[-1].text_T);
  (yyval.attrib_T)->second = *(yyvsp[0].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[-2].text_T);
  delete (yyvsp[-1].text_T);
  delete (yyvsp[0].text_T);
}
#line 1975 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 58: /* attrib: PLUS_TOKEN  */
#line 513 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
             {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = "+1";
}
#line 1986 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 59: /* attrib: PLUS_TOKEN DIGIT_TOKEN  */
#line 519 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                         {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = SLNParse::convertToString((yyvsp[0].ival_T));
}
#line 1997 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 60: /* attrib: MINUS_TOKEN  */
#line 525 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
              {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = "-1";
}
#line 2008 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 61: /* attrib: MINUS_TOKEN DIGIT_TOKEN  */
#line 531 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                          {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "charge";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = SLNParse::convertToString(-(yyvsp[0].ival_T));
}
#line 2019 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 62: /* attrib: ASTERIX_TOKEN  */
#line 537 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = "spin";
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = "d";
}
#line 2030 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 63: /* attrib: recursivequery  */
#line 543 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                 {
  (yyval.attrib_T) = (yyvsp[0].attrib_T);
}
#line 2038 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 64: /* recursivequery: RECURSE_TOKEN cmpd  */
#line 548 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                   {
   int sz = molList->size();
   RDKit::ROMol *mol=(*molList)[(yyvsp[0].mol_T)];
   molList->resize( sz-1 );
   SLNParse::finalizeQueryMol(mol,true);
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(mol);
   RDKit::ATOM_OR_QUERY *orq=new RDKit::ATOM_OR_QUERY();
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   (yyval.attrib_T) = new SLNParse::AttribType();
   (yyval.attrib_T)->first="is";
   (yyval.attrib_T)->op = "=";
   (yyval.attrib_T)->second = "";
   (yyval.attrib_T)->structQuery=static_cast<void *>(orq);
}
#line 2057 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 65: /* recursivequery: NEG_RECURSE_TOKEN cmpd  */
#line 562 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                         {
   int sz = molList->size();
   RDKit::ROMol *mol=(*molList)[(yyvsp[0].mol_T)];
   molList->resize( sz-1 );
   SLNParse::finalizeQueryMol(mol,true);
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(mol);
   RDKit::ATOM_OR_QUERY *orq=new RDKit::ATOM_OR_QUERY();
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   orq->setNegation(true);

   (yyval.attrib_T) = new SLNParse::AttribType();
   (yyval.attrib_T)->first="is";
   (yyval.attrib_T)->op = "=";
   (yyval.attrib_T)->second = "";
   (yyval.attrib_T)->structQuery=static_cast<void *>(orq);
}
#line 2078 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 66: /* recursivequery: recursivequery COMMA_TOKEN cmpd  */
#line 578 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                  {
   int sz = molList->size();
   RDKit::ROMol *mol=(*molList)[(yyvsp[0].mol_T)];
   molList->resize( sz-1 );
   SLNParse::finalizeQueryMol(mol,true);
   RDKit::RecursiveStructureQuery *rsq=new RDKit::RecursiveStructureQuery(mol);

   RDKit::ATOM_OR_QUERY *orq=static_cast<RDKit::ATOM_OR_QUERY *>((yyvsp[-2].attrib_T)->structQuery);
   orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
   (yyval.attrib_T)=(yyvsp[-2].attrib_T);
}
#line 2094 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 67: /* ctabattrib: TEXT_BLOCK  */
#line 592 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                       {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[0].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  (yyval.attrib_T)->op = "";
  (yyval.attrib_T)->second = "";
  delete (yyvsp[0].text_T);
}
#line 2107 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 68: /* ctabattrib: TEXT_BLOCK EQUALS_TOKEN TEXT_BLOCK  */
#line 600 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                     {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[-2].text_T);
  (yyval.attrib_T)->op = "=";
  (yyval.attrib_T)->second = *(yyvsp[0].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[-2].text_T);
  delete (yyvsp[0].text_T);
}
#line 2122 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 69: /* ctabattrib: TEXT_BLOCK COLON_EQUALS_TOKEN TEXT_BLOCK  */
#line 610 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                           {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[-2].text_T);
  (yyval.attrib_T)->op = ":=";
  (yyval.attrib_T)->second = *(yyvsp[0].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[-2].text_T);
  delete (yyvsp[0].text_T);
}
#line 2137 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 70: /* ctabattrib: TEXT_BLOCK CARET_EQUALS_TOKEN TEXT_BLOCK  */
#line 620 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                                           {
  (yyval.attrib_T) = new SLNParse::AttribType();
  (yyval.attrib_T)->first = *(yyvsp[-2].text_T);
  (yyval.attrib_T)->op = "^=";
  (yyval.attrib_T)->second = *(yyvsp[0].text_T);
  boost::to_lower((yyval.attrib_T)->first);
  boost::to_lower((yyval.attrib_T)->second);
  delete (yyvsp[-2].text_T);
  delete (yyvsp[0].text_T);
}
#line 2152 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;

  case 72: /* number: number DIGIT_TOKEN  */
#line 633 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"
                     { (yyval.ival_T) = (yyvsp[-1].ival_T)*10 + (yyvsp[0].ival_T); }
#line 2158 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"
    break;


#line 2162 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.tab.cpp"

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
      yyerror (input, molList, doQueries, scanner, YY_("syntax error"));
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
                      yytoken, &yylval, input, molList, doQueries, scanner);
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
                  YY_ACCESSING_SYMBOL (yystate), yyvsp, input, molList, doQueries, scanner);
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
  yyerror (input, molList, doQueries, scanner, YY_("memory exhausted"));
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
                  yytoken, &yylval, input, molList, doQueries, scanner);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp, input, molList, doQueries, scanner);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 635 "/scratch/RDKit_git/Code/GraphMol/SLNParse/sln.yy"



