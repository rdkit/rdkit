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
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Substitute the variable and function names.  */
#define yyparse yysln_parse
#define yylex yysln_lex
#define yyerror yysln_error
#define yydebug yysln_debug
#define yynerrs yysln_nerrs

/* Copy the first part of user declarations.  */
#line 1 "sln.yy" /* yacc.c:339  */

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

int yysln_lex(YYSTYPE *, void *);

namespace SLNParse = RDKit::SLNParse;

void yysln_error(const char *input, std::vector<RDKit::RWMol *> *ms, bool doQ,
                 void *scanner, const char *msg) {
  RDUNUSED_PARAM(ms);
  RDUNUSED_PARAM(doQ);
  RDUNUSED_PARAM(scanner);
  BOOST_LOG(rdErrorLog) << "SLN Parse Error: " << msg
                        << " while parsing: " << input << std::endl;

  for (auto &m : *ms) {
    SLNParse::CleanupAfterParse(m);
    delete m;
  }
  ms->clear();
  ms->resize(0);
}

#define YYPRINT(file, type, value) yyprint(file, type, value)

static void yyprint(FILE *file, int type, YYSTYPE value) {
  if (type == TEXT_BLOCK)
    fprintf(file, " %s", value.text_T->c_str());
  else
    fprintf(file, " %d", type);
}

#line 156 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:339  */

#ifndef YY_NULLPTR
#if defined __cplusplus && 201103L <= __cplusplus
#define YY_NULLPTR nullptr
#else
#define YY_NULLPTR 0
#endif
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
#undef YYERROR_VERBOSE
#define YYERROR_VERBOSE 1
#else
#define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "sln.tab.hpp".  */
#ifndef YY_YYSLN_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SLNPARSE_SLN_TAB_HPP_INCLUDED
#define YY_YYSLN_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SLNPARSE_SLN_TAB_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#if YYDEBUG
extern int yysln_debug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
#define YYTOKENTYPE
enum yytokentype {
  TEXT_BLOCK = 258,
  CHAR_TOKEN = 259,
  DIGIT_TOKEN = 260,
  H_TOKEN = 261,
  H_BRACKET_TOKEN = 262,
  H_ASTERIX_TOKEN = 263,
  AT_TOKEN = 264,
  ATOM_TOKEN = 265,
  COMPARE_TOKEN = 266,
  OPEN_PAREN_TOKEN = 267,
  CLOSE_PAREN_TOKEN = 268,
  OPEN_BRACKET_TOKEN = 269,
  CLOSE_BRACKET_TOKEN = 270,
  OPEN_ANGLE_TOKEN = 271,
  CLOSE_ANGLE_TOKEN = 272,
  SEPARATOR_TOKEN = 273,
  ASTERIX_TOKEN = 274,
  EOS_TOKEN = 275,
  PLUS_TOKEN = 276,
  MINUS_TOKEN = 277,
  COLON_TOKEN = 278,
  EQUALS_TOKEN = 279,
  TILDE_TOKEN = 280,
  HASH_TOKEN = 281,
  COMMA_TOKEN = 282,
  NOT_TOKEN = 283,
  AND_TOKEN = 284,
  OR_TOKEN = 285,
  SEMI_TOKEN = 286,
  CARET_EQUALS_TOKEN = 287,
  COLON_EQUALS_TOKEN = 288,
  RECURSE_TOKEN = 289,
  NEG_RECURSE_TOKEN = 290,
  ERROR_TOKEN = 291
};
#endif

/* Value type.  */
#if !defined YYSTYPE && !defined YYSTYPE_IS_DECLARED

union YYSTYPE {
#line 92 "sln.yy" /* yacc.c:355  */

  int mol_T;
  RDKit::Atom *atom_T;
  RDKit::Bond *bond_T;
  int ival_T;
  std::string *text_T;
  char char_T;
  RDKit::SLNParse::AttribType *attrib_T;
  RDKit::SLNParse::AttribListType *attriblist_T;

#line 244 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
#define YYSTYPE_IS_TRIVIAL 1
#define YYSTYPE_IS_DECLARED 1
#endif

int yysln_parse(const char *input, std::vector<RDKit::RWMol *> *molList,
                bool doQueries, void *scanner);

#endif /* !YY_YYSLN_HOME_RODRIGUE_DOCUMENTS_CODE_RDKIT_BUILDER_RDKIT_CODE_GRAPHMOL_SLNPARSE_SLN_TAB_HPP_INCLUDED \
        */

/* Copy the second part of user declarations.  */

#line 260 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:358  */

#ifdef short
#undef short
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
#ifdef __SIZE_TYPE__
#define YYSIZE_T __SIZE_TYPE__
#elif defined size_t
#define YYSIZE_T size_t
#elif !defined YYSIZE_T
#include <cstddef> /* INFRINGES ON USER NAME SPACE */
#define YYSIZE_T size_t
#else
#define YYSIZE_T unsigned
#endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T)-1)

#ifndef YY_
#if defined YYENABLE_NLS && YYENABLE_NLS
#if ENABLE_NLS
#include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#define YY_(Msgid) dgettext("bison-runtime", Msgid)
#endif
#endif
#ifndef YY_
#define YY_(Msgid) Msgid
#endif
#endif

#ifndef YY_ATTRIBUTE
#if (defined __GNUC__ &&                                           \
     (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__))) || \
    defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#define YY_ATTRIBUTE(Spec) __attribute__(Spec)
#else
#define YY_ATTRIBUTE(Spec) /* empty */
#endif
#endif

#ifndef YY_ATTRIBUTE_PURE
#define YY_ATTRIBUTE_PURE YY_ATTRIBUTE((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
#define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE((__unused__))
#endif

#if !defined _Noreturn && \
    (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
#if defined _MSC_VER && 1200 <= _MSC_VER
#define _Noreturn __declspec(noreturn)
#else
#define _Noreturn YY_ATTRIBUTE((__noreturn__))
#endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if !defined lint || defined __GNUC__
#define YYUSE(E) ((void)(E))
#else
#define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && !defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                 \
  _Pragma("GCC diagnostic push")                            \
      _Pragma("GCC diagnostic ignored \"-Wuninitialized\"") \
          _Pragma("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
#define YY_IGNORE_MAYBE_UNINITIALIZED_END _Pragma("GCC diagnostic pop")
#else
#define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
#define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
#define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
#define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if !defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

#ifdef YYSTACK_USE_ALLOCA
#if YYSTACK_USE_ALLOCA
#ifdef __GNUC__
#define YYSTACK_ALLOC __builtin_alloca
#elif defined __BUILTIN_VA_ARG_INCR
#include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#elif defined _AIX
#define YYSTACK_ALLOC __alloca
#elif defined _MSC_VER
#include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#define alloca _alloca
#else
#define YYSTACK_ALLOC alloca
#if !defined _ALLOCA_H && !defined EXIT_SUCCESS
#include <cstdlib> /* INFRINGES ON USER NAME SPACE */
/* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#endif
#endif
#endif
#endif

#ifdef YYSTACK_ALLOC
/* Pacify GCC's 'empty if-body' warning.  */
#define YYSTACK_FREE(Ptr) \
  do { /* empty */        \
    ;                     \
  } while (0)
#ifndef YYSTACK_ALLOC_MAXIMUM
/* The OS might guarantee only one guard page at the bottom of the stack,
   and a page size can be as small as 4096 bytes.  So we cannot safely
   invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
   to allow for a few compiler-allocated temporary stack slots.  */
#define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#endif
#else
#define YYSTACK_ALLOC YYMALLOC
#define YYSTACK_FREE YYFREE
#ifndef YYSTACK_ALLOC_MAXIMUM
#define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#endif
#if (defined __cplusplus && !defined EXIT_SUCCESS && \
     !((defined YYMALLOC || defined malloc) &&       \
       (defined YYFREE || defined free)))
#include <cstdlib> /* INFRINGES ON USER NAME SPACE */
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#endif
#ifndef YYMALLOC
#define YYMALLOC malloc
#if !defined malloc && !defined EXIT_SUCCESS
void *malloc(YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#endif
#endif
#ifndef YYFREE
#define YYFREE free
#if !defined free && !defined EXIT_SUCCESS
void free(void *);      /* INFRINGES ON USER NAME SPACE */
#endif
#endif
#endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */

#if (!defined yyoverflow &&   \
     (!defined __cplusplus || \
      (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc {
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
#define YYSTACK_GAP_MAXIMUM (sizeof(union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
#define YYSTACK_BYTES(N) \
  ((N) * (sizeof(yytype_int16) + sizeof(YYSTYPE)) + YYSTACK_GAP_MAXIMUM)

#define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
#define YYSTACK_RELOCATE(Stack_alloc, Stack)                         \
  do {                                                               \
    YYSIZE_T yynewbytes;                                             \
    YYCOPY(&yyptr->Stack_alloc, Stack, yysize);                      \
    Stack = &yyptr->Stack_alloc;                                     \
    yynewbytes = yystacksize * sizeof(*Stack) + YYSTACK_GAP_MAXIMUM; \
    yyptr += yynewbytes / sizeof(*yyptr);                            \
  } while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
#ifndef YYCOPY
#if defined __GNUC__ && 1 < __GNUC__
#define YYCOPY(Dst, Src, Count) \
  __builtin_memcpy(Dst, Src, (Count) * sizeof(*(Src)))
#else
#define YYCOPY(Dst, Src, Count)                                  \
  do {                                                           \
    YYSIZE_T yyi;                                                \
    for (yyi = 0; yyi < (Count); yyi++) (Dst)[yyi] = (Src)[yyi]; \
  } while (0)
#endif
#endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL 25
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST 219

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS 37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS 16
/* YYNRULES -- Number of rules.  */
#define YYNRULES 72
/* YYNSTATES -- Number of states.  */
#define YYNSTATES 116

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK 2
#define YYMAXUTOK 291

#define YYTRANSLATE(YYX) \
  ((unsigned)(YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] = {
    0,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
    30, 31, 32, 33, 34, 35, 36};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] = {
    0,   144, 144, 145, 153, 158, 161, 168, 169, 174, 184, 187, 191, 195, 199,
    204, 208, 214, 218, 222, 233, 245, 258, 270, 284, 296, 310, 311, 312, 316,
    322, 330, 338, 348, 359, 360, 363, 367, 373, 381, 382, 389, 394, 403, 404,
    417, 426, 435, 445, 457, 461, 465, 469, 476, 481, 487, 495, 499, 510, 516,
    522, 528, 534, 540, 545, 559, 575, 589, 597, 607, 617, 629, 630};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] = {"$end",
                                      "error",
                                      "$undefined",
                                      "TEXT_BLOCK",
                                      "CHAR_TOKEN",
                                      "DIGIT_TOKEN",
                                      "H_TOKEN",
                                      "H_BRACKET_TOKEN",
                                      "H_ASTERIX_TOKEN",
                                      "AT_TOKEN",
                                      "ATOM_TOKEN",
                                      "COMPARE_TOKEN",
                                      "OPEN_PAREN_TOKEN",
                                      "CLOSE_PAREN_TOKEN",
                                      "OPEN_BRACKET_TOKEN",
                                      "CLOSE_BRACKET_TOKEN",
                                      "OPEN_ANGLE_TOKEN",
                                      "CLOSE_ANGLE_TOKEN",
                                      "SEPARATOR_TOKEN",
                                      "ASTERIX_TOKEN",
                                      "EOS_TOKEN",
                                      "PLUS_TOKEN",
                                      "MINUS_TOKEN",
                                      "COLON_TOKEN",
                                      "EQUALS_TOKEN",
                                      "TILDE_TOKEN",
                                      "HASH_TOKEN",
                                      "COMMA_TOKEN",
                                      "NOT_TOKEN",
                                      "AND_TOKEN",
                                      "OR_TOKEN",
                                      "SEMI_TOKEN",
                                      "CARET_EQUALS_TOKEN",
                                      "COLON_EQUALS_TOKEN",
                                      "RECURSE_TOKEN",
                                      "NEG_RECURSE_TOKEN",
                                      "ERROR_TOKEN",
                                      "$accept",
                                      "cmpd",
                                      "mol",
                                      "primmol",
                                      "atom",
                                      "hatom",
                                      "primatom",
                                      "bond",
                                      "primbond",
                                      "onebond",
                                      "attriblist",
                                      "ctabattriblist",
                                      "attrib",
                                      "recursivequery",
                                      "ctabattrib",
                                      "number",
                                      YY_NULLPTR};
#endif

#ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] = {
    0,   256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267,
    268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280,
    281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291};
#endif

#define YYPACT_NINF -40

#define yypact_value_is_default(Yystate) (!!((Yystate) == (-40)))

#define YYTABLE_NINF -67

#define yytable_value_is_error(Yytable_value) 0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] = {
    194, -13, -40, 28,  -40, -40, 5,   -6,  154, -40, -40, 46,  -40, 3,   -40,
    -40, 12,  15,  45,  194, 194, 7,   -40, 17,  4,   -40, 26,  48,  -40, 209,
    50,  107, -40, -40, -40, 39,  -40, -40, 199, 170, -40, 54,  28,  -40, 65,
    -40, -40, -40, 69,  75,  -40, 45,  45,  45,  194, -40, -40, 45,  -40, 141,
    -2,  -40, 154, 18,  50,  127, 204, 45,  -40, 50,  -40, 45,  -40, -40, 57,
    11,  -40, -40, -40, -40, 92,  95,  68,  71,  72,  -40, 48,  -40, 64,  86,
    50,  134, 139, 76,  152, -40, -40, 45,  -40, -40, -40, -40, -40, 88,  -40,
    96,  91,  -40, -40, -40, 156, -40, 97,  -40, -40, -40};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] = {
    0,  0,  9,  0,  30, 34, 0,  2,  7,  10, 26, 27, 6,  55, 71, 62, 58,
    60, 0,  0,  0,  0,  52, 63, 0,  1,  0,  0,  5,  0,  0,  0,  45, 48,
    46, 41, 47, 11, 0,  39, 43, 28, 0,  35, 0,  59, 61, 56, 0,  0,  33,
    0,  0,  0,  0,  72, 31, 0,  4,  67, 0,  53, 8,  13, 0,  0,  0,  0,
    19, 0,  12, 0,  44, 29, 0,  0,  57, 49, 50, 51, 0,  0,  0,  0,  0,
    3,  0,  20, 0,  15, 0,  0,  0,  14, 0,  38, 36, 0,  32, 68, 70, 69,
    54, 17, 22, 0,  16, 42, 21, 40, 0,  24, 18, 23, 37, 25};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] = {-40, -9, -40, -27, -8,  -40, -40, 47,
                                      -40, 79, -39, -40, -10, -40, 25,  -29};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] = {-1, 6,  7,  8,  9,  10, 11, 38,
                                        39, 40, 21, 60, 22, 23, 61, 24};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] = {
    37,  63,  62,  74,  65,  25, 26,  12,  47,  55,  48,  49, 29,  75,  44,
    85,  55,  45,  81,  56,  46, 27,  50,  55,  87,  28,  96, 57,  92,  86,
    70,  13,  94,  14,  97,  88, 51,  52,  53,  91,  93,  77, 78,  79,  54,
    80,  58,  15,  13,  16,  17, 59,  41,  67,  37,  14,  18, 37,  110, 73,
    42,  105, 19,  20,  15,  43, 16,  17,  76,  55,  26,  99, 95,  18,  100,
    101, 26,  103, 66,  19,  20, 55,  108, 37,  -64, 27,  51, 52,  53,  28,
    -65, 27,  104, 26,  111, 28, -64, 113, -64, -64, -64, 55, -65, 115, -65,
    -65, -65, -66, 27,  112, 98, 102, 28,  2,   3,   4,   64, 5,   72,  -66,
    0,   -66, -66, -66, 51,  52, 53,  0,   0,   32,  33,  34, 35,  36,  3,
    4,   30,  5,   0,   31,  89, 3,   4,   30,  5,   0,   31, 106, 0,   32,
    33,  34,  35,  36,  107, 0,  32,  33,  34,  35,  36,  3,  4,   30,  5,
    82,  31,  109, 51,  52,  53, 114, 0,   83,  84,  0,   32, 33,  34,  35,
    36,  51,  52,  53,  71,  51, 52,  53,  0,   0,   0,   0,  32,  33,  34,
    1,   36,  0,   0,   0,   2,  3,   4,   0,   5,   68,  3,  4,   69,  5,
    2,   3,   4,   90,  5,   2,  3,   4,   0,   5};

static const yytype_int8 yycheck[] = {
    8,  30, 29, 42, 31, 0,  1,  20, 18, 5,  19, 20, 18, 42, 11, 17, 5,  5,  57,
    15, 5,  16, 15, 5,  6,  20, 15, 23, 67, 31, 38, 3,  71, 5,  23, 64, 29, 30,
    31, 66, 69, 51, 52, 53, 27, 54, 20, 19, 3,  21, 22, 3,  6,  14, 62, 5,  28,
    65, 97, 5,  14, 90, 34, 35, 19, 19, 21, 22, 3,  5,  1,  3,  15, 28, 3,  3,
    1,  13, 31, 34, 35, 5,  6,  91, 15, 16, 29, 30, 31, 20, 15, 16, 6,  1,  6,
    20, 27, 6,  29, 30, 31, 5,  27, 6,  29, 30, 31, 15, 16, 13, 15, 86, 20, 6,
    7,  8,  9,  10, 39, 27, -1, 29, 30, 31, 29, 30, 31, -1, -1, 22, 23, 24, 25,
    26, 7,  8,  9,  10, -1, 12, 13, 7,  8,  9,  10, -1, 12, 13, -1, 22, 23, 24,
    25, 26, 15, -1, 22, 23, 24, 25, 26, 7,  8,  9,  10, 24, 12, 15, 29, 30, 31,
    15, -1, 32, 33, -1, 22, 23, 24, 25, 26, 29, 30, 31, 14, 29, 30, 31, -1, -1,
    -1, -1, 22, 23, 24, 1,  26, -1, -1, -1, 6,  7,  8,  -1, 10, 6,  7,  8,  9,
    10, 6,  7,  8,  9,  10, 6,  7,  8,  -1, 10};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] = {
    0,  1,  6,  7,  8,  10, 38, 39, 40, 41, 42, 43, 20, 3,  5,  19, 21,
    22, 28, 34, 35, 47, 49, 50, 52, 0,  1,  16, 20, 18, 9,  12, 22, 23,
    24, 25, 26, 41, 44, 45, 46, 6,  14, 19, 11, 5,  5,  49, 38, 38, 15,
    29, 30, 31, 27, 5,  15, 23, 20, 3,  48, 51, 40, 52, 9,  40, 44, 14,
    6,  9,  41, 14, 46, 5,  47, 52, 3,  49, 49, 49, 38, 47, 24, 32, 33,
    17, 31, 6,  52, 13, 9,  40, 47, 52, 47, 15, 15, 23, 15, 3,  3,  3,
    51, 13, 6,  52, 13, 15, 6,  15, 47, 6,  13, 6,  15, 6};

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] = {
    0,  37, 38, 38, 38, 38, 38, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 42, 42, 42, 42, 43, 43, 43, 43,
    43, 44, 44, 44, 44, 45, 45, 46, 46, 46, 46, 47, 47, 47, 47, 48, 48, 49, 49,
    49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 51, 51, 51, 51, 52, 52};

/* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] = {
    0, 2, 1, 4, 3, 2, 2, 1, 3, 1, 1, 2, 3, 3, 4, 4, 5, 5, 6, 3, 4, 5, 5, 6, 6,
    7, 1, 1, 2, 3, 1, 3, 5, 3, 1, 2, 4, 6, 4, 1, 4, 1, 4, 1, 2, 1, 1, 1, 1, 3,
    3, 3, 1, 1, 3, 1, 2, 3, 1, 2, 1, 2, 1, 1, 2, 2, 3, 1, 3, 3, 3, 1, 2};

#define yyerrok (yyerrstatus = 0)
#define yyclearin (yychar = YYEMPTY)
#define YYEMPTY (-2)
#define YYEOF 0

#define YYACCEPT goto yyacceptlab
#define YYABORT goto yyabortlab
#define YYERROR goto yyerrorlab

#define YYRECOVERING() (!!yyerrstatus)

#define YYBACKUP(Token, Value)                      \
  do                                                \
    if (yychar == YYEMPTY) {                        \
      yychar = (Token);                             \
      yylval = (Value);                             \
      YYPOPSTACK(yylen);                            \
      yystate = *yyssp;                             \
      goto yybackup;                                \
    } else {                                        \
      yyerror(input, molList, doQueries, scanner,   \
              YY_("syntax error: cannot back up")); \
      YYERROR;                                      \
    }                                               \
  while (0)

/* Error token number */
#define YYTERROR 1
#define YYERRCODE 256

/* Enable debugging if requested.  */
#if YYDEBUG

#ifndef YYFPRINTF
#include <cstdio> /* INFRINGES ON USER NAME SPACE */
#define YYFPRINTF fprintf
#endif

#define YYDPRINTF(Args)          \
  do {                           \
    if (yydebug) YYFPRINTF Args; \
  } while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
#define YY_LOCATION_PRINT(File, Loc) ((void)0)
#endif

#define YY_SYMBOL_PRINT(Title, Type, Value, Location)                 \
  do {                                                                \
    if (yydebug) {                                                    \
      YYFPRINTF(stderr, "%s ", Title);                                \
      yy_symbol_print(stderr, Type, Value, input, molList, doQueries, \
                      scanner);                                       \
      YYFPRINTF(stderr, "\n");                                        \
    }                                                                 \
  } while (0)

/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void yy_symbol_value_print(FILE *yyoutput, int yytype,
                                  YYSTYPE const *const yyvaluep,
                                  const char *input,
                                  std::vector<RDKit::RWMol *> *molList,
                                  bool doQueries, void *scanner) {
  FILE *yyo = yyoutput;
  YYUSE(yyo);
  YYUSE(input);
  YYUSE(molList);
  YYUSE(doQueries);
  YYUSE(scanner);
  if (!yyvaluep) return;
#ifdef YYPRINT
  if (yytype < YYNTOKENS) YYPRINT(yyoutput, yytoknum[yytype], *yyvaluep);
#endif
  YYUSE(yytype);
}

/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void yy_symbol_print(FILE *yyoutput, int yytype,
                            YYSTYPE const *const yyvaluep, const char *input,
                            std::vector<RDKit::RWMol *> *molList,
                            bool doQueries, void *scanner) {
  YYFPRINTF(yyoutput, "%s %s (", yytype < YYNTOKENS ? "token" : "nterm",
            yytname[yytype]);

  yy_symbol_value_print(yyoutput, yytype, yyvaluep, input, molList, doQueries,
                        scanner);
  YYFPRINTF(yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void yy_stack_print(yytype_int16 *yybottom, yytype_int16 *yytop) {
  YYFPRINTF(stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++) {
    int yybot = *yybottom;
    YYFPRINTF(stderr, " %d", yybot);
  }
  YYFPRINTF(stderr, "\n");
}

#define YY_STACK_PRINT(Bottom, Top)               \
  do {                                            \
    if (yydebug) yy_stack_print((Bottom), (Top)); \
  } while (0)

/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void yy_reduce_print(yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule,
                            const char *input,
                            std::vector<RDKit::RWMol *> *molList,
                            bool doQueries, void *scanner) {
  unsigned long yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF(stderr, "Reducing stack by rule %d (line %lu):\n", yyrule - 1,
            yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++) {
    YYFPRINTF(stderr, "   $%d = ", yyi + 1);
    yy_symbol_print(stderr, yystos[yyssp[yyi + 1 - yynrhs]],
                    &(yyvsp[(yyi + 1) - (yynrhs)]), input, molList, doQueries,
                    scanner);
    YYFPRINTF(stderr, "\n");
  }
}

#define YY_REDUCE_PRINT(Rule)                                                  \
  do {                                                                         \
    if (yydebug)                                                               \
      yy_reduce_print(yyssp, yyvsp, Rule, input, molList, doQueries, scanner); \
  } while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
#define YYDPRINTF(Args)
#define YY_SYMBOL_PRINT(Title, Type, Value, Location)
#define YY_STACK_PRINT(Bottom, Top)
#define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

#if YYERROR_VERBOSE

#ifndef yystrlen
#if defined __GLIBC__ && defined _STRING_H
#define yystrlen strlen
#else
/* Return the length of YYSTR.  */
static YYSIZE_T yystrlen(const char *yystr) {
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++) continue;
  return yylen;
}
#endif
#endif

#ifndef yystpcpy
#if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#define yystpcpy stpcpy
#else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *yystpcpy(char *yydest, const char *yysrc) {
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0') continue;

  return yyd - 1;
}
#endif
#endif

#ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T yytnamerr(char *yyres, const char *yystr) {
  if (*yystr == '"') {
    YYSIZE_T yyn = 0;
    char const *yyp = yystr;

    for (;;) switch (*++yyp) {
        case '\'':
        case ',':
          goto do_not_strip_quotes;

        case '\\':
          if (*++yyp != '\\') goto do_not_strip_quotes;
          /* Fall through.  */
        default:
          if (yyres) yyres[yyn] = *yyp;
          yyn++;
          break;

        case '"':
          if (yyres) yyres[yyn] = '\0';
          return yyn;
      }
  do_not_strip_quotes:;
  }

  if (!yyres) return yystrlen(yystr);

  return yystpcpy(yyres, yystr) - yyres;
}
#endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int yysyntax_error(YYSIZE_T *yymsg_alloc, char **yymsg,
                          yytype_int16 *yyssp, int yytoken) {
  YYSIZE_T yysize0 = yytnamerr(YY_NULLPTR, yytname[yytoken]);
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
  if (yytoken != YYEMPTY) {
    int yyn = yypact[*yyssp];
    yyarg[yycount++] = yytname[yytoken];
    if (!yypact_value_is_default(yyn)) {
      /* Start YYX at -YYN if negative to avoid negative indexes in
         YYCHECK.  In other words, skip the first -YYN actions for
         this state because they are default actions.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;
      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yyx;

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
        if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR &&
            !yytable_value_is_error(yytable[yyx + yyn])) {
          if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM) {
            yycount = 1;
            yysize = yysize0;
            break;
          }
          yyarg[yycount++] = yytname[yyx];
          {
            YYSIZE_T yysize1 = yysize + yytnamerr(YY_NULLPTR, yytname[yyx]);
            if (!(yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
              return 2;
            yysize = yysize1;
          }
        }
    }
  }

  switch (yycount) {
#define YYCASE_(N, S) \
  case N:             \
    yyformat = S;     \
    break
    default: /* Avoid compiler warnings. */
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(
          5,
          YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
  }

  {
    YYSIZE_T yysize1 = yysize + yystrlen(yyformat);
    if (!(yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM)) return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize) {
    *yymsg_alloc = 2 * yysize;
    if (!(yysize <= *yymsg_alloc && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
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
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount) {
        yyp += yytnamerr(yyp, yyarg[yyi++]);
        yyformat += 2;
      } else {
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

static void yydestruct(const char *yymsg, int yytype, YYSTYPE *yyvaluep,
                       const char *input, std::vector<RDKit::RWMol *> *molList,
                       bool doQueries, void *scanner) {
  YYUSE(yyvaluep);
  YYUSE(input);
  YYUSE(molList);
  YYUSE(doQueries);
  YYUSE(scanner);
  if (!yymsg) yymsg = "Deleting";
  YY_SYMBOL_PRINT(yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  switch (yytype) {
    case 3:        /* TEXT_BLOCK  */
#line 136 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).text_T);
    }
#line 1196 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 10:       /* ATOM_TOKEN  */
#line 134 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).atom_T);
    }
#line 1202 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 11:       /* COMPARE_TOKEN  */
#line 136 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).text_T);
    }
#line 1208 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 41:       /* atom  */
#line 134 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).atom_T);
    }
#line 1214 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 42:       /* hatom  */
#line 134 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).atom_T);
    }
#line 1220 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 43:       /* primatom  */
#line 134 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).atom_T);
    }
#line 1226 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 44:       /* bond  */
#line 135 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).bond_T);
    }
#line 1232 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 45:       /* primbond  */
#line 135 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).bond_T);
    }
#line 1238 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 46:       /* onebond  */
#line 135 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).bond_T);
    }
#line 1244 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 47:       /* attriblist  */
#line 138 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).attriblist_T);
    }
#line 1250 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 48:       /* ctabattriblist  */
#line 138 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).attriblist_T);
    }
#line 1256 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 49:       /* attrib  */
#line 137 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).attrib_T);
    }
#line 1262 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 50:       /* recursivequery  */
#line 137 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).attrib_T);
    }
#line 1268 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    case 51:       /* ctabattrib  */
#line 137 "sln.yy" /* yacc.c:1258  */
    {
      delete ((*yyvaluep).attrib_T);
    }
#line 1274 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1258  */
    break;

    default:
      break;
  }
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}

/*----------.
| yyparse.  |
`----------*/

int yyparse(const char *input, std::vector<RDKit::RWMol *> *molList,
            bool doQueries, void *scanner) {
  /* The lookahead symbol.  */
  int yychar;

  /* The semantic value of the lookahead symbol.  */
  /* Default value used for initialization, for pacifying older GCCs
     or non-GCC compilers.  */
  YY_INITIAL_VALUE(static YYSTYPE yyval_default;)
  YYSTYPE yylval YY_INITIAL_VALUE(= yyval_default);

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

#define YYPOPSTACK(N) (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF((stderr, "Starting parse\n"));

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

  if (yyss + yystacksize - 1 <= yyssp) {
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
      yyoverflow(YY_("memory exhausted"), &yyss1, yysize * sizeof(*yyssp),
                 &yyvs1, yysize * sizeof(*yyvsp), &yystacksize);

      yyss = yyss1;
      yyvs = yyvs1;
    }
#else /* no yyoverflow */
#ifndef YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    /* Extend the stack our own way.  */
    if (YYMAXDEPTH <= yystacksize) goto yyexhaustedlab;
    yystacksize *= 2;
    if (YYMAXDEPTH < yystacksize) yystacksize = YYMAXDEPTH;

    {
      yytype_int16 *yyss1 = yyss;
      union yyalloc *yyptr =
          (union yyalloc *)YYSTACK_ALLOC(YYSTACK_BYTES(yystacksize));
      if (!yyptr) goto yyexhaustedlab;
      YYSTACK_RELOCATE(yyss_alloc, yyss);
      YYSTACK_RELOCATE(yyvs_alloc, yyvs);
#undef YYSTACK_RELOCATE
      if (yyss1 != yyssa) YYSTACK_FREE(yyss1);
    }
#endif
#endif /* no yyoverflow */

    yyssp = yyss + yysize - 1;
    yyvsp = yyvs + yysize - 1;

    YYDPRINTF(
        (stderr, "Stack size increased to %lu\n", (unsigned long)yystacksize));

    if (yyss + yystacksize - 1 <= yyssp) YYABORT;
  }

  YYDPRINTF((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL) YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default(yyn)) goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY) {
    YYDPRINTF((stderr, "Reading a token: "));
    yychar = yylex(&yylval, scanner);
  }

  if (yychar <= YYEOF) {
    yychar = yytoken = YYEOF;
    YYDPRINTF((stderr, "Now at end of input.\n"));
  } else {
    yytoken = YYTRANSLATE(yychar);
    YY_SYMBOL_PRINT("Next token is", yytoken, &yylval, &yylloc);
  }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken) goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0) {
    if (yytable_value_is_error(yyn)) goto yyerrlab;
    yyn = -yyn;
    goto yyreduce;
  }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus) yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT("Shifting", yytoken, &yylval, &yylloc);

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
  if (yyn == 0) goto yyerrlab;
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
  yyval = yyvsp[1 - yylen];

  YY_REDUCE_PRINT(yyn);
  switch (yyn) {
    case 3:
#line 145 "sln.yy" /* yacc.c:1651  */
    {
      // allowing mol<attrs><attrs> seems to be a NIBR thing, I don't
      // think it's standard SLN
      RDKit::ROMol *mol = (*molList)[(yyvsp[-3].mol_T)];
      SLNParse::parseMolAttribs(mol, *(yyvsp[-1].attriblist_T));
      delete (yyvsp[-1].attriblist_T);
      (yyval.mol_T) = (yyvsp[-3].mol_T);
    }
#line 1549 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 4:
#line 153 "sln.yy" /* yacc.c:1651  */
    {
      yyclearin;
      yyerrok;
      YYABORT;
    }
#line 1559 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 5:
#line 158 "sln.yy" /* yacc.c:1651  */
    {
      YYACCEPT;
    }
#line 1567 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 6:
#line 161 "sln.yy" /* yacc.c:1651  */
    {
      yyclearin;
      yyerrok;
      YYABORT;
    }
#line 1577 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 8:
#line 169 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.mol_T) =
          SLNParse::addFragToMol(*molList, (yyvsp[-2].mol_T), (yyvsp[0].mol_T));
    }
#line 1585 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 9:
#line 174 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }

      (yyval.mol_T) = SLNParse::startMol(*molList, newAtom, doQueries);
    }
#line 1600 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 10:
#line 184 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.mol_T) =
          SLNParse::startMol(*molList, (yyvsp[0].atom_T), doQueries);
    }
#line 1608 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 11:
#line 187 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), (yyvsp[0].atom_T),
                             doQueries);
      (yyval.mol_T) = (yyvsp[-1].mol_T);
    }
#line 1617 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 12:
#line 191 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), (yyvsp[0].atom_T),
                             (yyvsp[-1].bond_T), doQueries);
      (yyval.mol_T) = (yyvsp[-2].mol_T);
    }
#line 1626 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 13:
#line 195 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[0].ival_T));
      (yyval.mol_T) = (yyvsp[-2].mol_T);
    }
#line 1635 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 14:
#line 199 "sln.yy" /* yacc.c:1651  */
    {
      // closeRingBond() takes ownership of the bond
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[0].ival_T),
                              (yyvsp[-2].bond_T));
      (yyval.mol_T) = (yyvsp[-3].mol_T);
    }
#line 1645 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 15:
#line 204 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::addBranchToMol(*molList, (yyval.mol_T), (yyvsp[-1].mol_T));
      (yyval.mol_T) = (yyvsp[-3].mol_T);
    }
#line 1654 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 16:
#line 208 "sln.yy" /* yacc.c:1651  */
    {
      // addBranchToMol() takes ownership of the bond and deletes the
      // branch, so no leaks here'
      SLNParse::addBranchToMol(*molList, (yyval.mol_T), (yyvsp[-1].mol_T),
                               (yyvsp[-2].bond_T));
      (yyval.mol_T) = (yyvsp[-4].mol_T);
    }
#line 1665 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 17:
#line 214 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[-1].ival_T));
      (yyval.mol_T) = (yyvsp[-4].mol_T);
    }
#line 1674 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 18:
#line 218 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[-1].ival_T),
                              (yyvsp[-3].bond_T));
      (yyval.mol_T) = (yyvsp[-5].mol_T);
    }
#line 1683 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 19:
#line 222 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }

      SLNParse::addAtomToMol(*molList, (yyval.mol_T), newAtom,
                             (yyvsp[-1].bond_T), doQueries);
      (yyval.mol_T) = (yyvsp[-2].mol_T);
    }
#line 1699 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 20:
#line 233 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[-1].ival_T));
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), newAtom, doQueries);

      (yyval.mol_T) = (yyvsp[-3].mol_T);
    }
#line 1716 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 21:
#line 245 "sln.yy" /* yacc.c:1651  */
    {
      // closeRingBond() takes ownership of the bond
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[-1].ival_T),
                              (yyvsp[-3].bond_T));
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), newAtom, doQueries);

      (yyval.mol_T) = (yyvsp[-4].mol_T);
    }
#line 1734 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 22:
#line 258 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::addBranchToMol(*molList, (yyval.mol_T), (yyvsp[-2].mol_T));
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), newAtom, doQueries);

      (yyval.mol_T) = (yyvsp[-4].mol_T);
    }
#line 1751 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 23:
#line 270 "sln.yy" /* yacc.c:1651  */
    {
      // addBranchToMol() takes ownership of the bond and deletes the
      // branch, so no leaks here'
      SLNParse::addBranchToMol(*molList, (yyval.mol_T), (yyvsp[-2].mol_T),
                               (yyvsp[-3].bond_T));
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), newAtom, doQueries);

      (yyval.mol_T) = (yyvsp[-5].mol_T);
    }
#line 1770 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 24:
#line 284 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[-2].ival_T));
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), newAtom, doQueries);

      (yyval.mol_T) = (yyvsp[-5].mol_T);
    }
#line 1787 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 25:
#line 296 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::closeRingBond(*molList, (yyval.mol_T), (yyvsp[-2].ival_T),
                              (yyvsp[-4].bond_T));
      RDKit::Atom *newAtom;
      if (!doQueries) {
        newAtom = new RDKit::Atom(1);
      } else {
        newAtom = new RDKit::QueryAtom(1);
      }
      SLNParse::addAtomToMol(*molList, (yyval.mol_T), newAtom, doQueries);

      (yyval.mol_T) = (yyvsp[-6].mol_T);
    }
#line 1804 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 28:
#line 312 "sln.yy" /* yacc.c:1651  */
    {
      (yyvsp[-1].atom_T)->setNumExplicitHs(1);
      (yyval.atom_T) = (yyvsp[-1].atom_T);
    }
#line 1813 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 29:
#line 316 "sln.yy" /* yacc.c:1651  */
    {
      (yyvsp[-2].atom_T)->setNumExplicitHs((yyvsp[0].ival_T));
      (yyval.atom_T) = (yyvsp[-2].atom_T);
    }
#line 1822 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 30:
#line 322 "sln.yy" /* yacc.c:1651  */
    {
      if (!doQueries) {
        (yyval.atom_T) = new RDKit::Atom(1);
      } else {
        (yyval.atom_T) = new RDKit::QueryAtom(1);
      }
      (yyval.atom_T)->setProp(RDKit::common_properties::_starred, 1, true);
    }
#line 1835 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 31:
#line 330 "sln.yy" /* yacc.c:1651  */
    {
      if (!doQueries) {
        (yyval.atom_T) = new RDKit::Atom(1);
      } else {
        (yyval.atom_T) = new RDKit::QueryAtom(1);
      }
      (yyval.atom_T)
          ->setProp(RDKit::common_properties::_AtomID,
                    static_cast<unsigned int>((yyvsp[-1].ival_T)));
    }
#line 1848 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 32:
#line 338 "sln.yy" /* yacc.c:1651  */
    {
      if (!doQueries) {
        (yyval.atom_T) = new RDKit::Atom(1);
      } else {
        (yyval.atom_T) = new RDKit::QueryAtom(1);
      }
      (yyval.atom_T)
          ->setProp(RDKit::common_properties::_AtomID,
                    static_cast<unsigned int>((yyvsp[-3].ival_T)));
      SLNParse::parseAtomAttribs((yyval.atom_T), *(yyvsp[-1].attriblist_T),
                                 doQueries);
      delete (yyvsp[-1].attriblist_T);
    }
#line 1863 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 33:
#line 348 "sln.yy" /* yacc.c:1651  */
    {
      if (!doQueries) {
        (yyval.atom_T) = new RDKit::Atom(1);
      } else {
        (yyval.atom_T) = new RDKit::QueryAtom(1);
      }
      SLNParse::parseAtomAttribs((yyval.atom_T), *(yyvsp[-1].attriblist_T),
                                 doQueries);
      delete (yyvsp[-1].attriblist_T);
    }
#line 1877 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 35:
#line 360 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.atom_T)->setProp(RDKit::common_properties::_starred, 1, true);
    }
#line 1885 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 36:
#line 363 "sln.yy" /* yacc.c:1651  */
    {
      (yyvsp[-3].atom_T)
          ->setProp(RDKit::common_properties::_AtomID,
                    static_cast<unsigned int>((yyvsp[-1].ival_T)));
      (yyval.atom_T) = (yyvsp[-3].atom_T);
    }
#line 1894 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 37:
#line 367 "sln.yy" /* yacc.c:1651  */
    {
      (yyvsp[-5].atom_T)
          ->setProp(RDKit::common_properties::_AtomID,
                    static_cast<unsigned int>((yyvsp[-3].ival_T)));
      SLNParse::parseAtomAttribs((yyvsp[-5].atom_T), *(yyvsp[-1].attriblist_T),
                                 doQueries);
      delete (yyvsp[-1].attriblist_T);
      (yyval.atom_T) = (yyvsp[-5].atom_T);
    }
#line 1905 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 38:
#line 373 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::parseAtomAttribs((yyvsp[-3].atom_T), *(yyvsp[-1].attriblist_T),
                                 doQueries);
      delete (yyvsp[-1].attriblist_T);
      (yyval.atom_T) = (yyvsp[-3].atom_T);
    }
#line 1915 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 40:
#line 382 "sln.yy" /* yacc.c:1651  */
    {
      SLNParse::parseBondAttribs((yyvsp[-3].bond_T), *(yyvsp[-1].attriblist_T),
                                 doQueries);
      delete (yyvsp[-1].attriblist_T);
      (yyval.bond_T) = (yyvsp[-3].bond_T);
    }
#line 1925 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 41:
#line 389 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Bond *bond = new RDKit::QueryBond();
      bond->setQuery(RDKit::makeBondNullQuery());
      (yyval.bond_T) = bond;
    }
#line 1935 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 42:
#line 394 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Bond *bond = new RDKit::QueryBond();
      bond->setQuery(RDKit::makeBondNullQuery());
      SLNParse::parseBondAttribs(bond, *(yyvsp[-1].attriblist_T), doQueries);
      delete (yyvsp[-1].attriblist_T);
      (yyval.bond_T) = bond;
    }
#line 1947 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 44:
#line 404 "sln.yy" /* yacc.c:1651  */
    {
      if (!doQueries) {
        yysln_error(input, molList, doQueries, nullptr,
                    "sequential bonds not allowed in non-queries");
        YYABORT;
      } else {
        RDKit::QueryBond *b1 =
            static_cast<RDKit::QueryBond *>((yyvsp[-1].bond_T));
        RDKit::QueryBond *b2 =
            static_cast<RDKit::QueryBond *>((yyvsp[0].bond_T));
        b1->expandQuery(b2->getQuery()->copy(), Queries::COMPOSITE_OR, true);
        delete b2;
      }
    }
#line 1963 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 45:
#line 417 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Bond *bond;
      if (doQueries) {
        bond = new RDKit::QueryBond(RDKit::Bond::SINGLE);
      } else {
        bond = new RDKit::Bond(RDKit::Bond::SINGLE);
      }
      (yyval.bond_T) = bond;
    }
#line 1977 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 46:
#line 426 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Bond *bond;
      if (doQueries) {
        bond = new RDKit::QueryBond(RDKit::Bond::DOUBLE);
      } else {
        bond = new RDKit::Bond(RDKit::Bond::DOUBLE);
      }
      (yyval.bond_T) = bond;
    }
#line 1991 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 47:
#line 435 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Bond *bond;
      if (doQueries) {
        bond = new RDKit::QueryBond(RDKit::Bond::TRIPLE);
      } else {
        bond = new RDKit::Bond(RDKit::Bond::TRIPLE);
      }
      (yyval.bond_T) = bond;

    }
#line 2006 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 48:
#line 445 "sln.yy" /* yacc.c:1651  */
    {
      RDKit::Bond *bond;
      if (doQueries) {
        bond = new RDKit::QueryBond(RDKit::Bond::AROMATIC);
      } else {
        bond = new RDKit::Bond(RDKit::Bond::AROMATIC);
      }
      (yyval.bond_T) = bond;
    }
#line 2020 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 49:
#line 457 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attriblist_T)
          ->push_back(std::make_pair(
              SLNParse::AttribAnd,
              boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
    }
#line 2029 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 50:
#line 461 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attriblist_T)
          ->push_back(std::make_pair(
              SLNParse::AttribOr,
              boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
    }
#line 2038 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 51:
#line 465 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attriblist_T)
          ->push_back(std::make_pair(
              SLNParse::AttribLowPriAnd,
              boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
    }
#line 2047 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 52:
#line 469 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attriblist_T) = new SLNParse::AttribListType();
      (yyval.attriblist_T)
          ->push_back(std::make_pair(
              SLNParse::AttribLowPriAnd,
              boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
    }
#line 2057 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 53:
#line 476 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attriblist_T) = new SLNParse::AttribListType();
      (yyval.attriblist_T)
          ->push_back(std::make_pair(
              SLNParse::AttribAnd,
              boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
    }
#line 2067 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 54:
#line 481 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attriblist_T)
          ->push_back(std::make_pair(
              SLNParse::AttribAnd,
              boost::shared_ptr<SLNParse::AttribType>((yyvsp[0].attrib_T))));
    }
#line 2076 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 55:
#line 487 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = *(yyvsp[0].text_T);
      boost::to_lower((yyval.attrib_T)->first);
      (yyval.attrib_T)->op = "";
      (yyval.attrib_T)->second = "";
      delete (yyvsp[0].text_T);
    }
#line 2089 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 56:
#line 495 "sln.yy" /* yacc.c:1651  */
    {
      (yyvsp[0].attrib_T)->negated = true;
      (yyval.attrib_T) = (yyvsp[0].attrib_T);
    }
#line 2098 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 57:
#line 499 "sln.yy" /* yacc.c:1651  */
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
#line 2114 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 58:
#line 510 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = "charge";
      (yyval.attrib_T)->op = "=";
      (yyval.attrib_T)->second = "+1";
    }
#line 2125 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 59:
#line 516 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = "charge";
      (yyval.attrib_T)->op = "=";
      (yyval.attrib_T)->second = SLNParse::convertToString((yyvsp[0].ival_T));
    }
#line 2136 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 60:
#line 522 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = "charge";
      (yyval.attrib_T)->op = "=";
      (yyval.attrib_T)->second = "-1";
    }
#line 2147 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 61:
#line 528 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = "charge";
      (yyval.attrib_T)->op = "=";
      (yyval.attrib_T)->second = SLNParse::convertToString(-(yyvsp[0].ival_T));
    }
#line 2158 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 62:
#line 534 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = "spin";
      (yyval.attrib_T)->op = "=";
      (yyval.attrib_T)->second = "d";
    }
#line 2169 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 63:
#line 540 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = (yyvsp[0].attrib_T);
    }
#line 2177 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 64:
#line 545 "sln.yy" /* yacc.c:1651  */
    {
      int sz = molList->size();
      RDKit::ROMol *mol = (*molList)[(yyvsp[0].mol_T)];
      molList->resize(sz - 1);
      SLNParse::finalizeQueryMol(mol, true);
      RDKit::RecursiveStructureQuery *rsq =
          new RDKit::RecursiveStructureQuery(mol);
      RDKit::ATOM_OR_QUERY *orq = new RDKit::ATOM_OR_QUERY();
      orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = "is";
      (yyval.attrib_T)->op = "=";
      (yyval.attrib_T)->second = "";
      (yyval.attrib_T)->structQuery = static_cast<void *>(orq);
    }
#line 2196 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 65:
#line 559 "sln.yy" /* yacc.c:1651  */
    {
      int sz = molList->size();
      RDKit::ROMol *mol = (*molList)[(yyvsp[0].mol_T)];
      molList->resize(sz - 1);
      SLNParse::finalizeQueryMol(mol, true);
      RDKit::RecursiveStructureQuery *rsq =
          new RDKit::RecursiveStructureQuery(mol);
      RDKit::ATOM_OR_QUERY *orq = new RDKit::ATOM_OR_QUERY();
      orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
      orq->setNegation(true);

      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = "is";
      (yyval.attrib_T)->op = "=";
      (yyval.attrib_T)->second = "";
      (yyval.attrib_T)->structQuery = static_cast<void *>(orq);
    }
#line 2217 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 66:
#line 575 "sln.yy" /* yacc.c:1651  */
    {
      int sz = molList->size();
      RDKit::ROMol *mol = (*molList)[(yyvsp[0].mol_T)];
      molList->resize(sz - 1);
      SLNParse::finalizeQueryMol(mol, true);
      RDKit::RecursiveStructureQuery *rsq =
          new RDKit::RecursiveStructureQuery(mol);

      RDKit::ATOM_OR_QUERY *orq = static_cast<RDKit::ATOM_OR_QUERY *>(
          (yyvsp[-2].attrib_T)->structQuery);
      orq->addChild(RDKit::ATOM_OR_QUERY::CHILD_TYPE(rsq));
      (yyval.attrib_T) = (yyvsp[-2].attrib_T);
    }
#line 2233 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 67:
#line 589 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.attrib_T) = new SLNParse::AttribType();
      (yyval.attrib_T)->first = *(yyvsp[0].text_T);
      boost::to_lower((yyval.attrib_T)->first);
      (yyval.attrib_T)->op = "";
      (yyval.attrib_T)->second = "";
      delete (yyvsp[0].text_T);
    }
#line 2246 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 68:
#line 597 "sln.yy" /* yacc.c:1651  */
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
#line 2261 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 69:
#line 607 "sln.yy" /* yacc.c:1651  */
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
#line 2276 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 70:
#line 617 "sln.yy" /* yacc.c:1651  */
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
#line 2291 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

    case 72:
#line 630 "sln.yy" /* yacc.c:1651  */
    {
      (yyval.ival_T) = (yyvsp[-1].ival_T) * 10 + (yyvsp[0].ival_T);
    }
#line 2297 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    break;

#line 2301 "/home/rodrigue/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SLNParse/sln.tab.cpp" /* yacc.c:1651  */
    default:
      break;
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
  YY_SYMBOL_PRINT("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK(yylen);
  yylen = 0;
  YY_STACK_PRINT(yyss, yyssp);

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
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE(yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus) {
    ++yynerrs;
#if !YYERROR_VERBOSE
    yyerror(input, molList, doQueries, scanner, YY_("syntax error"));
#else
#define YYSYNTAX_ERROR yysyntax_error(&yymsg_alloc, &yymsg, yyssp, yytoken)
    {
      char const *yymsgp = YY_("syntax error");
      int yysyntax_error_status;
      yysyntax_error_status = YYSYNTAX_ERROR;
      if (yysyntax_error_status == 0)
        yymsgp = yymsg;
      else if (yysyntax_error_status == 1) {
        if (yymsg != yymsgbuf) YYSTACK_FREE(yymsg);
        yymsg = (char *)YYSTACK_ALLOC(yymsg_alloc);
        if (!yymsg) {
          yymsg = yymsgbuf;
          yymsg_alloc = sizeof yymsgbuf;
          yysyntax_error_status = 2;
        } else {
          yysyntax_error_status = YYSYNTAX_ERROR;
          yymsgp = yymsg;
        }
      }
      yyerror(input, molList, doQueries, scanner, yymsgp);
      if (yysyntax_error_status == 2) goto yyexhaustedlab;
    }
#undef YYSYNTAX_ERROR
#endif
  }

  if (yyerrstatus == 3) {
    /* If just tried and failed to reuse lookahead token after an
       error, discard it.  */

    if (yychar <= YYEOF) {
      /* Return failure if at end of input.  */
      if (yychar == YYEOF) YYABORT;
    } else {
      yydestruct("Error: discarding", yytoken, &yylval, input, molList,
                 doQueries, scanner);
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
  if (/*CONSTCOND*/ 0) goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK(yylen);
  yylen = 0;
  YY_STACK_PRINT(yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;

/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3; /* Each real token shifted decrements this.  */

  for (;;) {
    yyn = yypact[yystate];
    if (!yypact_value_is_default(yyn)) {
      yyn += YYTERROR;
      if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR) {
        yyn = yytable[yyn];
        if (0 < yyn) break;
      }
    }

    /* Pop the current state because it cannot handle the error token.  */
    if (yyssp == yyss) YYABORT;

    yydestruct("Error: popping", yystos[yystate], yyvsp, input, molList,
               doQueries, scanner);
    YYPOPSTACK(1);
    yystate = *yyssp;
    YY_STACK_PRINT(yyss, yyssp);
  }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Shift the error token.  */
  YY_SYMBOL_PRINT("Shifting", yystos[yyn], yyvsp, yylsp);

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
  yyerror(input, molList, doQueries, scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY) {
    /* Make sure we have latest lookahead translation.  See comments at
       user semantic actions for why this is necessary.  */
    yytoken = YYTRANSLATE(yychar);
    yydestruct("Cleanup: discarding lookahead", yytoken, &yylval, input,
               molList, doQueries, scanner);
  }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK(yylen);
  YY_STACK_PRINT(yyss, yyssp);
  while (yyssp != yyss) {
    yydestruct("Cleanup: popping", yystos[*yyssp], yyvsp, input, molList,
               doQueries, scanner);
    YYPOPSTACK(1);
  }
#ifndef yyoverflow
  if (yyss != yyssa) YYSTACK_FREE(yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf) YYSTACK_FREE(yymsg);
#endif
  return yyresult;
}
#line 632 "sln.yy" /* yacc.c:1910  */
