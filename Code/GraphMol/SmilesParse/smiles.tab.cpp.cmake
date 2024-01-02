// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.


// Take the name prefix into account.
#define yylex   yysmiles_lex

// First part of user prologue.
#line 1 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"


  // $Id$
  //
  //  Copyright (C) 2001-2016 Randal Henne, Greg Landrum and Rational Discovery LLC
  //
  //   @@ All Rights Reserved  @@
  //

#line 53 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"


#include "smiles.tab.hpp"


// Unqualified %code blocks.
#line 40 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"

#include <charconv>
#include <utility>

#include "SmilesScanner.h"
#include "SmilesASTParser.h"

#undef yylex
#define yylex token_scanner.lex


#line 72 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif


// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
# if defined __GNUC__ && !defined __EXCEPTIONS
#  define YY_EXCEPTIONS 0
# else
#  define YY_EXCEPTIONS 1
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (false)
# endif


// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << '\n';                       \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yy_stack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YY_USE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

#line 16 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
namespace smiles_parser { namespace ast_parser {
#line 165 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"

  /// Build a parser object.
  Parser::Parser (Scanner& token_scanner_yyarg, Builder& ast_yyarg)
#if YYDEBUG
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
#else
    :
#endif
      token_scanner (token_scanner_yyarg),
      ast (ast_yyarg)
  {}

  Parser::~Parser ()
  {}

  Parser::syntax_error::~syntax_error () YY_NOEXCEPT YY_NOTHROW
  {}

  /*---------.
  | symbol.  |
  `---------*/

  // basic_symbol.
  template <typename Base>
  Parser::basic_symbol<Base>::basic_symbol (const basic_symbol& that)
    : Base (that)
    , value ()
    , location (that.location)
  {
    switch (this->kind ())
    {
      case symbol_kind::S_plus_signs: // plus_signs
      case symbol_kind::S_minus_signs: // minus_signs
        value.copy< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_HYDROGEN: // HYDROGEN
      case symbol_kind::S_ORGANIC_ATOM: // ORGANIC_ATOM
      case symbol_kind::S_AROMATIC_SYMBOL: // AROMATIC_SYMBOL
      case symbol_kind::S_ELEMENT_SYMBOL: // ELEMENT_SYMBOL
      case symbol_kind::S_CHIRAL_TAG: // CHIRAL_TAG
      case symbol_kind::S_NUMBER: // NUMBER
        value.copy< std::string_view > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

  }




  template <typename Base>
  Parser::symbol_kind_type
  Parser::basic_symbol<Base>::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }


  template <typename Base>
  bool
  Parser::basic_symbol<Base>::empty () const YY_NOEXCEPT
  {
    return this->kind () == symbol_kind::S_YYEMPTY;
  }

  template <typename Base>
  void
  Parser::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move (s);
    switch (this->kind ())
    {
      case symbol_kind::S_plus_signs: // plus_signs
      case symbol_kind::S_minus_signs: // minus_signs
        value.move< int > (YY_MOVE (s.value));
        break;

      case symbol_kind::S_HYDROGEN: // HYDROGEN
      case symbol_kind::S_ORGANIC_ATOM: // ORGANIC_ATOM
      case symbol_kind::S_AROMATIC_SYMBOL: // AROMATIC_SYMBOL
      case symbol_kind::S_ELEMENT_SYMBOL: // ELEMENT_SYMBOL
      case symbol_kind::S_CHIRAL_TAG: // CHIRAL_TAG
      case symbol_kind::S_NUMBER: // NUMBER
        value.move< std::string_view > (YY_MOVE (s.value));
        break;

      default:
        break;
    }

    location = YY_MOVE (s.location);
  }

  // by_kind.
  Parser::by_kind::by_kind () YY_NOEXCEPT
    : kind_ (symbol_kind::S_YYEMPTY)
  {}

#if 201103L <= YY_CPLUSPLUS
  Parser::by_kind::by_kind (by_kind&& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {
    that.clear ();
  }
#endif

  Parser::by_kind::by_kind (const by_kind& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {}

  Parser::by_kind::by_kind (token_kind_type t) YY_NOEXCEPT
    : kind_ (yytranslate_ (t))
  {}



  void
  Parser::by_kind::clear () YY_NOEXCEPT
  {
    kind_ = symbol_kind::S_YYEMPTY;
  }

  void
  Parser::by_kind::move (by_kind& that)
  {
    kind_ = that.kind_;
    that.clear ();
  }

  Parser::symbol_kind_type
  Parser::by_kind::kind () const YY_NOEXCEPT
  {
    return kind_;
  }


  Parser::symbol_kind_type
  Parser::by_kind::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }



  // by_state.
  Parser::by_state::by_state () YY_NOEXCEPT
    : state (empty_state)
  {}

  Parser::by_state::by_state (const by_state& that) YY_NOEXCEPT
    : state (that.state)
  {}

  void
  Parser::by_state::clear () YY_NOEXCEPT
  {
    state = empty_state;
  }

  void
  Parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  Parser::by_state::by_state (state_type s) YY_NOEXCEPT
    : state (s)
  {}

  Parser::symbol_kind_type
  Parser::by_state::kind () const YY_NOEXCEPT
  {
    if (state == empty_state)
      return symbol_kind::S_YYEMPTY;
    else
      return YY_CAST (symbol_kind_type, yystos_[+state]);
  }

  Parser::stack_symbol_type::stack_symbol_type ()
  {}

  Parser::stack_symbol_type::stack_symbol_type (YY_RVREF (stack_symbol_type) that)
    : super_type (YY_MOVE (that.state), YY_MOVE (that.location))
  {
    switch (that.kind ())
    {
      case symbol_kind::S_plus_signs: // plus_signs
      case symbol_kind::S_minus_signs: // minus_signs
        value.YY_MOVE_OR_COPY< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_HYDROGEN: // HYDROGEN
      case symbol_kind::S_ORGANIC_ATOM: // ORGANIC_ATOM
      case symbol_kind::S_AROMATIC_SYMBOL: // AROMATIC_SYMBOL
      case symbol_kind::S_ELEMENT_SYMBOL: // ELEMENT_SYMBOL
      case symbol_kind::S_CHIRAL_TAG: // CHIRAL_TAG
      case symbol_kind::S_NUMBER: // NUMBER
        value.YY_MOVE_OR_COPY< std::string_view > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  Parser::stack_symbol_type::stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) that)
    : super_type (s, YY_MOVE (that.location))
  {
    switch (that.kind ())
    {
      case symbol_kind::S_plus_signs: // plus_signs
      case symbol_kind::S_minus_signs: // minus_signs
        value.move< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_HYDROGEN: // HYDROGEN
      case symbol_kind::S_ORGANIC_ATOM: // ORGANIC_ATOM
      case symbol_kind::S_AROMATIC_SYMBOL: // AROMATIC_SYMBOL
      case symbol_kind::S_ELEMENT_SYMBOL: // ELEMENT_SYMBOL
      case symbol_kind::S_CHIRAL_TAG: // CHIRAL_TAG
      case symbol_kind::S_NUMBER: // NUMBER
        value.move< std::string_view > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

    // that is emptied.
    that.kind_ = symbol_kind::S_YYEMPTY;
  }

#if YY_CPLUSPLUS < 201103L
  Parser::stack_symbol_type&
  Parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    switch (that.kind ())
    {
      case symbol_kind::S_plus_signs: // plus_signs
      case symbol_kind::S_minus_signs: // minus_signs
        value.copy< int > (that.value);
        break;

      case symbol_kind::S_HYDROGEN: // HYDROGEN
      case symbol_kind::S_ORGANIC_ATOM: // ORGANIC_ATOM
      case symbol_kind::S_AROMATIC_SYMBOL: // AROMATIC_SYMBOL
      case symbol_kind::S_ELEMENT_SYMBOL: // ELEMENT_SYMBOL
      case symbol_kind::S_CHIRAL_TAG: // CHIRAL_TAG
      case symbol_kind::S_NUMBER: // NUMBER
        value.copy< std::string_view > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    return *this;
  }

  Parser::stack_symbol_type&
  Parser::stack_symbol_type::operator= (stack_symbol_type& that)
  {
    state = that.state;
    switch (that.kind ())
    {
      case symbol_kind::S_plus_signs: // plus_signs
      case symbol_kind::S_minus_signs: // minus_signs
        value.move< int > (that.value);
        break;

      case symbol_kind::S_HYDROGEN: // HYDROGEN
      case symbol_kind::S_ORGANIC_ATOM: // ORGANIC_ATOM
      case symbol_kind::S_AROMATIC_SYMBOL: // AROMATIC_SYMBOL
      case symbol_kind::S_ELEMENT_SYMBOL: // ELEMENT_SYMBOL
      case symbol_kind::S_CHIRAL_TAG: // CHIRAL_TAG
      case symbol_kind::S_NUMBER: // NUMBER
        value.move< std::string_view > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void
  Parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);
  }

#if YYDEBUG
  template <typename Base>
  void
  Parser::yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YY_USE (yyoutput);
    if (yysym.empty ())
      yyo << "empty symbol";
    else
      {
        symbol_kind_type yykind = yysym.kind ();
        yyo << (yykind < YYNTOKENS ? "token" : "nterm")
            << ' ' << yysym.name () << " ("
            << yysym.location << ": ";
        YY_USE (yykind);
        yyo << ')';
      }
  }
#endif

  void
  Parser::yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT (m, sym);
    yystack_.push (YY_MOVE (sym));
  }

  void
  Parser::yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_ (m, stack_symbol_type (s, std::move (sym)));
#else
    stack_symbol_type ss (s, sym);
    yypush_ (m, ss);
#endif
  }

  void
  Parser::yypop_ (int n) YY_NOEXCEPT
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  Parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  Parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  Parser::debug_level_type
  Parser::debug_level () const
  {
    return yydebug_;
  }

  void
  Parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  Parser::state_type
  Parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - YYNTOKENS] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - YYNTOKENS];
  }

  bool
  Parser::yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yypact_ninf_;
  }

  bool
  Parser::yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yytable_ninf_;
  }

  int
  Parser::operator() ()
  {
    return parse ();
  }

  int
  Parser::parse ()
  {
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
      {
    YYCDEBUG << "Starting parse\n";


    // User initialization code.
#line 38 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
{ yyla.location.begin.column = 0; }

#line 609 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, YY_MOVE (yyla));

  /*-----------------------------------------------.
  | yynewstate -- push a new symbol on the stack.  |
  `-----------------------------------------------*/
  yynewstate:
    YYCDEBUG << "Entering state " << int (yystack_[0].state) << '\n';
    YY_STACK_PRINT ();

    // Accept?
    if (yystack_[0].state == yyfinal_)
      YYACCEPT;

    goto yybackup;


  /*-----------.
  | yybackup.  |
  `-----------*/
  yybackup:
    // Try to take a decision without lookahead.
    yyn = yypact_[+yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token\n";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
          {
            yyla.kind_ = yytranslate_ (yylex (&yyla.value, &yyla.location));
          }
#if YY_EXCEPTIONS
        catch (const syntax_error& yyexc)
          {
            YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
            error (yyexc);
            goto yyerrlab1;
          }
#endif // YY_EXCEPTIONS
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    if (yyla.kind () == symbol_kind::S_YYerror)
    {
      // The scanner already issued an error message, process directly
      // to error recovery.  But do not keep the error token as
      // lookahead, it is too special and may lead us to an endless
      // loop in error recovery. */
      yyla.kind_ = symbol_kind::S_YYUNDEF;
      goto yyerrlab1;
    }

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.kind ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.kind ())
      {
        goto yydefault;
      }

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", state_type (yyn), YY_MOVE (yyla));
    goto yynewstate;


  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[+yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;


  /*-----------------------------.
  | yyreduce -- do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_ (yystack_[yylen].state, yyr1_[yyn]);
      /* Variants are always initialized to an empty instance of the
         correct type. The default '$$ = $1' action is NOT applied
         when using variants.  */
      switch (yyr1_[yyn])
    {
      case symbol_kind::S_plus_signs: // plus_signs
      case symbol_kind::S_minus_signs: // minus_signs
        yylhs.value.emplace< int > ();
        break;

      case symbol_kind::S_HYDROGEN: // HYDROGEN
      case symbol_kind::S_ORGANIC_ATOM: // ORGANIC_ATOM
      case symbol_kind::S_AROMATIC_SYMBOL: // AROMATIC_SYMBOL
      case symbol_kind::S_ELEMENT_SYMBOL: // ELEMENT_SYMBOL
      case symbol_kind::S_CHIRAL_TAG: // CHIRAL_TAG
      case symbol_kind::S_NUMBER: // NUMBER
        yylhs.value.emplace< std::string_view > ();
        break;

      default:
        break;
    }


      // Default location.
      {
        stack_type::slice range (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, range, yylen);
        yyerror_range[1].location = yylhs.location;
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
#if YY_EXCEPTIONS
      try
#endif // YY_EXCEPTIONS
        {
          switch (yyn)
            {
  case 11: // dot: '.'
#line 69 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
         { ast.add_dot(); }
#line 760 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 12: // branch_open: '('
#line 70 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                  { ast.open_branch(); }
#line 766 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 13: // branch_close: ')'
#line 71 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                   { ast.close_branch(); }
#line 772 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 14: // bond: '-'
#line 73 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
          { ast.add_bond("-"); }
#line 778 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 15: // bond: '='
#line 74 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
          { ast.add_bond("="); }
#line 784 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 16: // bond: '#'
#line 75 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
          { ast.add_bond("#"); }
#line 790 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 17: // bond: ':'
#line 76 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
          { ast.add_bond(":"); }
#line 796 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 18: // bond: '$'
#line 77 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
          { ast.add_bond("$"); }
#line 802 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 19: // bond: '~'
#line 78 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
          { ast.add_bond("~"); }
#line 808 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 20: // bond: '/'
#line 79 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
          { ast.add_bond("/"); }
#line 814 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 21: // bond: '\\'
#line 80 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
              { ast.add_bond("\\"); }
#line 820 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 22: // bond: '-' '>'
#line 81 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
              { ast.add_bond("->"); }
#line 826 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 23: // bond: '<' '-'
#line 82 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
              { ast.add_bond("<-"); }
#line 832 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 24: // bond: '\\' '\\'
#line 83 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                { ast.add_bond("\\\\"); }
#line 838 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 27: // atom: bracket_atom
#line 85 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                                          { ast.set_no_implicit_hs(); }
#line 844 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 31: // isotope: NUMBER symbol
#line 90 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                                { ast.add_isotope_num(yystack_[1].value.as < std::string_view > ()); }
#line 850 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 35: // organic: ORGANIC_ATOM
#line 94 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                      { ast.add_atom(yystack_[0].value.as < std::string_view > ()); }
#line 856 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 36: // element_symbol: ELEMENT_SYMBOL
#line 96 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                               { ast.add_atom(yystack_[0].value.as < std::string_view > ()); }
#line 862 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 37: // element_symbol: HYDROGEN
#line 97 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                         { ast.add_atom(yystack_[0].value.as < std::string_view > ()); }
#line 868 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 38: // element_symbol: '\'' ELEMENT_SYMBOL '\''
#line 99 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                                         { ast.add_atom(yystack_[1].value.as < std::string_view > ()); }
#line 874 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 39: // element_symbol: AROMATIC_SYMBOL
#line 100 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                                 { ast.add_atom(yystack_[0].value.as < std::string_view > ()); }
#line 880 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 40: // element_symbol: '#' NUMBER
#line 101 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                           { ast.add_atom(yystack_[0].value.as < std::string_view > ()); }
#line 886 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 41: // dummy_atom: '*'
#line 104 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                { ast.add_atom("*"); }
#line 892 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 43: // chirality: '@'
#line 107 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
               { ast.add_chirality_tag("@"); }
#line 898 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 44: // chirality: '@' '@'
#line 108 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                   { ast.add_chirality_tag("@@"); }
#line 904 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 45: // chirality: CHIRAL_TAG
#line 109 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                      { ast.add_chirality_class(yystack_[0].value.as < std::string_view > ()); }
#line 910 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 46: // chirality: CHIRAL_TAG NUMBER
#line 110 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                             { ast.add_chirality_class(yystack_[1].value.as < std::string_view > (), yystack_[0].value.as < std::string_view > ()); }
#line 916 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 48: // hcount: HYDROGEN
#line 114 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                 { ast.add_explicit_h(1); }
#line 922 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 49: // hcount: HYDROGEN NUMBER
#line 115 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                        { ast.add_explicit_h(yystack_[0].value.as < std::string_view > ()); }
#line 928 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 51: // charge: plus_signs
#line 118 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                   { ast.add_atom_charge(yystack_[0].value.as < int > ()); }
#line 934 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 52: // charge: minus_signs
#line 119 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                    { ast.add_atom_charge(yystack_[0].value.as < int > ()); }
#line 940 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 53: // charge: '+' NUMBER
#line 120 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                   { ast.add_atom_charge(yystack_[0].value.as < std::string_view > ()); }
#line 946 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 54: // charge: '-' NUMBER
#line 121 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                   { ast.add_atom_charge(yystack_[0].value.as < std::string_view > ()); }
#line 952 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 55: // plus_signs: '+'
#line 123 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                { yylhs.value.as < int > () = 1; }
#line 958 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 56: // plus_signs: plus_signs '+'
#line 123 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                                             { yylhs.value.as < int > () = yystack_[1].value.as < int > () + 1; }
#line 964 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 57: // minus_signs: '-'
#line 124 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                 { yylhs.value.as < int > () = -1; }
#line 970 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 58: // minus_signs: minus_signs '-'
#line 124 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                                                { yylhs.value.as < int > () = yystack_[1].value.as < int > () - 1; }
#line 976 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 59: // map_number: NUMBER
#line 126 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                   { ast.add_atom_map_number(yystack_[0].value.as < std::string_view > ()); }
#line 982 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 62: // ring_number: NUMBER
#line 130 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                    { ast.add_ring(yystack_[0].value.as < std::string_view > ()); }
#line 988 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 63: // ring_number: '%' NUMBER
#line 131 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                        { ast.add_ring(yystack_[0].value.as < std::string_view > ()); }
#line 994 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;

  case 64: // ring_number: '%' '(' NUMBER ')'
#line 132 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
                                { ast.add_ring(yystack_[1].value.as < std::string_view > ()); }
#line 1000 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"
    break;


#line 1004 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"

            default:
              break;
            }
        }
#if YY_EXCEPTIONS
      catch (const syntax_error& yyexc)
        {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error (yyexc);
          YYERROR;
        }
#endif // YY_EXCEPTIONS
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, YY_MOVE (yylhs));
    }
    goto yynewstate;


  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        std::string msg = YY_("syntax error");
        error (yyla.location, YY_MOVE (msg));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.kind () == symbol_kind::S_YYEOF)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:
    /* Pacify compilers when the user code never invokes YYERROR and
       the label yyerrorlab therefore never appears in user code.  */
    if (false)
      YYERROR;

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    YY_STACK_PRINT ();
    goto yyerrlab1;


  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    // Pop stack until we find a state that shifts the error token.
    for (;;)
      {
        yyn = yypact_[+yystack_[0].state];
        if (!yy_pact_value_is_default_ (yyn))
          {
            yyn += symbol_kind::S_YYerror;
            if (0 <= yyn && yyn <= yylast_
                && yycheck_[yyn] == symbol_kind::S_YYerror)
              {
                yyn = yytable_[yyn];
                if (0 < yyn)
                  break;
              }
          }

        // Pop the current state because it cannot handle the error token.
        if (yystack_.size () == 1)
          YYABORT;

        yyerror_range[1].location = yystack_[0].location;
        yy_destroy_ ("Error: popping", yystack_[0]);
        yypop_ ();
        YY_STACK_PRINT ();
      }
    {
      stack_symbol_type error_token;

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = state_type (yyn);
      yypush_ ("Shifting", YY_MOVE (error_token));
    }
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


  /*-----------------------------------------------------.
  | yyreturn -- parsing is finished, return the result.  |
  `-----------------------------------------------------*/
  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    YY_STACK_PRINT ();
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
#if YY_EXCEPTIONS
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
        // Do not try to display the values of the reclaimed symbols,
        // as their printers might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
#endif // YY_EXCEPTIONS
  }

  void
  Parser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what ());
  }

#if YYDEBUG || 0
  const char *
  Parser::symbol_name (symbol_kind_type yysymbol)
  {
    return yytname_[yysymbol];
  }
#endif // #if YYDEBUG || 0









  const signed char Parser::yypact_ninf_ = -14;

  const signed char Parser::yytable_ninf_ = -4;

  const signed char
  Parser::yypact_[] =
  {
      16,   -14,    19,   -14,     1,     0,   -14,   -14,   -14,   -14,
     -14,   -14,   -14,    26,    -3,    31,    38,   -14,   -14,   -14,
     -14,   -14,   -14,   -14,   -14,   -13,   -14,   -14,   -14,   -14,
     -14,   -14,    20,    30,    -6,    16,   -14,    61,    -2,   -14,
     -14,   -14,   -14,    24,    58,    41,    65,   -14,   -14,   -14,
     -14,    63,   -14,    16,   -14,   -14,   -14,   -14,   -14,    73,
      22,    72,    44,   -14,    76,    77,   -12,    42,    74,   -14,
     -14,   -14,   -14,   -14,    79,   -14,   -14,   -14,   -14,    66,
     -14
  };

  const signed char
  Parser::yydefact_[] =
  {
       2,    35,     0,    41,     0,     8,     4,    27,    25,    26,
      37,    39,    36,     0,     0,     0,    42,    30,    32,    33,
      34,     1,    62,    11,    12,    14,    15,    16,    17,    18,
      19,    20,    21,     0,     0,     0,     9,     8,    10,     6,
      60,    31,    40,     0,    45,    43,    47,    22,    24,    23,
      63,     0,     5,     0,    10,    61,    38,    46,    44,    48,
      50,     0,     8,    49,    57,    55,     0,    51,    52,    64,
      13,     7,    54,    53,     0,    28,    56,    58,    59,     0,
      29
  };

  const signed char
  Parser::yypgoto_[] =
  {
     -14,   -14,    35,    53,   -14,   -14,   -14,    54,    57,   -14,
     -14,    80,    33,   -14,    34,   -14,   -14,   -14,   -14,   -14,
     -14,   -14,    56
  };

  const signed char
  Parser::yydefgoto_[] =
  {
       0,     4,     5,    35,    36,    37,    71,    38,     6,     7,
      16,    17,     8,    19,     9,    46,    60,    66,    67,    68,
      79,    39,    40
  };

  const signed char
  Parser::yytable_[] =
  {
      -3,    21,    50,    74,    51,    42,    22,    47,    22,    23,
      24,    75,    25,    26,    27,    28,    29,    30,    31,    32,
       1,    33,    10,     1,    11,    12,    34,    13,    34,    10,
       1,    11,    12,    14,    64,    18,    20,    43,     2,    48,
      14,     3,    49,    15,     3,    44,    18,    20,    56,    65,
      15,     3,    22,    23,    24,    70,    25,    26,    27,    28,
      29,    30,    31,    32,    45,    33,    57,    58,    59,    76,
      23,    61,    34,    25,    26,    27,    28,    29,    30,    31,
      32,    63,    33,    69,    72,    73,    77,    78,    62,    80,
      53,    54,    52,    41,    55
  };

  const signed char
  Parser::yycheck_[] =
  {
       0,     0,     8,    15,    10,     8,     8,    20,     8,     9,
      10,    23,    12,    13,    14,    15,    16,    17,    18,    19,
       4,    21,     3,     4,     5,     6,    28,     8,    28,     3,
       4,     5,     6,    14,    12,     2,     2,     6,    22,    19,
      14,    25,    12,    24,    25,     7,    13,    13,    24,    27,
      24,    25,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    26,    21,     8,    26,     3,    27,
       9,     8,    28,    12,    13,    14,    15,    16,    17,    18,
      19,     8,    21,    11,     8,     8,    12,     8,    53,    23,
      37,    37,    35,    13,    38
  };

  const signed char
  Parser::yystos_[] =
  {
       0,     4,    22,    25,    30,    31,    37,    38,    41,    43,
       3,     5,     6,     8,    14,    24,    39,    40,    41,    42,
      43,     0,     8,     9,    10,    12,    13,    14,    15,    16,
      17,    18,    19,    21,    28,    32,    33,    34,    36,    50,
      51,    40,     8,     6,     7,    26,    44,    20,    19,    12,
       8,    10,    37,    32,    36,    51,    24,     8,    26,     3,
      45,     8,    31,     8,    12,    27,    46,    47,    48,    11,
      11,    35,     8,     8,    15,    23,    27,    12,     8,    49,
      23
  };

  const signed char
  Parser::yyr1_[] =
  {
       0,    29,    30,    30,    31,    31,    31,    31,    32,    32,
      32,    33,    34,    35,    36,    36,    36,    36,    36,    36,
      36,    36,    36,    36,    36,    37,    37,    37,    38,    38,
      39,    39,    40,    40,    40,    41,    42,    42,    42,    42,
      42,    43,    44,    44,    44,    44,    44,    45,    45,    45,
      46,    46,    46,    46,    46,    47,    47,    48,    48,    49,
      50,    50,    51,    51,    51
  };

  const signed char
  Parser::yyr2_[] =
  {
       0,     2,     0,     1,     1,     3,     2,     5,     0,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     2,     2,     2,     1,     1,     1,     6,     8,
       1,     2,     1,     1,     1,     1,     1,     1,     3,     1,
       2,     1,     0,     1,     2,     1,     2,     0,     1,     2,
       0,     1,     1,     2,     2,     1,     2,     1,     2,     1,
       1,     2,     1,     2,     4
  };


#if YYDEBUG
  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
  const char*
  const Parser::yytname_[] =
  {
  "\"end of file\"", "error", "\"invalid token\"", "HYDROGEN",
  "ORGANIC_ATOM", "AROMATIC_SYMBOL", "ELEMENT_SYMBOL", "CHIRAL_TAG",
  "NUMBER", "'.'", "'('", "')'", "'-'", "'='", "'#'", "':'", "'$'", "'~'",
  "'/'", "'\\\\'", "'>'", "'<'", "'['", "']'", "'\\''", "'*'", "'@'",
  "'+'", "'%'", "$accept", "mol", "chain", "connector", "dot",
  "branch_open", "branch_close", "bond", "atom", "bracket_atom", "isotope",
  "symbol", "organic", "element_symbol", "dummy_atom", "chirality",
  "hcount", "charge", "plus_signs", "minus_signs", "map_number",
  "ring_bond", "ring_number", YY_NULLPTR
  };
#endif


#if YYDEBUG
  const unsigned char
  Parser::yyrline_[] =
  {
       0,    60,    60,    60,    62,    63,    64,    65,    67,    67,
      67,    69,    70,    71,    73,    74,    75,    76,    77,    78,
      79,    80,    81,    82,    83,    85,    85,    85,    87,    88,
      90,    90,    92,    92,    92,    94,    96,    97,    99,   100,
     101,   104,   106,   107,   108,   109,   110,   113,   114,   115,
     117,   118,   119,   120,   121,   123,   123,   124,   124,   126,
     128,   128,   130,   131,   132
  };

  void
  Parser::yy_stack_print_ () const
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << int (i->state);
    *yycdebug_ << '\n';
  }

  void
  Parser::yy_reduce_print_ (int yyrule) const
  {
    int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG

  Parser::symbol_kind_type
  Parser::yytranslate_ (int t) YY_NOEXCEPT
  {
    // YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to
    // TOKEN-NUM as returned by yylex.
    static
    const signed char
    translate_table[] =
    {
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,    14,    16,    28,     2,    24,
      10,    11,    25,    27,     2,    12,     9,    18,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    15,     2,
      21,    13,    20,     2,    26,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    22,    19,    23,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,    17,     2,     2,     2,
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
       5,     6,     7,     8
    };
    // Last valid token kind.
    const int code_max = 263;

    if (t <= 0)
      return symbol_kind::S_YYEOF;
    else if (t <= code_max)
      return static_cast <symbol_kind_type> (translate_table[t]);
    else
      return symbol_kind::S_YYUNDEF;
  }

#line 16 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"
} } // smiles_parser::ast_parser
#line 1423 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.tab.cpp"

#line 133 "/Users/faara/Documents/code/rdkit_builder/rdkit/Code/GraphMol/SmilesParse/smiles.yy"


void smiles_parser::ast_parser::Parser::error(const location& loc, const std::string& msg) {
    auto bad_token_position = loc.begin.column - 1;
    ast.save_error_message(bad_token_position, msg);
}
