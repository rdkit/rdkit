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

#include "smiles.tab.hpp"
#include "SmilesScanner.h"

#undef yylex
#define yylex token_scanner.lex




#include "smiles.tab.hpp"




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

namespace RDKit {

  /// Build a parser object.
  SmilesParser::SmilesParser (RDKit::SmilesScanner& token_scanner_yyarg, std::vector<RDKit::RWMol *>& molList_yyarg, RDKit::Atom* &lastAtom_yyarg, RDKit::Bond* &lastBond_yyarg, unsigned &numBondsParsed_yyarg, std::list<unsigned int>& branchPoints_yyarg, int& start_token_yyarg)
#if YYDEBUG
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
#else
    :
#endif
      token_scanner (token_scanner_yyarg),
      molList (molList_yyarg),
      lastAtom (lastAtom_yyarg),
      lastBond (lastBond_yyarg),
      numBondsParsed (numBondsParsed_yyarg),
      branchPoints (branchPoints_yyarg),
      start_token (start_token_yyarg)
  {}

  SmilesParser::~SmilesParser ()
  {}

  SmilesParser::syntax_error::~syntax_error () YY_NOEXCEPT YY_NOTHROW
  {}

  /*---------.
  | symbol.  |
  `---------*/

  // basic_symbol.
  template <typename Base>
  SmilesParser::basic_symbol<Base>::basic_symbol (const basic_symbol& that)
    : Base (that)
    , value (that.value)
    , location (that.location)
  {}


  /// Constructor for valueless symbols.
  template <typename Base>
  SmilesParser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, YY_MOVE_REF (location_type) l)
    : Base (t)
    , value ()
    , location (l)
  {}

  template <typename Base>
  SmilesParser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, YY_RVREF (value_type) v, YY_RVREF (location_type) l)
    : Base (t)
    , value (YY_MOVE (v))
    , location (YY_MOVE (l))
  {}


  template <typename Base>
  SmilesParser::symbol_kind_type
  SmilesParser::basic_symbol<Base>::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }


  template <typename Base>
  bool
  SmilesParser::basic_symbol<Base>::empty () const YY_NOEXCEPT
  {
    return this->kind () == symbol_kind::S_YYEMPTY;
  }

  template <typename Base>
  void
  SmilesParser::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move (s);
    value = YY_MOVE (s.value);
    location = YY_MOVE (s.location);
  }

  // by_kind.
  SmilesParser::by_kind::by_kind () YY_NOEXCEPT
    : kind_ (symbol_kind::S_YYEMPTY)
  {}

#if 201103L <= YY_CPLUSPLUS
  SmilesParser::by_kind::by_kind (by_kind&& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {
    that.clear ();
  }
#endif

  SmilesParser::by_kind::by_kind (const by_kind& that) YY_NOEXCEPT
    : kind_ (that.kind_)
  {}

  SmilesParser::by_kind::by_kind (token_kind_type t) YY_NOEXCEPT
    : kind_ (yytranslate_ (t))
  {}



  void
  SmilesParser::by_kind::clear () YY_NOEXCEPT
  {
    kind_ = symbol_kind::S_YYEMPTY;
  }

  void
  SmilesParser::by_kind::move (by_kind& that)
  {
    kind_ = that.kind_;
    that.clear ();
  }

  SmilesParser::symbol_kind_type
  SmilesParser::by_kind::kind () const YY_NOEXCEPT
  {
    return kind_;
  }


  SmilesParser::symbol_kind_type
  SmilesParser::by_kind::type_get () const YY_NOEXCEPT
  {
    return this->kind ();
  }



  // by_state.
  SmilesParser::by_state::by_state () YY_NOEXCEPT
    : state (empty_state)
  {}

  SmilesParser::by_state::by_state (const by_state& that) YY_NOEXCEPT
    : state (that.state)
  {}

  void
  SmilesParser::by_state::clear () YY_NOEXCEPT
  {
    state = empty_state;
  }

  void
  SmilesParser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  SmilesParser::by_state::by_state (state_type s) YY_NOEXCEPT
    : state (s)
  {}

  SmilesParser::symbol_kind_type
  SmilesParser::by_state::kind () const YY_NOEXCEPT
  {
    if (state == empty_state)
      return symbol_kind::S_YYEMPTY;
    else
      return YY_CAST (symbol_kind_type, yystos_[+state]);
  }

  SmilesParser::stack_symbol_type::stack_symbol_type ()
  {}

  SmilesParser::stack_symbol_type::stack_symbol_type (YY_RVREF (stack_symbol_type) that)
    : super_type (YY_MOVE (that.state), YY_MOVE (that.value), YY_MOVE (that.location))
  {
#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  SmilesParser::stack_symbol_type::stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) that)
    : super_type (s, YY_MOVE (that.value), YY_MOVE (that.location))
  {
    // that is emptied.
    that.kind_ = symbol_kind::S_YYEMPTY;
  }

#if YY_CPLUSPLUS < 201103L
  SmilesParser::stack_symbol_type&
  SmilesParser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    value = that.value;
    location = that.location;
    return *this;
  }

  SmilesParser::stack_symbol_type&
  SmilesParser::stack_symbol_type::operator= (stack_symbol_type& that)
  {
    state = that.state;
    value = that.value;
    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void
  SmilesParser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);

    // User destructor.
    switch (yysym.kind ())
    {
      case symbol_kind::S_AROMATIC_ATOM_TOKEN: // AROMATIC_ATOM_TOKEN
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_ATOM_TOKEN: // ATOM_TOKEN
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_ORGANIC_ATOM_TOKEN: // ORGANIC_ATOM_TOKEN
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_BOND_TOKEN: // BOND_TOKEN
                    { delete (yysym.value.bond); }
        break;

      case symbol_kind::S_bondd: // bondd
                    { delete (yysym.value.bond); }
        break;

      case symbol_kind::S_atomd: // atomd
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_charge_element: // charge_element
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_h_element: // h_element
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_chiral_element: // chiral_element
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_element: // element
                    { delete (yysym.value.atom); }
        break;

      case symbol_kind::S_simple_atom: // simple_atom
                    { delete (yysym.value.atom); }
        break;

      default:
        break;
    }
  }

#if YYDEBUG
  template <typename Base>
  void
  SmilesParser::yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const
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
  SmilesParser::yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT (m, sym);
    yystack_.push (YY_MOVE (sym));
  }

  void
  SmilesParser::yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_ (m, stack_symbol_type (s, std::move (sym)));
#else
    stack_symbol_type ss (s, sym);
    yypush_ (m, ss);
#endif
  }

  void
  SmilesParser::yypop_ (int n) YY_NOEXCEPT
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  SmilesParser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  SmilesParser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  SmilesParser::debug_level_type
  SmilesParser::debug_level () const
  {
    return yydebug_;
  }

  void
  SmilesParser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  SmilesParser::state_type
  SmilesParser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - YYNTOKENS] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - YYNTOKENS];
  }

  bool
  SmilesParser::yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yypact_ninf_;
  }

  bool
  SmilesParser::yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yytable_ninf_;
  }

  int
  SmilesParser::operator() ()
  {
    return parse ();
  }

  int
  SmilesParser::parse ()
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
{ yyla.location.begin.column = 0; }



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
            yyla.kind_ = yytranslate_ (yylex (&yyla.value, &yyla.location, start_token));
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
      /* If YYLEN is nonzero, implement the default value of the
         action: '$$ = $1'.  Otherwise, use the top of the stack.

         Otherwise, the following line sets YYLHS.VALUE to garbage.
         This behavior is undocumented and Bison users should not rely
         upon it.  */
      if (yylen)
        yylhs.value = yystack_[yylen - 1].value;
      else
        yylhs.value = yystack_[0].value;

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
  case 2: // meta_start: START_MOL mol
              {
// the molList has already been updated, no need to do anything
}
    break;

  case 3: // meta_start: START_ATOM atomd EOS_TOKEN
                             {
  lastAtom = (yystack_[1].value.atom);
  YYACCEPT;
}
    break;

  case 4: // meta_start: START_ATOM bad_atom_def
                          {
  YYABORT;
}
    break;

  case 5: // meta_start: START_BOND bondd EOS_TOKEN
                             {
  lastBond = (yystack_[1].value.bond);
  YYACCEPT;
}
    break;

  case 6: // meta_start: START_BOND bondd
                   {
  delete (yystack_[0].value.bond);
  YYABORT;
}
    break;

  case 7: // meta_start: START_BOND
             {
  YYABORT;
}
    break;

  case 8: // meta_start: meta_start error EOS_TOKEN
                            {
  yyerrok;
  YYABORT;
}
    break;

  case 9: // meta_start: meta_start EOS_TOKEN
                       {
  YYACCEPT;
}
    break;

  case 10: // meta_start: error EOS_TOKEN
                  {
  yyerrok;
  YYABORT;
}
    break;

  case 14: // bad_atom_def: charge_element
                 {
  delete (yystack_[0].value.atom);
  YYABORT;
}
    break;

  case 15: // mol: atomd
           {
  int sz     = molList.size();
  molList.resize( sz + 1);
  molList[ sz ] = new RWMol();
  RDKit::RWMol *curMol = molList[ sz ];
  (yystack_[0].value.atom)->setProp(RDKit::common_properties::_SmilesStart,1);
  curMol->addAtom((yystack_[0].value.atom), true, true);
  //delete $1;
  (yylhs.value.moli) = sz;
}
    break;

  case 16: // mol: mol atomd
                  {
  RWMol *mp = molList[(yylhs.value.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yystack_[0].value.atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $2;
}
    break;

  case 17: // mol: mol BOND_TOKEN atomd
                        {
  RWMol *mp = molList[(yylhs.value.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yystack_[0].value.atom),true,true);
  if( (yystack_[1].value.bond)->getBondType() == Bond::DATIVER ){
    (yystack_[1].value.bond)->setBeginAtomIdx(atomIdx1);
    (yystack_[1].value.bond)->setEndAtomIdx(atomIdx2);
    (yystack_[1].value.bond)->setBondType(Bond::DATIVE);
  }else if ( (yystack_[1].value.bond)->getBondType() == Bond::DATIVEL ){
    (yystack_[1].value.bond)->setBeginAtomIdx(atomIdx2);
    (yystack_[1].value.bond)->setEndAtomIdx(atomIdx1);
    (yystack_[1].value.bond)->setBondType(Bond::DATIVE);
  } else {
    (yystack_[1].value.bond)->setBeginAtomIdx(atomIdx1);
    (yystack_[1].value.bond)->setEndAtomIdx(atomIdx2);
  }
  (yystack_[1].value.bond)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond((yystack_[1].value.bond),true);
  //delete $3;
}
    break;

  case 18: // mol: mol MINUS_TOKEN atomd
                        {
  RWMol *mp = molList[(yylhs.value.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yystack_[0].value.atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  //delete $3;
}
    break;

  case 19: // mol: mol SEPARATOR_TOKEN atomd
                            {
  RWMol *mp = molList[(yylhs.value.moli)];
  (yystack_[0].value.atom)->setProp(RDKit::common_properties::_SmilesStart,1,true);
  mp->addAtom((yystack_[0].value.atom),true,true);
}
    break;

  case 20: // mol: mol ring_number
                  {
  RWMol * mp = molList[(yylhs.value.moli)];
  Atom *atom=mp->getActiveAtom();
  mp->setAtomBookmark(atom,(yystack_[0].value.ival));

  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::UNSPECIFIED);
  mp->setBondBookmark(newB,(yystack_[0].value.ival));
  newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  if(!(mp->getAllBondsWithBookmark((yystack_[0].value.ival)).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yystack_[0].value.ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
}
    break;

  case 21: // mol: mol BOND_TOKEN ring_number
                             {
  RWMol * mp = molList[(yylhs.value.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     (yystack_[1].value.bond)->getBondType());
  if((yystack_[1].value.bond)->hasProp(RDKit::common_properties::_unspecifiedOrder)){
    newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
  }
  newB->setBondDir((yystack_[1].value.bond)->getBondDir());
  mp->setAtomBookmark(atom,(yystack_[0].value.ival));
  mp->setBondBookmark(newB,(yystack_[0].value.ival));
  if(!(mp->getAllBondsWithBookmark((yystack_[0].value.ival)).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yystack_[0].value.ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
  delete (yystack_[1].value.bond);
}
    break;

  case 22: // mol: mol MINUS_TOKEN ring_number
                              {
  RWMol * mp = molList[(yylhs.value.moli)];
  Atom *atom=mp->getActiveAtom();
  Bond *newB = mp->createPartialBond(atom->getIdx(),
				     Bond::SINGLE);
  mp->setAtomBookmark(atom,(yystack_[0].value.ival));
  mp->setBondBookmark(newB,(yystack_[0].value.ival));
  if(!(mp->getAllBondsWithBookmark((yystack_[0].value.ival)).size()%2)){
    newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
  }

  SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);

  INT_VECT tmp;
  atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
  tmp.push_back(-((yystack_[0].value.ival)+1));
  atom->setProp(RDKit::common_properties::_RingClosures,tmp);
}
    break;

  case 23: // mol: mol GROUP_OPEN_TOKEN atomd
                             {
  RWMol *mp = molList[(yylhs.value.moli)];
  Atom *a1 = mp->getActiveAtom();
  int atomIdx1=a1->getIdx();
  int atomIdx2=mp->addAtom((yystack_[0].value.atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,
	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  branchPoints.push_back(atomIdx1);
}
    break;

  case 24: // mol: mol GROUP_OPEN_TOKEN BOND_TOKEN atomd
                                         {
  RWMol *mp = molList[(yylhs.value.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yystack_[0].value.atom),true,true);
  if( (yystack_[1].value.bond)->getBondType() == Bond::DATIVER ){
    (yystack_[1].value.bond)->setBeginAtomIdx(atomIdx1);
    (yystack_[1].value.bond)->setEndAtomIdx(atomIdx2);
    (yystack_[1].value.bond)->setBondType(Bond::DATIVE);
  }else if ( (yystack_[1].value.bond)->getBondType() == Bond::DATIVEL ){
    (yystack_[1].value.bond)->setBeginAtomIdx(atomIdx2);
    (yystack_[1].value.bond)->setEndAtomIdx(atomIdx1);
    (yystack_[1].value.bond)->setBondType(Bond::DATIVE);
  } else {
    (yystack_[1].value.bond)->setBeginAtomIdx(atomIdx1);
    (yystack_[1].value.bond)->setEndAtomIdx(atomIdx2);
  }
  (yystack_[1].value.bond)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  mp->addBond((yystack_[1].value.bond),true);

  branchPoints.push_back(atomIdx1);
}
    break;

  case 25: // mol: mol GROUP_OPEN_TOKEN MINUS_TOKEN atomd
                                         {
  RWMol *mp = molList[(yylhs.value.moli)];
  int atomIdx1 = mp->getActiveAtom()->getIdx();
  int atomIdx2 = mp->addAtom((yystack_[0].value.atom),true,true);
  mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
  mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
  branchPoints.push_back(atomIdx1);
}
    break;

  case 26: // mol: mol GROUP_CLOSE_TOKEN
                        {
  if(branchPoints.empty()){
     SmilesParser::error(yystack_[0].location, "extra close parentheses");
     YYABORT;
  }
  RWMol *mp = molList[(yylhs.value.moli)];
  mp->setActiveAtom(branchPoints.back());
  branchPoints.pop_back();
}
    break;

  case 27: // bondd: BOND_TOKEN
            { (yylhs.value.bond) = (yystack_[0].value.bond); }
    break;

  case 28: // bondd: MINUS_TOKEN
                        {
          (yylhs.value.bond) = new Bond(Bond::SINGLE);
          }
    break;

  case 29: // atomd: simple_atom
        { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 30: // atomd: ATOM_OPEN_TOKEN charge_element COLON_TOKEN number ATOM_CLOSE_TOKEN
{
  (yylhs.value.atom) = (yystack_[3].value.atom);
  (yylhs.value.atom)->setNoImplicit(true);
  (yylhs.value.atom)->setProp(RDKit::common_properties::molAtomMapNumber,(yystack_[1].value.ival));
}
    break;

  case 31: // atomd: ATOM_OPEN_TOKEN charge_element ATOM_CLOSE_TOKEN
{
  (yylhs.value.atom) = (yystack_[1].value.atom);
  (yystack_[1].value.atom)->setNoImplicit(true);
}
    break;

  case 32: // charge_element: h_element
                { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 33: // charge_element: h_element PLUS_TOKEN
                       { (yystack_[1].value.atom)->setFormalCharge(1); }
    break;

  case 34: // charge_element: h_element PLUS_TOKEN PLUS_TOKEN
                                  { (yystack_[2].value.atom)->setFormalCharge(2); }
    break;

  case 35: // charge_element: h_element PLUS_TOKEN number
                              { (yystack_[2].value.atom)->setFormalCharge((yystack_[0].value.ival)); }
    break;

  case 36: // charge_element: h_element MINUS_TOKEN
                        { (yystack_[1].value.atom)->setFormalCharge(-1); }
    break;

  case 37: // charge_element: h_element MINUS_TOKEN MINUS_TOKEN
                                    { (yystack_[2].value.atom)->setFormalCharge(-2); }
    break;

  case 38: // charge_element: h_element MINUS_TOKEN number
                               { (yystack_[2].value.atom)->setFormalCharge(-(yystack_[0].value.ival)); }
    break;

  case 39: // h_element: H_TOKEN
                        { (yylhs.value.atom) = new Atom(1); }
    break;

  case 40: // h_element: number H_TOKEN
                                 { (yylhs.value.atom) = new Atom(1); (yylhs.value.atom)->setIsotope((yystack_[1].value.ival)); }
    break;

  case 41: // h_element: H_TOKEN H_TOKEN
                                  { (yylhs.value.atom) = new Atom(1); (yylhs.value.atom)->setNumExplicitHs(1); }
    break;

  case 42: // h_element: number H_TOKEN H_TOKEN
                                         { (yylhs.value.atom) = new Atom(1); (yylhs.value.atom)->setIsotope((yystack_[2].value.ival)); (yylhs.value.atom)->setNumExplicitHs(1);}
    break;

  case 43: // h_element: H_TOKEN H_TOKEN number
                                         { (yylhs.value.atom) = new Atom(1); (yylhs.value.atom)->setNumExplicitHs((yystack_[0].value.ival)); }
    break;

  case 44: // h_element: number H_TOKEN H_TOKEN number
                                                { (yylhs.value.atom) = new Atom(1); (yylhs.value.atom)->setIsotope((yystack_[3].value.ival)); (yylhs.value.atom)->setNumExplicitHs((yystack_[0].value.ival));}
    break;

  case 45: // h_element: chiral_element
                  { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 46: // h_element: chiral_element H_TOKEN
                                                        { (yylhs.value.atom) = (yystack_[1].value.atom); (yystack_[1].value.atom)->setNumExplicitHs(1);}
    break;

  case 47: // h_element: chiral_element H_TOKEN number
                                                { (yylhs.value.atom) = (yystack_[2].value.atom); (yystack_[2].value.atom)->setNumExplicitHs((yystack_[0].value.ival));}
    break;

  case 48: // chiral_element: element
                 { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 49: // chiral_element: element AT_TOKEN
                   { (yystack_[1].value.atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
    break;

  case 50: // chiral_element: element AT_TOKEN AT_TOKEN
                            { (yystack_[2].value.atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
    break;

  case 51: // chiral_element: element CHI_CLASS_TOKEN
                          { (yystack_[1].value.atom)->setChiralTag((yystack_[0].value.chiraltype)); (yystack_[1].value.atom)->setProp(common_properties::_chiralPermutation,0); }
    break;

  case 52: // chiral_element: element CHI_CLASS_TOKEN number
                                 { (yystack_[2].value.atom)->setChiralTag((yystack_[1].value.chiraltype)); (yystack_[2].value.atom)->setProp(common_properties::_chiralPermutation,(yystack_[0].value.ival)); }
    break;

  case 53: // element: simple_atom
                { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 54: // element: number simple_atom
                                           { (yystack_[0].value.atom)->setIsotope( (yystack_[1].value.ival) ); (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 55: // element: ATOM_TOKEN
                        { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 56: // element: number ATOM_TOKEN
                                                   { (yystack_[0].value.atom)->setIsotope( (yystack_[1].value.ival) ); (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 57: // element: HASH_TOKEN number
                                                 { (yylhs.value.atom) = new Atom((yystack_[0].value.ival)); }
    break;

  case 58: // element: number HASH_TOKEN number
                                                         { (yylhs.value.atom) = new Atom((yystack_[0].value.ival)); (yylhs.value.atom)->setIsotope((yystack_[2].value.ival)); }
    break;

  case 59: // simple_atom: ORGANIC_ATOM_TOKEN
                   { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 60: // simple_atom: AROMATIC_ATOM_TOKEN
                  { (yylhs.value.atom) = (yystack_[0].value.atom); }
    break;

  case 61: // ring_number: digit
              { (yylhs.value.ival) = (yystack_[0].value.ival); }
    break;

  case 62: // ring_number: PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit
                                          { (yylhs.value.ival) = (yystack_[1].value.ival)*10+(yystack_[0].value.ival); }
    break;

  case 63: // ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN
                                                         { (yylhs.value.ival) = (yystack_[1].value.ival); }
    break;

  case 64: // ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN
                                                               { (yylhs.value.ival) = (yystack_[2].value.ival)*10+(yystack_[1].value.ival); }
    break;

  case 65: // ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN
                                                                     { (yylhs.value.ival) = (yystack_[3].value.ival)*100+(yystack_[2].value.ival)*10+(yystack_[1].value.ival); }
    break;

  case 66: // ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN
                                                                           { (yylhs.value.ival) = (yystack_[4].value.ival)*1000+(yystack_[3].value.ival)*100+(yystack_[2].value.ival)*10+(yystack_[1].value.ival); }
    break;

  case 67: // ring_number: PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN
                                                                                 { (yylhs.value.ival) = (yystack_[5].value.ival)*10000+(yystack_[4].value.ival)*1000+(yystack_[3].value.ival)*100+(yystack_[2].value.ival)*10+(yystack_[1].value.ival); }
    break;

  case 68: // number: ZERO_TOKEN
         { (yylhs.value.ival) = (yystack_[0].value.ival); }
    break;

  case 69: // number: nonzero_number
  { (yylhs.value.ival) = (yystack_[0].value.ival); }
    break;

  case 70: // nonzero_number: NONZERO_DIGIT_TOKEN
                 { (yylhs.value.ival) = (yystack_[0].value.ival); }
    break;

  case 71: // nonzero_number: nonzero_number digit
                       {
  if((yystack_[1].value.ival) >= std::numeric_limits<std::int32_t>::max()/10 ||
     (yystack_[1].value.ival)*10 >= std::numeric_limits<std::int32_t>::max()-(yystack_[0].value.ival) ){
     SmilesParser::error(yystack_[1].location, "number too large");
     YYABORT;
  }
  (yylhs.value.ival) = (yystack_[1].value.ival)*10 + (yystack_[0].value.ival);
  }
    break;

  case 72: // digit: NONZERO_DIGIT_TOKEN
       { (yylhs.value.ival) = (yystack_[0].value.ival); }
    break;

  case 73: // digit: ZERO_TOKEN
  { (yylhs.value.ival) = (yystack_[0].value.ival); }
    break;



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
  SmilesParser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what ());
  }

#if YYDEBUG || 0
  const char *
  SmilesParser::symbol_name (symbol_kind_type yysymbol)
  {
    return yytname_[yysymbol];
  }
#endif // #if YYDEBUG || 0









  const signed char SmilesParser::yypact_ninf_ = -30;

  const signed char SmilesParser::yytable_ninf_ = -30;

  const short
  SmilesParser::yypact_[] =
  {
     135,   -25,    33,    70,   -12,     5,   -30,   -30,   -30,    96,
     110,   -30,   -30,   -30,   -30,   -30,   -10,    90,    13,    90,
      90,   -30,   -13,   -30,    77,    21,    14,     4,   120,    99,
     -30,   -30,    26,   -30,    47,   -30,    39,   -30,   -30,   -30,
      18,   -30,    33,    10,   137,    10,   -30,   -30,   -30,    13,
      90,   -30,   -30,   -30,    39,   -30,   -30,    52,    -1,    13,
      63,    13,   -30,    48,    13,   -30,   -30,   -30,   -30,    13,
     -30,    33,    33,   -30,   -30,   -30,   -30,    99,    99,   -30,
     -30,   -30,   -30,   -30,   -30,   -30,   -30,   -30,   -30,    13,
     -30,    64,   -30,   -30,   -30,    59,   -30,   -30,   -30,    76,
     -30,   121,   -30,   133,   -30,    72,   -30
  };

  const signed char
  SmilesParser::yydefact_[] =
  {
       0,     0,     0,     0,     7,     0,    10,    60,    59,     0,
       2,    15,    29,    55,    70,    68,    39,     0,     0,     0,
       0,     4,     0,    14,    32,    45,    48,    53,     0,    69,
      28,    27,     6,     1,     0,     9,     0,    53,    72,    73,
       0,    26,     0,     0,     0,     0,    16,    20,    61,    41,
       0,    13,    57,    11,    14,    12,     3,    36,    33,    46,
      49,    51,    56,    40,     0,    54,    71,     5,     8,     0,
      31,     0,     0,    23,    19,    18,    22,     0,     0,    17,
      21,    43,    37,    38,    34,    35,    47,    50,    52,    42,
      58,     0,    25,    24,    62,     0,    44,    30,    63,     0,
      64,     0,    65,     0,    66,     0,    67
  };

  const signed char
  SmilesParser::yypgoto_[] =
  {
     -30,   -30,    43,   -30,   -30,    11,    -7,   -30,   -30,   -30,
       8,   104,   -14,   -30,   -29
  };

  const signed char
  SmilesParser::yydefgoto_[] =
  {
       0,     5,    53,    10,    32,    11,    23,    24,    25,    26,
      12,    47,    28,    29,    48
  };

  const signed char
  SmilesParser::yytable_[] =
  {
      66,     6,    36,    30,    52,    33,    34,    49,    14,    15,
      31,    27,    54,    56,    22,    84,     7,    37,     8,    38,
      39,    46,    14,    15,     7,    37,     8,    37,    37,    44,
     -29,    35,    60,    71,     9,    81,    65,    61,    59,     7,
      72,     8,     9,    83,    85,    86,    21,    88,    94,    95,
      90,    73,    67,    74,    75,    91,    79,     9,    37,    69,
      51,    14,    15,    55,    70,    89,    99,    82,    38,    39,
     101,    98,   103,    68,   105,    96,     7,    13,     8,    14,
      15,    87,    92,    93,   106,    38,    39,    16,   100,    97,
      17,    18,    57,    58,    19,    20,     7,    13,     8,    14,
      15,     0,     7,    13,     8,    14,    15,    16,    38,    39,
      17,    18,     0,    16,    50,    20,     7,    18,     8,    38,
      39,    40,    41,    42,     0,    43,     7,    62,     8,    44,
      38,    39,    45,   102,     9,     0,     1,    63,     2,     3,
       4,    64,    38,    39,     0,   104,    77,    76,    78,    80
  };

  const signed char
  SmilesParser::yycheck_[] =
  {
      29,    26,     9,    15,    18,     0,     1,    17,     9,    10,
      22,     3,    19,    26,     3,    16,     6,     9,     8,     9,
      10,    10,     9,    10,     6,    17,     8,    19,    20,    19,
      26,    26,    18,    15,    24,    49,    28,    23,    17,     6,
      22,     8,    24,    57,    58,    59,     3,    61,    77,    78,
      64,    40,    26,    42,    43,    69,    45,    24,    50,    20,
      17,     9,    10,    20,    25,    17,    95,    15,     9,    10,
      99,    12,   101,    26,   103,    89,     6,     7,     8,     9,
      10,    18,    71,    72,    12,     9,    10,    17,    12,    25,
      20,    21,    15,    16,    24,    25,     6,     7,     8,     9,
      10,    -1,     6,     7,     8,     9,    10,    17,     9,    10,
      20,    21,    -1,    17,    24,    25,     6,    21,     8,     9,
      10,    11,    12,    13,    -1,    15,     6,     7,     8,    19,
       9,    10,    22,    12,    24,    -1,     1,    17,     3,     4,
       5,    21,     9,    10,    -1,    12,     9,    43,    11,    45
  };

  const signed char
  SmilesParser::yystos_[] =
  {
       0,     1,     3,     4,     5,    28,    26,     6,     8,    24,
      30,    32,    37,     7,     9,    10,    17,    20,    21,    24,
      25,    29,    32,    33,    34,    35,    36,    37,    39,    40,
      15,    22,    31,     0,     1,    26,    33,    37,     9,    10,
      11,    12,    13,    15,    19,    22,    32,    38,    41,    17,
      24,    29,    39,    29,    33,    29,    26,    15,    16,    17,
      18,    23,     7,    17,    21,    37,    41,    26,    26,    20,
      25,    15,    22,    32,    32,    32,    38,     9,    11,    32,
      38,    39,    15,    39,    16,    39,    39,    18,    39,    17,
      39,    39,    32,    32,    41,    41,    39,    25,    12,    41,
      12,    41,    12,    41,    12,    41,    12
  };

  const signed char
  SmilesParser::yyr1_[] =
  {
       0,    27,    28,    28,    28,    28,    28,    28,    28,    28,
      28,    29,    29,    29,    29,    30,    30,    30,    30,    30,
      30,    30,    30,    30,    30,    30,    30,    31,    31,    32,
      32,    32,    33,    33,    33,    33,    33,    33,    33,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    35,    35,
      35,    35,    35,    36,    36,    36,    36,    36,    36,    37,
      37,    38,    38,    38,    38,    38,    38,    38,    39,    39,
      40,    40,    41,    41
  };

  const signed char
  SmilesParser::yyr2_[] =
  {
       0,     2,     2,     3,     2,     3,     2,     1,     3,     2,
       2,     2,     2,     2,     1,     1,     2,     3,     3,     3,
       2,     3,     3,     3,     4,     4,     2,     1,     1,     1,
       5,     3,     1,     2,     3,     3,     2,     3,     3,     1,
       2,     2,     3,     3,     4,     1,     2,     3,     1,     2,
       3,     2,     3,     1,     2,     1,     2,     2,     3,     1,
       1,     1,     3,     4,     5,     6,     7,     8,     1,     1,
       1,     2,     1,     1
  };


#if YYDEBUG
  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a YYNTOKENS, nonterminals.
  const char*
  const SmilesParser::yytname_[] =
  {
  "\"end of file\"", "error", "\"invalid token\"", "START_MOL",
  "START_ATOM", "START_BOND", "AROMATIC_ATOM_TOKEN", "ATOM_TOKEN",
  "ORGANIC_ATOM_TOKEN", "NONZERO_DIGIT_TOKEN", "ZERO_TOKEN",
  "GROUP_OPEN_TOKEN", "GROUP_CLOSE_TOKEN", "SEPARATOR_TOKEN",
  "LOOP_CONNECTOR_TOKEN", "MINUS_TOKEN", "PLUS_TOKEN", "H_TOKEN",
  "AT_TOKEN", "PERCENT_TOKEN", "COLON_TOKEN", "HASH_TOKEN", "BOND_TOKEN",
  "CHI_CLASS_TOKEN", "ATOM_OPEN_TOKEN", "ATOM_CLOSE_TOKEN", "EOS_TOKEN",
  "$accept", "meta_start", "bad_atom_def", "mol", "bondd", "atomd",
  "charge_element", "h_element", "chiral_element", "element",
  "simple_atom", "ring_number", "number", "nonzero_number", "digit", YY_NULLPTR
  };
#endif


#if YYDEBUG
  const short
  SmilesParser::yyrline_[] =
  {
       0,    89,    89,    92,    96,    99,   103,   107,   110,   114,
     117,   124,   125,   126,   127,   135,   146,   157,   178,   187,
     193,   214,   238,   257,   267,   288,   296,   308,   309,   315,
     316,   322,   330,   331,   332,   333,   334,   335,   336,   340,
     341,   342,   343,   344,   345,   346,   347,   348,   352,   353,
     354,   355,   356,   360,   361,   362,   363,   364,   365,   369,
     370,   374,   375,   376,   377,   378,   379,   380,   384,   385,
     389,   390,   400,   401
  };

  void
  SmilesParser::yy_stack_print_ () const
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
  SmilesParser::yy_reduce_print_ (int yyrule) const
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

  SmilesParser::symbol_kind_type
  SmilesParser::yytranslate_ (int t) YY_NOEXCEPT
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
    // Last valid token kind.
    const int code_max = 281;

    if (t <= 0)
      return symbol_kind::S_YYEOF;
    else if (t <= code_max)
      return static_cast <symbol_kind_type> (translate_table[t]);
    else
      return symbol_kind::S_YYUNDEF;
  }

} // RDKit



static void free_allocated_mols(std::vector<RDKit::RWMol *>& molList) {
    for (auto& mol : molList) {
        if (mol) {
            SmilesParseOps::CleanupAfterParseError(mol);
            delete mol;
        }
    }
    molList.clear();
    molList.resize(0);
}

void RDKit::SmilesParser::error(const location& loc, const std::string& msg) {
    // we should try to free allocated memory in error conditions
    free_allocated_mols(molList);

    // we should attempt to point to the general area of the bad token
    auto bad_token_position = loc.begin.column;
    BOOST_LOG(rdErrorLog) << "SMILES Parse Error: " << msg << " while parsing input:" << std::endl;
    // NOTE: This may not be very useful for very long inputs
    BOOST_LOG(rdErrorLog) << token_scanner.input() << std::endl;
    BOOST_LOG(rdErrorLog) << std::string(bad_token_position - 1, '-') << '^' << std::endl;
}

