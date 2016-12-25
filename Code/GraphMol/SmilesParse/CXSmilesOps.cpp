//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <iostream>
#include "SmilesParseOps.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>



namespace SmilesParseOps {
  using namespace RDKit;

  namespace parser {
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    template <typename Iterator>
    bool parse_it(Iterator first, Iterator last)
    {
      using qi::double_;
      using qi::phrase_parse;
      using ascii::space;
      std::cerr << "PARSE: " << last-first << std::endl;
      bool r = phrase_parse(
        first,                          /*< start iterator >*/
        last,                           /*< end iterator >*/
        ('|' >>
          *('(' >>
            -(double_) >> ',' >> -(double_) >> ',' >> -(double_) >>
            *( ';' >> -(double_) >> ',' >> -(double_) >> ',' >> -(double_) ) >>
          ')') >> 
          '|'
          ),   /*< the parser >*/
        space                           /*< the skip-parser >*/
      );
      std::cerr << "DONE: " << (last - first) << "? " << r << std::endl;
      if (first != last) // fail if we did not get a full match
        return false;
      return r;
    }

  } // end of namespace parser


  void parseCXNExtensions(RDKit::RWMol & mol, const std::string & extText) {
    std::cerr << "parseCXNExtensions: " << extText << std::endl;
    if (!extText.size() || extText[0] != '|') return;
    bool ok = parser::parse_it(extText.begin(), extText.end());
    POSTCONDITION(ok,"parse failed");



  }
} // end of namespace SmilesParseOps