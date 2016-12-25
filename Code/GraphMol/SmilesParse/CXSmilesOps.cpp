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
//#include <boost/config/warning_disable.hpp>
//#include <boost/spirit/include/qi.hpp>



namespace SmilesParseOps {
  using namespace RDKit;

  namespace parser {
    template <typename Iterator>
    std::string read_text_to(Iterator &first, Iterator last, const char sep, const char blockend) {
      // EFF: there are certainly faster ways to do this
      std::string res = "";
      while (first != last && *first != sep && *first != blockend) {
        res += *first;
        ++first;
      }
      std::cerr << "   tkn: " << res << std::endl;
      return res;
    }

    template <typename Iterator>
    bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol) {
      if (first >= last || *first != '$') return false;
      ++first;
      std::cerr << "atom labels" << std::endl;
      unsigned int atIdx = 0;
      while (first != last && *first != '$') {
        std::string tkn = read_text_to(first, last, ';', '$');
        if (tkn != "") {
          mol.getAtomWithIdx(atIdx)->setProp("_atomLabel", tkn);
        }
        ++atIdx;
        if(first != last && *first!='$') ++first;
      }
      if (first == last || *first != '$') return false;
      return true;
    }

    template <typename Iterator>
    bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol) {
      if (first >= last || *first != '(') return false;
      ++first;
      std::cerr << "atom coords" << std::endl;
      while (first != last && *first != ')') {
        ++first;
      }
      if (first == last || *first != ')') return false;
      return true;

    }
    template <typename Iterator>
    bool parse_it(Iterator first, Iterator last,RDKit::RWMol &mol) {
      if (first >= last || *first != '|') return false;
      ++first;
      while (first < last && *first != '|') {
        if (*first == '(') {
          if (!parse_coords(first, last, mol)) return false;
        }
        else if (*first == '$') {
          if (!parse_atom_labels(first, last, mol)) return false;
        }
        ++first;
      }
      if (first >= last || *first != '|') return false;
      return true;
    }
  } // end of namespace parser

  


  void parseCXNExtensions(RDKit::RWMol & mol, const std::string & extText) {
    std::cerr << "parseCXNExtensions: " << extText << std::endl;
    if (!extText.size() || extText[0] != '|') return;
    bool ok = parser::parse_it(extText.begin(), extText.end(),mol);
    POSTCONDITION(ok,"parse failed");



  }
} // end of namespace SmilesParseOps