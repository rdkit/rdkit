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
#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>



namespace SmilesParseOps {
  using namespace RDKit;

  namespace parser {
    template <typename Iterator>
    std::string read_text_to(Iterator &first, Iterator last, const char sep, const char blockend) {
      Iterator start = first; 
      // EFF: there are certainly faster ways to do this
      while (first != last && *first != sep && *first != blockend) {
        ++first;
      }
      std::string res(start, first);
      return res;
    }

    template <typename Iterator>
    bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol) {
      if (first >= last || *first != '$') return false;
      ++first;
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

      RDKit::Conformer *conf = new Conformer(mol.getNumAtoms());
      mol.addConformer(conf);
      ++first;
      unsigned int atIdx = 0;
      while (first != last && *first != ')') {
        RDGeom::Point3D pt;
        std::string tkn = read_text_to(first, last, ';', ')');
        if (tkn != "") {
          std::vector<std::string> tokens;
          boost::split(tokens, tkn, boost::is_any_of(","));
          if(tokens.size()>=1 && tokens[0].size())
             pt.x = boost::lexical_cast<double>(tokens[0]);
          if (tokens.size() >= 2 && tokens[1].size())
            pt.y = boost::lexical_cast<double>(tokens[1]);
          if (tokens.size() >= 3 && tokens[2].size())
            pt.z = boost::lexical_cast<double>(tokens[2]);
        }

        conf->setAtomPos(atIdx, pt);
        ++atIdx;
        if (first != last && *first != ')') ++first;
      }
      if (first == last || *first != ')') return false;
      return true;

    }
    template <typename Iterator>
    bool parse_it(Iterator &first, Iterator last,RDKit::RWMol &mol) {
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
      ++first; // step past the last '|'
      return true;
    }
  } // end of namespace parser

  void parseCXExtensions(RDKit::RWMol & mol, const std::string & extText, std::string::const_iterator &first) {
    //std::cerr << "parseCXNExtensions: " << extText << std::endl;
    if (!extText.size() || extText[0] != '|') return;
    first = extText.begin();
    bool ok = parser::parse_it(first, extText.end(),mol);
    POSTCONDITION(ok,"parse failed");
  }
} // end of namespace SmilesParseOps