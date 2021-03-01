//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>

#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/format.hpp>
#include <algorithm>

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
namespace MolEnumerator {

struct LinkNode {
  unsigned int minRep = 0;
  unsigned int maxRep = 0;
  unsigned int nBonds = 0;
  std::vector<std::pair<unsigned int, unsigned int>> bondAtoms;
};

namespace utils {
inline std::vector<LinkNode> getMolLinkNodes(
    const ROMol &mol, bool strict = true,
    const std::map<unsigned, Atom *> *atomIdxMap = nullptr) {
  std::vector<LinkNode> res;
  std::string pval;
  if (!mol.getPropIfPresent(common_properties::molFileLinkNodes, pval)) {
    return res;
  }
  std::vector<int> mapping;

  boost::char_separator<char> pipesep("|");
  boost::char_separator<char> spacesep(" ");
  for (auto linknodetext : tokenizer(pval, pipesep)) {
    LinkNode node;
    tokenizer tokens(linknodetext, spacesep);
    std::vector<unsigned int> data;
    try {
      std::transform(tokens.begin(), tokens.end(), std::back_inserter(data),
                     [](const std::string &token) -> unsigned int {
                       return boost::lexical_cast<unsigned int>(token);
                     });
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Cannot convert values in LINKNODE '" << linknodetext
             << "' to unsigned ints";
      if (strict) {
        throw ValueErrorException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        continue;
      }
    }
    // the second test here is for the atom-pairs defining the bonds
    // data[2] contains the number of bonds
    if (data.size() < 5 || data.size() < 3 + 2 * data[2]) {
      std::ostringstream errout;
      errout << "not enough values in LINKNODE '" << linknodetext << "'";
      if (strict) {
        throw ValueErrorException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        continue;
      }
    }

    node.minRep = data[0];
    node.maxRep = data[1];
    if (node.minRep == 0 || node.maxRep < node.minRep) {
      std::ostringstream errout;
      errout << "bad counts in LINKNODE '" << linknodetext << "'";
      if (strict) {
        throw ValueErrorException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        continue;
      }
    }
    node.nBonds = data[2];
    if (node.nBonds != 2) {
      if (strict) {
        UNDER_CONSTRUCTION(
            "only link nodes with 2 bonds are currently supported");
      } else {
        BOOST_LOG(rdWarningLog)
            << "only link nodes with 2 bonds are currently supported"
            << std::endl;
        continue;
      }
    }
    // both bonds must start from the same atom:
    if (data[3] != data[5]) {
      std::ostringstream errout;
      errout << "bonds don't start at the same atom for LINKNODE '"
             << linknodetext << "'";
      if (strict) {
        throw ValueErrorException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        continue;
      }
    }

    if (atomIdxMap) {
      // map the indices back to the original atom numbers
      for (unsigned int i = 3; i <= 6; ++i) {
        const auto aidx = atomIdxMap->find(data[i] - 1);
        if (aidx == atomIdxMap->end()) {
          std::ostringstream errout;
          errout << "atom index " << data[i]
                 << " cannot be found in molecule for LINKNODE '"
                 << linknodetext << "'";
          if (strict) {
            throw ValueErrorException(errout.str());
          } else {
            BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
            continue;
          }
        } else {
          data[i] = aidx->second->getIdx();
        }
      }
    } else {
      for (unsigned int i = 3; i <= 6; ++i) {
        --data[i];
      }
    }
    node.bondAtoms.push_back(std::make_pair(data[3], data[4]));
    node.bondAtoms.push_back(std::make_pair(data[5], data[6]));
    if (!mol.getBondBetweenAtoms(data[4], data[3]) ||
        !mol.getBondBetweenAtoms(data[6], data[5])) {
      std::ostringstream errout;
      errout << "bond not found between atoms in LINKNODE '" << linknodetext
             << "'";
      if (strict) {
        throw ValueErrorException(errout.str());
      } else {
        BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
        continue;
      }
    }
    res.push_back(std::move(node));
  }
  return res;
}

}  // namespace utils
}  // namespace MolEnumerator

}  // namespace RDKit
