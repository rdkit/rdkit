//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <map>
#include "StructChecker.h"

namespace RDKit {
namespace StructureCheck {

static const char *flags[] = {
    "BAD_MOLECULE",
    "ALIAS_CONVERSION_FAILED",
    "STEREO_ERROR",
    "STEREO_FORCED_BAD",
    "ATOM_CLASH",
    "ATOM_CHECK_FAILED",
    "SIZE_CHECK_FAILED",
    "",  // reserved error = 0x0080,
    "TRANSFORMED",
    "FRAGMENTS_FOUND",
    "EITHER_WARNING",
    "DUBIOUS_STEREO_REMOVED",
    "RECHARGED",
    "STEREO_TRANSFORMED",
    "TEMPLATE_TRANSFORMED",
    "TAUTOMER_TRANSFORMED",
};

// Converts structure property flags to a comma separated string
std::string StructChecker::StructureFlagsToString(unsigned f) {
  std::string s;
  for (unsigned bit = 0; bit < 16; bit++) {
    if (0 != (f & (1 << bit))) {
      if (!s.empty()) s += ",";
      s += flags[bit];
    }
  }
  return s;
}

// Converts a comma separated string to a StructureFlag unsigned integer
class FMap : public std::map<std::string, unsigned> {
 public:
  FMap() {
    for (unsigned bit = 0; bit < 16; bit++)
      if (*flags[bit]) (*this)[std::string(flags[bit])] = (1 << bit);
  }
};

unsigned StructChecker::StringToStructureFlags(const std::string &str) {
  static const FMap fmap;  // map name string to StructureFlags enum value
  unsigned int f = 0;
  const char *token = str.c_str();
  while (*token) {
    while (*token && *token <= ' ')  // skip whitespaces (<tab>|<space>...)
      token++;
    unsigned len = 0;
    while (token[len] && !(token[len] == ',' || token[len] <= ' ')) len++;
    if (0 == len) continue;
    std::string name(token, len);
    auto it = fmap.find(name);
    if (fmap.end() != it) f |= it->second;
    while (token[len] &&
           (token[len] == ',' || token[len] <= ' '))  // skip delimeter
      len++;
    token += len;
  }
  return f;
  // there is no way to return syntax error in input string
}

}  // namespace StructureCheck
}  // namespace RDKit
