//
//  Copyright (C) 2015 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SEQUENCEPARSE_H_
#define _RD_SEQUENCEPARSE_H_
#include <string>

namespace RDKit{
  class RWMol;

  RWMol *SequenceToMol(const char *seq, bool sanitize=true, bool lowerD=false);
  RWMol *SequenceToMol(const std::string &seq, bool sanitize=true,
                       bool lowerD=false);

  RWMol *FASTAToMol(const char *seq, bool sanitize=true, bool lowerD=false);
  RWMol *FASTAToMol(const std::string &seq, bool sanitize=true,
                    bool lowerD=false);

  RWMol *HELMToMol(const char *helm, bool sanitize=true);
  RWMol *HELMToMol(const std::string &helm, bool sanitize=true);
}

#endif

