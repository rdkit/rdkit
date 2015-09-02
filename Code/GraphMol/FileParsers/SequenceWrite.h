//
//  Copyright (C) 2015 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SEQUENCEWRITE_H_
#define _RD_SEQUENCEWRITE_H_
#include <string>

namespace RDKit{
  class ROMol;

  std::string MolToSequence(const ROMol &mol);
  std::string MolToFASTA(const ROMol &mol);
  std::string MolToHELM(const ROMol &mol);
}

#endif

