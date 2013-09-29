//
//  Copyright (C) 2013 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_PDBPARSER_H_
#define _RD_PDBPARSER_H_
#include <GraphMol/FileParsers/FileParsers.h>

namespace RDKit {

  RWMol *PDBBlockToMol(const char *str, bool sanitize=true,
                       bool removeHs=true, unsigned int flavor=0);

  RWMol *PDBBlockToMol(const std::string &str, bool sanitize=true,
                       bool removeHs=true, unsigned int flavor=0);
  RWMol *PDBDataStreamToMol(std::istream *inStream, bool sanitize=true,
                            bool removeHs=true, unsigned int flavor=0);
  RWMol *PDBFileToMol(const std::string &fname, bool sanitize=true,
                      bool removeHs=true, unsigned int flavor=0);

}

#endif  // _RD_PDBPARSER_H_

