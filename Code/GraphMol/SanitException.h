//
//  Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef _RD_SANITEXCEPTION_H
#define _RD_SANITEXCEPTION_H

#include <RDGeneral/types.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

#include <string>
#include <vector>
#include <exception>

namespace RDKit {

//! class for flagging sanitization errors
class RDKIT_GRAPHMOL_EXPORT MolSanitizeException : public std::exception {
 public:
  MolSanitizeException(const char *msg) : _msg(msg){};
  MolSanitizeException(const std::string &msg) : _msg(msg){};
  const char *message() const { return _msg.c_str(); };
  ~MolSanitizeException() throw(){};

 private:
  std::string _msg;
};
}

#endif
