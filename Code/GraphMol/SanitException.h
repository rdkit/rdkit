//
//  Copyright (C) 2002-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RD_SANITEXCEPTION_H
#define RD_SANITEXCEPTION_H

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
  MolSanitizeException(const char *msg) : d_msg(msg){};
  MolSanitizeException(const std::string &msg) : d_msg(msg){};
  virtual const char *message() const { return d_msg.c_str(); };
  ~MolSanitizeException() throw(){};

 protected:
  std::string d_msg;
};

class RDKIT_GRAPHMOL_EXPORT AtomValenceException : public MolSanitizeException {
 public:
  AtomValenceException(const char *msg, unsigned int atomIdx)
      : MolSanitizeException(msg), d_atomIdx(atomIdx){};
  AtomValenceException(const std::string &msg, unsigned int atomIdx)
      : MolSanitizeException(msg), d_atomIdx(atomIdx){};
  unsigned int getAtomIdx() const { return d_atomIdx; };
  ~AtomValenceException() throw(){};

 protected:
  unsigned int d_atomIdx;
};

}  // namespace RDKit

#endif
