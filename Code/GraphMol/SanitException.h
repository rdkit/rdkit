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
  MolSanitizeException(const MolSanitizeException &other)
      : d_msg(other.d_msg){};
  virtual const char *what() const noexcept override { return d_msg.c_str(); };
  virtual ~MolSanitizeException() noexcept {};
  virtual MolSanitizeException *copy() const {
    return new MolSanitizeException(*this);
  };
  virtual std::string getType() const { return "MolSanitizeException"; };

 protected:
  std::string d_msg;
};

class RDKIT_GRAPHMOL_EXPORT AtomSanitizeException
    : public MolSanitizeException {
 public:
  AtomSanitizeException(const char *msg, unsigned int atomIdx)
      : MolSanitizeException(msg), d_atomIdx(atomIdx){};
  AtomSanitizeException(const std::string &msg, unsigned int atomIdx)
      : MolSanitizeException(msg), d_atomIdx(atomIdx){};
  AtomSanitizeException(const AtomSanitizeException &other)
      : MolSanitizeException(other), d_atomIdx(other.d_atomIdx){};
  unsigned int getAtomIdx() const { return d_atomIdx; };
  virtual ~AtomSanitizeException() noexcept {};
  virtual MolSanitizeException *copy() const {
    return new AtomSanitizeException(*this);
  };
  virtual std::string getType() const { return "AtomSanitizeException"; };

 protected:
  unsigned int d_atomIdx;
};

class RDKIT_GRAPHMOL_EXPORT AtomValenceException
    : public AtomSanitizeException {
 public:
  AtomValenceException(const char *msg, unsigned int atomIdx)
      : AtomSanitizeException(msg, atomIdx){};
  AtomValenceException(const std::string &msg, unsigned int atomIdx)
      : AtomSanitizeException(msg, atomIdx){};
  AtomValenceException(const AtomValenceException &other)
      : AtomSanitizeException(other){};
  virtual ~AtomValenceException() noexcept {};
  MolSanitizeException *copy() const {
    return new AtomValenceException(*this);
  };
  std::string getType() const { return "AtomValenceException"; };
};

class RDKIT_GRAPHMOL_EXPORT AtomKekulizeException
    : public AtomSanitizeException {
 public:
  AtomKekulizeException(const char *msg, unsigned int atomIdx)
      : AtomSanitizeException(msg, atomIdx){};
  AtomKekulizeException(const std::string &msg, unsigned int atomIdx)
      : AtomSanitizeException(msg, atomIdx){};
  AtomKekulizeException(const AtomKekulizeException &other)
      : AtomSanitizeException(other){};
  virtual ~AtomKekulizeException() noexcept {};
  MolSanitizeException *copy() const {
    return new AtomKekulizeException(*this);
  };
  std::string getType() const { return "AtomKekulizeException"; };
};

class RDKIT_GRAPHMOL_EXPORT KekulizeException : public MolSanitizeException {
 public:
  KekulizeException(const char *msg, const std::vector<unsigned int> &indices)
      : MolSanitizeException(msg), d_atomIndices(indices){};
  KekulizeException(const std::string &msg,
                    const std::vector<unsigned int> &indices)
      : MolSanitizeException(msg), d_atomIndices(indices){};
  KekulizeException(const KekulizeException &other)
      : MolSanitizeException(other), d_atomIndices(other.d_atomIndices){};
  const std::vector<unsigned int> &getAtomIndices() const {
    return d_atomIndices;
  };
  virtual ~KekulizeException() noexcept {};
  MolSanitizeException *copy() const { return new KekulizeException(*this); };
  std::string getType() const { return "KekulizeException"; };

 protected:
  std::vector<unsigned int> d_atomIndices;
};

}  // namespace RDKit

#endif
