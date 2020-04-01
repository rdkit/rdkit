//
// Copyright 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FILEPARSEEXCEPTION_H
#define RD_FILEPARSEEXCEPTION_H

#include <string>
#include <stdexcept>

namespace RDKit {
//! used by various file parsing classes to indicate a parse error
class FileParseException : public std::runtime_error {
 public:
  //! construct with an error message
  explicit FileParseException(const char *msg)
      : std::runtime_error("FileParseException"), _msg(msg){};
  //! construct with an error message
  explicit FileParseException(const std::string msg)
      : std::runtime_error("FileParseException"), _msg(msg){};
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); };
  ~FileParseException() noexcept {};

 private:
  std::string _msg;
};
}  // namespace RDKit

#endif
