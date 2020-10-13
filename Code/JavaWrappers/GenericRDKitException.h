// $Id: GenericRDKitException.h
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// A generic exception (to cause a corresponding one to be created in Java)

#include <RDGeneral/export.h>
#include <exception>

namespace RDKit {

class GenericRDKitException : public std::exception {
 public:
  GenericRDKitException(const std::string &i) : _value(i){};
  GenericRDKitException(const char *msg) : _value(msg){};
  const char *what() const noexcept override { return _value.c_str(); };
  ~GenericRDKitException() noexcept {};

 private:
  std::string _value;
};
}  // namespace RDKit
