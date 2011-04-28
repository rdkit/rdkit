// $Id: GenericRDKitException.h
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// A generic exception (to cause a corresponding one to be created in Java)

#include <exception>

namespace RDKit {

class GenericRDKitException : public std::exception
{
public:
  GenericRDKitException(const std::string i) : _value(i) {};
  GenericRDKitException(const char *msg) : _value(msg) {};
  std::string message () const { return _value; };
  ~GenericRDKitException () throw () {};
private:
  std::string _value;
};

}


