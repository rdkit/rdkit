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
#include <string>

// RDKIT_JAVAWRAPPERS_EXPORT does not get defined in RDGeneral/export.h,
// and we only use it here for non-windows builds, so just define it based
// on RDKIT_RDGENERAL_EXPORT
#if defined(RDKIT_DYN_LINK) && defined(WIN32) && defined(BOOST_HAS_DECLSPEC)
#define RDKIT_JAVAWRAPPERS_EXPORT
#else
#define RDKIT_JAVAWRAPPERS_EXPORT RDKIT_RDGENERAL_EXPORT
#endif

namespace RDKit {

class RDKIT_JAVAWRAPPERS_EXPORT GenericRDKitException : public std::exception {
 public:
  GenericRDKitException(const std::string &i) : _value(i) {}
  GenericRDKitException(const char *msg) : _value(msg) {}
  const char *what() const noexcept override { return _value.c_str(); }
  ~GenericRDKitException() noexcept = default;

 private:
  std::string _value;
};
}  // namespace RDKit
