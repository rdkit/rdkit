//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_CATALOGPARAMS_H__
#define __RD_CATALOGPARAMS_H__

#include <string>

namespace RDCatalog {
  //! abstract base class for the container used to create a catalog
  class CatalogParams {
  public:
    virtual ~CatalogParams() = 0;

    //! returns our type string
    std::string getTypeStr() const { return d_typeStr; };

    //! sets our type string
    void setTypeStr(const std::string &typeStr) { d_typeStr=typeStr; };

    //! serializes (pickles) to a stream
    virtual void toStream(std::ostream &) const = 0;
    //! returns a string with a serialized (pickled) representation
    virtual std::string Serialize() const = 0;
    //! initializes from a stream pickle
    virtual void initFromStream(std::istream &ss) = 0;
    //! initializes from a string pickle
    virtual void initFromString(const std::string &text) = 0;

  protected:
    std::string d_typeStr; //!< our type string
  };
}

#endif
