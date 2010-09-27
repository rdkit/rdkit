//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_CATALOGENTRY_H__
#define __RD_CATALOGENTRY_H__

#include <iostream>
#include <string>

namespace RDCatalog {
  
  //! Abstract base class to be used to represent an entry in a Catalog
  class CatalogEntry {
  public:
    virtual ~CatalogEntry() = 0;

    //! sets our bit Id
    void setBitId(int bid) {d_bitId = bid;}; 

    //! returns our bit Id
    int getBitId() const {return d_bitId;};

    //! returns a text description of this entry
    virtual std::string getDescription() const = 0;

    //! serializes (pickles) to a stream
    virtual void toStream(std::ostream &ss) const = 0;
    //! returns a string with a serialized (pickled) representation
    virtual std::string Serialize() const = 0;
    //! initializes from a stream pickle
    virtual void initFromStream(std::istream &ss) = 0;
    //! initializes from a string pickle
    virtual void initFromString(const std::string &text) = 0;


    
  private:
    int d_bitId; //!< our bit Id. This needs to be signed so that we can mark uninitialized entries.
  };
}

#endif
    
