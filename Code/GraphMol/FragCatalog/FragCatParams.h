//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_FRAG_CAT_PARAMS_H_
#define _RD_FRAG_CAT_PARAMS_H_

#include <Catalogs/CatalogParams.h>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>


namespace RDKit {
  class ROMol;
  typedef std::vector< boost::shared_ptr<ROMol> > MOL_SPTR_VECT;

  //! container for user parameters used to create a fragment catalog
  class FragCatParams : public RDCatalog::CatalogParams {
    // FIX: this container is still missing all the CASE-type functional groups stuff
  public:
    FragCatParams() {
      d_typeStr = "Fragment Catalog Parameters";
      d_lowerFragLen = 0;
      d_upperFragLen = 0;
      d_tolerance = 1e-8;
      d_funcGroups.clear();
    }
    //! construct from a function-group file
    /*!
      \param lLen       the lower limit on fragment size
      \param uLen       the upper limit on fragment size
      \param fgroupFile the name of the function-group file
      \param tol        (optional) the eigenvalue tolerance to be used
                        when comparing fragments
    */
    FragCatParams(unsigned int lLen, unsigned int uLen, std::string fgroupFile, double tol=1e-08);
    //! copy constructor
    FragCatParams(const FragCatParams &other);
    //! construct from a pickle string (serialized representation)
    FragCatParams(const std::string &pickle);

    ~FragCatParams();

    //! returns our lower fragment length
    unsigned int getLowerFragLength() const { return d_lowerFragLen;}
    //! sets our lower fragment length
    void setLowerFragLength(unsigned int lFrLen) {d_lowerFragLen = lFrLen;}
    
    //! returns our upper fragment length
    unsigned int getUpperFragLength() const { return d_upperFragLen;}
    //! sets our upper fragment length
    void setUpperFragLength(unsigned int uFrLen) { d_upperFragLen = uFrLen;}

    //! returns our fragment-comparison tolerance
    double getTolerance() const {return d_tolerance;}
    //! sets our fragment-comparison tolerance
    void setTolerance(double val) { d_tolerance = val;}

    //! returns our number of functional groups
    unsigned int getNumFuncGroups() const {return d_funcGroups.size();}

    //! returns our std::vector of functional groups
    const MOL_SPTR_VECT &getFuncGroups() const;

    //! returns a pointer to a specific functional group
    const ROMol *getFuncGroup(unsigned int fid) const; 

    void toStream(std::ostream &) const;
    std::string Serialize() const;
    void initFromStream(std::istream &ss);
    void initFromString(const std::string &text);
  private:
    unsigned int d_lowerFragLen; 
    unsigned int d_upperFragLen;

    double d_tolerance; //!< tolerance value used when comparing subgraph discriminators

    MOL_SPTR_VECT d_funcGroups;
  };
}

#endif
