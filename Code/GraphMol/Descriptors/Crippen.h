//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

/*! \file Crippen.h

  \brief Use MolDescriptors.h in client code.

*/
#ifndef __RD_CRIPPEN_H__
#define __RD_CRIPPEN_H__

#include <string>
#include <vector>
#include <boost/smart_ptr.hpp>

namespace RDKit {
  class ROMol;
  namespace Descriptors {
    const std::string crippenVersion="1.0.0";

    //! a class used to store Crippen parameters
    class CrippenParams {
    public:
      boost::shared_ptr<const ROMol> dp_pattern;
      std::string label;
      std::string smarts;
      double logp;
      double mr;
      ~CrippenParams();
    };

    //! singleton class for retrieving Crippen parameters
    /*!
      Use the singleton like this:
      
      \verbatim
      CrippenParamCollection *params=CrippenParamCollection::getParams();
      \endverbatim

      If you have your own parameter data, it can be supplied as a string:
      \verbatim
      CrippenParamCollection *params=CrippenParamCollection::getParams(myParamData);
      \endverbatim
      You are responsible for making sure that the data is in the correct
      format (see Crippen.cpp for an example).
      
    */
    class CrippenParamCollection {
    public:
      typedef std::vector<CrippenParams> ParamsVect;
      static CrippenParamCollection *getParams(const std::string &paramData="");
      ParamsVect::const_iterator begin() const { return d_params.begin(); };
      ParamsVect::const_iterator end() const { return d_params.end(); };
      
    private:
      //! to force this to be a singleton, the constructor must be private
      CrippenParamCollection(const std::string &paramData);
      static class CrippenParamCollection *ds_instance;    //!< the singleton
      ParamsVect d_params;                                 //!< the parameters
    };
  } // end of namespace Descriptors
}

#endif
