//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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
    const std::string crippenVersion="1.2.0";

    //! generate atomic contributions to the Wildman-Crippen LogP and MR
    //! estimates for a molecule
    /*!
      Uses an atom-based scheme based on the values in the paper:
        S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)

      \param mol           the molecule of interest
      \param logpContribs  used to return the logp contributions, must
                           be equal in length to the number of atoms
      \param mrContribs    used to return the MR contributions, must
                           be equal in length to the number of atoms
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
      \param atomTypes     if provided will be used to return the indices
                           of the atom types, should be as long as the
                           number of atoms
      \param atomTypeLabels   if provided will be used to return the labels
                              of the atom types, should be as long as the
                              number of atoms
	  
    */
    void getCrippenAtomContribs(const ROMol &mol,
				std::vector<double> &logpContribs,
				std::vector<double> &mrContribs,
				bool force=false,
                                std::vector<unsigned int> *atomTypes=0,
                                std::vector<std::string> *atomTypeLabels=0
                                );
    

    //! generate Wildman-Crippen LogP and MR estimates for a molecule
    /*!
      Uses an atom-based scheme based on the values in the paper:
        S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)

      \param mol        the molecule of interest
      \param logp       used to return the logp estimate
      \param mr         used to return the MR estimate
      \param includeHs  (optional) if this is true (the default), a
          copy of \c mol is made and Hs are added to it.  If false,
	  Hs that are not explicitly present in the graph will not
	  be included.
      \param force      forces the value to be recalculated instead of
                        pulled from the cache
	  
    */
    void calcCrippenDescriptors(const ROMol &mol,double &logp,double &mr,
				bool includeHs=true,bool force=false);
    
    //! a class used to store Crippen parameters
    class CrippenParams {
    public:
      boost::shared_ptr<const ROMol> dp_pattern;
      unsigned int idx;
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
      static const CrippenParamCollection *getParams(const std::string &paramData="");
      ParamsVect::const_iterator begin() const { return d_params.begin(); };
      ParamsVect::const_iterator end() const { return d_params.end(); };
      
      CrippenParamCollection(const std::string &paramData);
    private:
      ParamsVect d_params;                                 //!< the parameters
    };
  } // end of namespace Descriptors
}

#endif
