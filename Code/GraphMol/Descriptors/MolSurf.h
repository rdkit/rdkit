//
//  Copyright (C) 2007-2011 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file MolSurf.h

  \brief Use MolDescriptors.h in client code.

*/
#ifndef __RD_MOLSURF_H__
#define __RD_MOLSURF_H__

#include<vector>

namespace RDKit{
  class ROMol;
  namespace Descriptors {
    const std::string labuteASAVersion="1.0.2";

    //! calculates atomic contributions to Labute's Approximate Surface Area
    /*!  
      Definition from P. Labute's article in the Journal of the
      Chemical Computing Group and J. Mol. Graph. Mod.  _18_ 464-477
      (2000)

      \param mol        the molecule of interest
      \param Vi         used to return the explict atom contribs
      \param hContrib   used to return the H contributions (if calculated)
      \param includeHs  (optional) if this is true (the default),
          the contribution of H atoms to the ASA will be included.
      \param force      (optional) calculate the values even if they are cached.

      \return the sum of the atomic contributions

    */
    double getLabuteAtomContribs(const ROMol &mol,
				 std::vector<double> &Vi,
				 double &hContrib,
				 bool includeHs=true,
				 bool force=false);

    //! calculates Labute's Approximate Surface Area (ASA from MOE)
    /*!  
      Definition from P. Labute's article in the Journal of the
      Chemical Computing Group and J. Mol. Graph. Mod.  _18_ 464-477
      (2000)

      \param mol        the molecule of interest
      \param includeHs  (optional) if this is true (the default),
          the contribution of H atoms to the ASA will be included.
      \param force      (optional) calculate the value even if it's cached.
	  
    */
    double calcLabuteASA(const ROMol &mol,bool includeHs=true,
			 bool force=false);

    const std::string tpsaVersion="1.1.0";
    //! calculates atomic contributions to the TPSA value
    /*!  
      The TPSA definition is from:
      P. Ertl, B. Rohde, P. Selzer
       Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
       Contributions and Its Application to the Prediction of Drug Transport 
       Properties, J.Med.Chem. 43, 3714-3717, 2000

      \param mol        the molecule of interest
      \param Vi         used to return the atom contribs
      \param force      (optional) calculate the values even if they are cached.

      \return the sum of the atomic contributions

    */
    double getTPSAAtomContribs(const ROMol &mol,
                               std::vector<double> &Vi,
                               bool force=false);

    //! calculates the TPSA value for a molecule
    /*!  
      The TPSA definition is from:
      P. Ertl, B. Rohde, P. Selzer
       Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
       Contributions and Its Application to the Prediction of Drug Transport 
       Properties, J.Med.Chem. 43, 3714-3717, 2000

      \param mol        the molecule of interest
      \param force      (optional) calculate the value even if it's cached.
	  
    */
    double calcTPSA(const ROMol &mol,
                    bool force=false);

    std::vector<double> calcSlogP_VSA(const ROMol &mol,std::vector<double> *bins=0,
                                      bool force=false);
    std::vector<double> calcSMR_VSA(const ROMol &mol,std::vector<double> *bins=0,
                                      bool force=false);
    std::vector<double> calcPEOE_VSA(const ROMol &mol,std::vector<double> *bins=0,
                                      bool force=false);
    
  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
