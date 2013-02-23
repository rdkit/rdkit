//
//  Copyright (C) 2013 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file MQN.h

  \brief Use MolDescriptors.h in client code.

*/
#ifndef __RD_MQN_H__
#define __RD_MQN_H__

#include<vector>

namespace RDKit{
  class ROMol;
  namespace Descriptors {
    const std::string MQNVersion="1.0.0";

    //! calculates MQN descriptors
    /*!  
      Definition from 
      Nguyen, K. T., Blum, L. C., Van Deursen, R. & Reymond, J.-L. "Classification of Organic Molecules by Molecular Quantum Numbers." 
      ChemMedChem 4, 1803–1805 (2009).


      \param mol        the molecule of interest
      \param force      (optional) calculate the values even if they are cached.

      \return a vector with the MQNs

    */
    std::vector<unsigned int>  calcMQNs(const ROMol &mol,
                                       bool force=false);
    
  } // end of namespace Descriptors
} //end of namespace RDKit

#endif
