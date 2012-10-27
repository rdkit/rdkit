//
//  Copyright (C) 2012 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file ConnectivityDescriptors.h

  \brief Use MolDescriptors.h in client code.

*/
#ifndef __RD_CONNECTIVITYDESCRIPTORS_H__
#define __RD_CONNECTIVITYDESCRIPTORS_H__

#include <string>
#include <vector>
#include <boost/smart_ptr.hpp>

namespace RDKit {
  class ROMol;
  namespace Descriptors {

    //! 
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
	  
    */
    double calcChi0v(const ROMol &mol,bool force=false);
    double calcChi1v(const ROMol &mol,bool force=false);
    double calcChi2v(const ROMol &mol,bool force=false);    
    double calcChi3v(const ROMol &mol,bool force=false);    
    double calcChi4v(const ROMol &mol,bool force=false);    

    double calcChi0n(const ROMol &mol,bool force=false);
    double calcChi1n(const ROMol &mol,bool force=false);
    double calcChi2n(const ROMol &mol,bool force=false);    
    double calcChi3n(const ROMol &mol,bool force=false);    
    double calcChi4n(const ROMol &mol,bool force=false);    

    double calcHallKierAlpha(const ROMol &mol);

  } // end of namespace Descriptors
}

#endif
