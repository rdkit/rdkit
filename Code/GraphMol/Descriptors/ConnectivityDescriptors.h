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

    //! From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991) 
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi0v(const ROMol &mol,bool force=false);
    const std::string chi0vVersion="1.1.0";
    //! From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991) 
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi1v(const ROMol &mol,bool force=false);
    const std::string chi1vVersion="1.1.0";
    //! From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991) 
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi2v(const ROMol &mol,bool force=false);    
    const std::string chi2vVersion="1.1.0";
    //! From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991) 
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi3v(const ROMol &mol,bool force=false);    
    const std::string chi3vVersion="1.1.0";
    //! From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991) 
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi4v(const ROMol &mol,bool force=false);    
    const std::string chi4vVersion="1.1.0";
    //! From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991) 
    /*!
      \param mol           the molecule of interest
      \param n             the order of the connectivity index
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChiNv(const ROMol &mol,unsigned int n,bool force=false);    
    const std::string chiNvVersion="1.1.0";

    //! Similar to Hall Kier ChiXv, but uses nVal instead of valence
    //!   This makes a big difference after we get out of the first row.
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi0n(const ROMol &mol,bool force=false);
    const std::string chi0nVersion="1.1.0";
    //! Similar to Hall Kier ChiXv, but uses nVal instead of valence
    //!   This makes a big difference after we get out of the first row.
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi1n(const ROMol &mol,bool force=false);
    const std::string chi1nVersion="1.1.0";
    //! Similar to Hall Kier ChiXv, but uses nVal instead of valence
    //!   This makes a big difference after we get out of the first row.
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi2n(const ROMol &mol,bool force=false);    
    const std::string chi2nVersion="1.1.0";
    //! Similar to Hall Kier ChiXv, but uses nVal instead of valence
    //!   This makes a big difference after we get out of the first row.
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi3n(const ROMol &mol,bool force=false);    
    const std::string chi3nVersion="1.1.0";
    //! Similar to Hall Kier ChiXv, but uses nVal instead of valence
    //!   This makes a big difference after we get out of the first row.
    /*!
      \param mol           the molecule of interest
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChi4n(const ROMol &mol,bool force=false);    
    const std::string chi4nVersion="1.1.0";
    //! Similar to Hall Kier ChiXv, but uses nVal instead of valence
    //!   This makes a big difference after we get out of the first row.
    /*!
      \param mol           the molecule of interest
      \param n             the order of the connectivity index
      \param force         forces the value to be recalculated instead
                           of pulled from the cache
    */
    double calcChiNn(const ROMol &mol,unsigned int n,bool force=false);    
    const std::string chiNnVersion="1.1.0";

    //! calculate the Hall-Kier alpha value for a molecule
    //! From equation (58) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    /*!
      \param mol           the molecule of interest
      \param atomContribs  if provided, this will be used to return the contributions
                           of the individual atoms to the value. These do not
                           neccessarily sum to the full value.
                           Note: this can be a time-consuming calculation.
    */
    double calcHallKierAlpha(const ROMol &mol,std::vector<double> *atomContribs=0);
    const std::string hallKierAlphaVersion="1.2.0";

    //! calculate the Hall-Kier kappa1 value for a molecule
    //! From equations (58) and (59) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    /*!
      \param mol           the molecule of interest
    */
    double calcKappa1(const ROMol &mol);
    const std::string kappa1Version="1.1.0";

    //! calculate the Hall-Kier kappa2 value for a molecule
    //! From equations (58) and (60) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    /*!
      \param mol           the molecule of interest
    */
    double calcKappa2(const ROMol &mol);
    const std::string kappa2Version="1.1.0";

    //! calculate the Hall-Kier kappa3 value for a molecule
    //! From equations (58), (61) and (62) of Rev. Comp. Chem. vol 2, 367-422, (1991)
    /*!
      \param mol           the molecule of interest
    */
    double calcKappa3(const ROMol &mol);
    const std::string kappa3Version="1.1.0";

  } // end of namespace Descriptors
}

#endif
