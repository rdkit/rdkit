//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef PMI_H_SEPT2016
#define PMI_H_SEPT2016

#ifdef RDK_BUILD_DESCRIPTORS3D
namespace RDKit {
  class ROMol;
  namespace Descriptors {
    //! Normalized principal moments ratio 1 (=I1/I3)
    //!  from Sauer and Schwarz JCIM 43:987-1003 (2003)
    //!  https://dx.doi.org/10.1021/ci025599w
    double NPR1(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string NPR1Version = "1.0.0";
    //! Normalized principal moments ratio 2 (=I2/I3)
    //!  from Sauer and Schwarz JCIM 43:987-1003 (2003)
    //!  https://dx.doi.org/10.1021/ci025599w
    double NPR2(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string NPR2Version = "1.0.0";

    //! First (smallest) principal moment of inertia
    double PMI1(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string PMI1Version = "1.0.0";
    //! second principal moment of inertia
    double PMI2(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string PMI2Version = "1.0.0";
    //! Third (largest) principal moment of inertia
    double PMI3(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string PMI3Version = "1.0.0";

    /*!
     Radius of gyration
       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

     Definition:
        for planar molecules: sqrt( sqrt(pm3*pm2)/MW )
        for nonplanar molecules: sqrt( 2*pi*pow(pm3*pm2*pm1,1/3)/MW )
    */
    double radiusOfGyration(const ROMol&,int confId=-1,
      bool useAtomicMasses=true);
    const std::string radiusOfGyrationVersion = "1.0.0";
    /*!
     Inertial shape factor
       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

     Definition:
       pm2 / (pm1*pm3)
    */
    double inertialShapeFactor(const ROMol&,int confId=-1,
      bool useAtomicMasses=true);
    const std::string inertialShapeFactorVersion = "1.0.0";
    /*!
     Molecular eccentricity
       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

     Definition:
     sqrt(pm3**2 -pm1**2) / pm3**2
    */
    double eccentricity(const ROMol&,int confId=-1,
      bool useAtomicMasses=true);
    const std::string eccentricityVersion = "1.0.0";
  /*!
     molecular asphericity
       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

     Definition:
     0.5 * ((pm3-pm2)**2 + (pm3-pm1)**2 + (pm2-pm1)**2)/(pm1**2+pm2**2+pm3**2)
    */
    double asphericity(const ROMol&,int confId=-1,
      bool useAtomicMasses=true);
    const std::string asphericityVersion = "1.0.0";
    /*!
     Spherocity index
       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

     Definition:
     3 * pm1 / (pm1+pm2+pm3) where the moments are calculated without weights
    */
    double spherocityIndex(const ROMol&,int confId=-1);
    const std::string spherocityIndexVersion = "1.0.0";
  }
}
#endif
#endif
