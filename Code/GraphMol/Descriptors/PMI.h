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
    const std::string NPR1Version = "1.0.0";
    //! Normalized principal moments ratio 1 (=I1/I3)
    //!  from Sauer and Schwarz JCIM 43:987-1003 (2003)
    //!  https://dx.doi.org/10.1021/ci025599w
    double NPR1(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string NPR2Version = "1.0.0";
    //! Normalized principal moments ratio 2 (=I2/I3)
    //!  from Sauer and Schwarz JCIM 43:987-1003 (2003)
    //!  https://dx.doi.org/10.1021/ci025599w
    double NPR2(const ROMol&, int confId=-1, bool useAtomicMasses=true);

    const std::string PMI1Version = "1.0.0";
    //! First (smallest) principal moment of inertia
    double PMI1(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string PMI2Version = "1.0.0";
    //! second principal moment of inertia
    double PMI2(const ROMol&, int confId=-1, bool useAtomicMasses=true);
    const std::string PMI3Version = "1.0.0";
    //! Third (largest) principal moment of inertia
    double PMI3(const ROMol&, int confId=-1, bool useAtomicMasses=true);

  }
}
#endif
#endif
