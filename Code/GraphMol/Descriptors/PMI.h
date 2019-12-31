//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef PMI_H_SEPT2016
#define PMI_H_SEPT2016

#ifdef RDK_BUILD_DESCRIPTORS3D
namespace RDKit {
class ROMol;
namespace Descriptors {
//! Normalized principal moments ratio 1 (=I1/I3)
//!  from Sauer and Schwarz JCIM 43:987-1003 (2003)
//!  https://doi.org/10.1021/ci025599w
RDKIT_DESCRIPTORS_EXPORT double NPR1(const ROMol&, int confId = -1,
                                     bool useAtomicMasses = true,
                                     bool force = false);
const std::string NPR1Version = "1.0.0";
//! Normalized principal moments ratio 2 (=I2/I3)
//!  from Sauer and Schwarz JCIM 43:987-1003 (2003)
//!  https://doi.org/10.1021/ci025599w
RDKIT_DESCRIPTORS_EXPORT double NPR2(const ROMol&, int confId = -1,
                                     bool useAtomicMasses = true,
                                     bool force = false);
const std::string NPR2Version = "1.0.0";

//! First (smallest) principal moment of inertia
RDKIT_DESCRIPTORS_EXPORT double PMI1(const ROMol&, int confId = -1,
                                     bool useAtomicMasses = true,
                                     bool force = false);
const std::string PMI1Version = "1.0.0";
//! second principal moment of inertia
RDKIT_DESCRIPTORS_EXPORT RDKIT_DESCRIPTORS_EXPORT double PMI2(
    const ROMol&, int confId = -1, bool useAtomicMasses = true,
    bool force = false);
const std::string PMI2Version = "1.0.0";
//! Third (largest) principal moment of inertia
RDKIT_DESCRIPTORS_EXPORT double PMI3(const ROMol&, int confId = -1,
                                     bool useAtomicMasses = true,
                                     bool force = false);
const std::string PMI3Version = "1.0.0";

/*!
 Radius of gyration
   from G. A. Arteca "Molecular Shape Descriptors"
   Reviews in Computational Chemistry vol 9
   https://doi.org/10.1002/9780470125861.ch5

 Definition (eq: A4):
    sqrt(t_1 + t_2 + t_3) where t_i is the ith moment from the gyration matrix
*/
RDKIT_DESCRIPTORS_EXPORT double radiusOfGyration(const ROMol&, int confId = -1,
                                                 bool useAtomicMasses = true,
                                                 bool force = false);
const std::string radiusOfGyrationVersion = "1.0.0";
/*!
 Inertial shape factor
   from Todeschini and Consoni "Descriptors from Molecular Geometry"
   Handbook of Chemoinformatics
   https://doi.org/10.1002/9783527618279.ch37

 Definition:
   pm2 / (pm1*pm3)
*/
RDKIT_DESCRIPTORS_EXPORT double inertialShapeFactor(const ROMol&,
                                                    int confId = -1,
                                                    bool useAtomicMasses = true,
                                                    bool force = false);
const std::string inertialShapeFactorVersion = "1.0.0";
/*!
 Molecular eccentricity
   from G. A. Arteca "Molecular Shape Descriptors"
   Reviews in Computational Chemistry vol 9
   https://doi.org/10.1002/9780470125861.ch5

 Definition (eq 4):
 sqrt(pm_3**2 -pm_1**2) / pm_3**2  where pm_i is the ith moment of inertia
*/
RDKIT_DESCRIPTORS_EXPORT double eccentricity(const ROMol&, int confId = -1,
                                             bool useAtomicMasses = true,
                                             bool force = false);
const std::string eccentricityVersion = "1.0.0";
/*!
 molecular asphericity
   from A. Baumgaertner, "Shapes of flexible vesicles"
   J. Chem. Phys. 98:7496 (1993)
   https://doi.org/10.1063/1.464689

   Definition (eq 11):
   0.5 * ((t_3-t_2)**2 + (t_3-t_1)**2 + (t_2-t_1)**2)/(t_1+t_2+t_3)**2
                  where t_i is the ith moment from the gyration matrix

   Some explanation of that definition: the text of the paper mentions axes
   of inertia, but then the definition of the radius of gyration in eq.9 clearly
   uses the moments of the gyration matrix. The text under equation 11 has the
   appropriate inequalities and limits for the moments of gyration, but the
   description of the geometry provided corresponds to the moments of inertia.
   The definition here corresponds to what Dragon generates and seem logical

*/
RDKIT_DESCRIPTORS_EXPORT double asphericity(const ROMol&, int confId = -1,
                                            bool useAtomicMasses = true,
                                            bool force = false);
const std::string asphericityVersion = "1.0.0";
/*!
 Spherocity index
   from Todeschini and Consoni "Descriptors from Molecular Geometry"
   Handbook of Chemoinformatics
   https://doi.org/10.1002/9783527618279.ch37

 Definition:
 3 * t_1 / (t_1+t_2+t_3)
    where the moments of the gyration matrix are calculated without weights
*/
RDKIT_DESCRIPTORS_EXPORT double spherocityIndex(const ROMol&, int confId = -1,
                                                bool force = false);
const std::string spherocityIndexVersion = "1.0.0";
}  // namespace Descriptors
}  // namespace RDKit
#endif
#endif
