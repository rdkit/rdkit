#
# Copyright (C) 2016 greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Descriptors derived from a molecule's 3D structure

"""
from __future__ import print_function
from rdkit.Chem import rdMolDescriptors

<<<<<<< HEAD
<<<<<<< HEAD
if hasattr(rdMolDescriptors,'CalcPMI1'):
    PMI1 = lambda *x,**y:rdMolDescriptors.CalcPMI1(*x,**y)
    PMI1.version = rdMolDescriptors._CalcPMI1_version
    PMI1.__doc__= """ First (smallest) principal moment of inertia
=======
=======
>>>>>>> rdkit/master
if hasattr(rdMolDescriptors, 'CalcPMI1'):
  PMI1 = lambda *x, **y: rdMolDescriptors.CalcPMI1(*x, **y)
  PMI1.version = rdMolDescriptors._CalcPMI1_version
  PMI1.__doc__ = """ First (smallest) principal moment of inertia
<<<<<<< HEAD
>>>>>>> rdkit/master
=======
>>>>>>> rdkit/master

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    PMI2 = lambda *x,**y:rdMolDescriptors.CalcPMI2(*x,**y)
    PMI2.version = rdMolDescriptors._CalcPMI2_version
    PMI2.__doc__= """ Second principal moment of inertia
=======
  PMI2 = lambda *x, **y: rdMolDescriptors.CalcPMI2(*x, **y)
  PMI2.version = rdMolDescriptors._CalcPMI2_version
  PMI2.__doc__ = """ Second principal moment of inertia
>>>>>>> rdkit/master
=======
  PMI2 = lambda *x, **y: rdMolDescriptors.CalcPMI2(*x, **y)
  PMI2.version = rdMolDescriptors._CalcPMI2_version
  PMI2.__doc__ = """ Second principal moment of inertia
>>>>>>> rdkit/master

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    PMI3 = lambda *x,**y:rdMolDescriptors.CalcPMI3(*x,**y)
    PMI3.version = rdMolDescriptors._CalcPMI3_version
    PMI3.__doc__= """ Third (largest) principal moment of inertia
=======
  PMI3 = lambda *x, **y: rdMolDescriptors.CalcPMI3(*x, **y)
  PMI3.version = rdMolDescriptors._CalcPMI3_version
  PMI3.__doc__ = """ Third (largest) principal moment of inertia
>>>>>>> rdkit/master
=======
  PMI3 = lambda *x, **y: rdMolDescriptors.CalcPMI3(*x, **y)
  PMI3.version = rdMolDescriptors._CalcPMI3_version
  PMI3.__doc__ = """ Third (largest) principal moment of inertia
>>>>>>> rdkit/master

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    NPR1 = lambda *x,**y:rdMolDescriptors.CalcNPR1(*x,**y)
    NPR1.version = rdMolDescriptors._CalcNPR1_version
    NPR1.__doc__= """ Normalized principal moments ratio 1 (=I1/I3)
=======
  NPR1 = lambda *x, **y: rdMolDescriptors.CalcNPR1(*x, **y)
  NPR1.version = rdMolDescriptors._CalcNPR1_version
  NPR1.__doc__ = """ Normalized principal moments ratio 1 (=I1/I3)
>>>>>>> rdkit/master
=======
  NPR1 = lambda *x, **y: rdMolDescriptors.CalcNPR1(*x, **y)
  NPR1.version = rdMolDescriptors._CalcNPR1_version
  NPR1.__doc__ = """ Normalized principal moments ratio 1 (=I1/I3)
>>>>>>> rdkit/master

        from Sauer and Schwarz JCIM 43:987-1003 (2003)
        https://dx.doi.org/10.1021/ci025599w


    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    NPR2 = lambda *x,**y:rdMolDescriptors.CalcNPR2(*x,**y)
    NPR2.version = rdMolDescriptors._CalcNPR2_version
    NPR2.__doc__= """ Normalized principal moments ratio 2 (=I2/I3)
=======
  NPR2 = lambda *x, **y: rdMolDescriptors.CalcNPR2(*x, **y)
  NPR2.version = rdMolDescriptors._CalcNPR2_version
  NPR2.__doc__ = """ Normalized principal moments ratio 2 (=I2/I3)
>>>>>>> rdkit/master
=======
  NPR2 = lambda *x, **y: rdMolDescriptors.CalcNPR2(*x, **y)
  NPR2.version = rdMolDescriptors._CalcNPR2_version
  NPR2.__doc__ = """ Normalized principal moments ratio 2 (=I2/I3)
>>>>>>> rdkit/master

        from Sauer and Schwarz JCIM 43:987-1003 (2003)
        https://dx.doi.org/10.1021/ci025599w


    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    RadiusOfGyration = lambda *x,**y:rdMolDescriptors.CalcRadiusOfGyration(*x,**y)
    RadiusOfGyration.version = rdMolDescriptors._CalcRadiusOfGyration_version
    RadiusOfGyration.__doc__= """ Radius of gyration
=======
  RadiusOfGyration = lambda *x, **y: rdMolDescriptors.CalcRadiusOfGyration(*x, **y)
  RadiusOfGyration.version = rdMolDescriptors._CalcRadiusOfGyration_version
  RadiusOfGyration.__doc__ = """ Radius of gyration
>>>>>>> rdkit/master
=======
  RadiusOfGyration = lambda *x, **y: rdMolDescriptors.CalcRadiusOfGyration(*x, **y)
  RadiusOfGyration.version = rdMolDescriptors._CalcRadiusOfGyration_version
  RadiusOfGyration.__doc__ = """ Radius of gyration
>>>>>>> rdkit/master

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

       Definition:
         for planar molecules: sqrt( sqrt(pm3*pm2)/MW )
         for nonplanar molecules: sqrt( 2*pi*pow(pm3*pm2*pm1,1/3)/MW )

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    InertialShapeFactor = lambda *x,**y:rdMolDescriptors.CalcInertialShapeFactor(*x,**y)
    InertialShapeFactor.version = rdMolDescriptors._CalcInertialShapeFactor_version
    InertialShapeFactor.__doc__= """ Inertial shape factor
=======
  InertialShapeFactor = lambda *x, **y: rdMolDescriptors.CalcInertialShapeFactor(*x, **y)
  InertialShapeFactor.version = rdMolDescriptors._CalcInertialShapeFactor_version
  InertialShapeFactor.__doc__ = """ Inertial shape factor
>>>>>>> rdkit/master
=======
  InertialShapeFactor = lambda *x, **y: rdMolDescriptors.CalcInertialShapeFactor(*x, **y)
  InertialShapeFactor.version = rdMolDescriptors._CalcInertialShapeFactor_version
  InertialShapeFactor.__doc__ = """ Inertial shape factor
>>>>>>> rdkit/master

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

       Definition:
         pm2 / (pm1*pm3)

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    Eccentricity = lambda *x,**y:rdMolDescriptors.CalcEccentricity(*x,**y)
    Eccentricity.version = rdMolDescriptors._CalcEccentricity_version
    Eccentricity.__doc__= """ molecular eccentricity
=======
  Eccentricity = lambda *x, **y: rdMolDescriptors.CalcEccentricity(*x, **y)
  Eccentricity.version = rdMolDescriptors._CalcEccentricity_version
  Eccentricity.__doc__ = """ molecular eccentricity
>>>>>>> rdkit/master
=======
  Eccentricity = lambda *x, **y: rdMolDescriptors.CalcEccentricity(*x, **y)
  Eccentricity.version = rdMolDescriptors._CalcEccentricity_version
  Eccentricity.__doc__ = """ molecular eccentricity
>>>>>>> rdkit/master

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

       Definition:
         sqrt(pm3**2 -pm1**2) / pm3**2

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    Asphericity = lambda *x,**y:rdMolDescriptors.CalcAsphericity(*x,**y)
    Asphericity.version = rdMolDescriptors._CalcAsphericity_version
    Asphericity.__doc__= """ molecular asphericity
=======
  Asphericity = lambda *x, **y: rdMolDescriptors.CalcAsphericity(*x, **y)
  Asphericity.version = rdMolDescriptors._CalcAsphericity_version
  Asphericity.__doc__ = """ molecular asphericity
>>>>>>> rdkit/master
=======
  Asphericity = lambda *x, **y: rdMolDescriptors.CalcAsphericity(*x, **y)
  Asphericity.version = rdMolDescriptors._CalcAsphericity_version
  Asphericity.__doc__ = """ molecular asphericity
>>>>>>> rdkit/master

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

       Definition:
         0.5 * ((pm3-pm2)**2 + (pm3-pm1)**2 + (pm2-pm1)**2)/(pm1**2+pm2**2+pm3**2)

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
<<<<<<< HEAD
<<<<<<< HEAD
    SpherocityIndex = lambda *x,**y:rdMolDescriptors.CalcSpherocityIndex(*x,**y)
    SpherocityIndex.version = rdMolDescriptors._CalcSpherocityIndex_version
    SpherocityIndex.__doc__= """ Molecular spherocityIndex
=======
  SpherocityIndex = lambda *x, **y: rdMolDescriptors.CalcSpherocityIndex(*x, **y)
  SpherocityIndex.version = rdMolDescriptors._CalcSpherocityIndex_version
  SpherocityIndex.__doc__ = """ Molecular spherocityIndex
>>>>>>> rdkit/master
=======
  SpherocityIndex = lambda *x, **y: rdMolDescriptors.CalcSpherocityIndex(*x, **y)
  SpherocityIndex.version = rdMolDescriptors._CalcSpherocityIndex_version
  SpherocityIndex.__doc__ = """ Molecular spherocityIndex
>>>>>>> rdkit/master

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       http://dx.doi.org/10.1002/9783527618279.ch37

       Definition:
         3 * pm1 / (pm1+pm2+pm3) where the moments are calculated without weights

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

    """
