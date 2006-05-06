# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'libs'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'DistGeomHelpers'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'ForceFieldHelpers'),
         os.path.join(RDConfig.RDBaseDir,'Code','ForceField'),
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'EigenSolvers'),
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'Optimizer'),
         os.path.join(RDConfig.RDBaseDir,'Code','DistGeom'),
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         ]

libraries += ["DistGeomHelpers", "GraphMol","ForceFieldHelpers", "DistGeom", "ForceFields",
              "Optimizer", "EigenSolvers", "RDGeometry", "RDGeneral"]

setup(name="Chem.rdDistGeom", version="2.0",
      ext_modules=[Extension("Chem.rdDistGeom",
                             ["rdDistGeom.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
