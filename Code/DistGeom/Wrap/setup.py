# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'Optimizer'),
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'EigenSolvers'),
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         os.path.join(RDConfig.RDBaseDir,'Code','DistGeom'),
         os.path.join(RDConfig.RDBaseDir,'Code','ForceField'),
         ]

libraries = ["DistGeom", "ForceFields", "EigenSolvers", "Optimizer",
              "RDGeometry"] + libraries

setup(name="DistanceGeometry.DistGeom", version="2.0",
      ext_modules=[Extension("DistanceGeometry.DistGeom",
                             ["DistGeom.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
