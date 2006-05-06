# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'libs'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'ShapeHelpers'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'MolTransforms'),
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'EigenSolvers'),
         os.path.join(RDConfig.RDBaseDir,'Code',"DataStructs"),
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         ]

libraries += ["ShapeHelpers", "MolTransforms", "GraphMol", "EigenSolvers", 
              "RDGeometry", "DataStructs", "RDGeneral"]

setup(name="Chem.rdShapeHelpers", version="2.0",
      ext_modules=[Extension("Chem.rdShapeHelpers",
                             ["rdShapeHelpers.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
