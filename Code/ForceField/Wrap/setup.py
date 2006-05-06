# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'libs'),
         os.path.join(RDConfig.RDBaseDir,'Code','ForceField'),
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'Optimizer'),
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         ]

libraries += ["ForceFields", "Optimizer", "RDGeometry", "RDGeneral"]

setup(name="ForceField.rdForceField", version="2.0",
      ext_modules=[Extension("ForceField.rdForceField",
                             ["ForceField.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
