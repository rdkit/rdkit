# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension
import os.path
import RDConfig

libDirs.extend([rdLibDir,
		os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','PartialCharges'),
                rdGeometryDir,
		]
	       )

libraries.extend(["PartialCharges", "GraphMol", "RDGeometry"])

setup(name="Chem.rdPartialCharges", version="2.0",
      ext_modules=[Extension("Chem.rdPartialCharges",
                             ["rdPartialCharges.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_compile_args=compileArgs)])
