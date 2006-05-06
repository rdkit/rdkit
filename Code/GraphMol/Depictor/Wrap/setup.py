# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
  rdLibDir,
  rdBitDir,
  rdGeometryDir,
         ]

libraries += ["Depictor", "GraphMol","DataStructs", "RDGeometry"]

setup(name="Chem.rdDepictor", version="2.0",
      ext_modules=[Extension("Chem.rdDepictor",
                             ["rdDepictor.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
