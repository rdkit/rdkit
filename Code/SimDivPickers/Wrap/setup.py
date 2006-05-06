# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension
import RDConfig


libDirs.extend([
         os.path.join(RDConfig.RDBaseDir,'Code','SimDivPickers'),
         os.path.join(RDConfig.RDBaseDir, 'Code','ML','Cluster','Murtagh'),
         ])

libraries.extend(["SimDivPickers", "hc"])

setup(name="SimDivFilters.rdSimDivPickers", version="2.0",
      ext_modules=[Extension("SimDivFilters.rdSimDivPickers",
                             ["MaxMinPicker.cpp",
                              "HierarchicalClusterPicker.cpp",
                              "rdSimDivPickers.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_compile_args=compileArgs)])
