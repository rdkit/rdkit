# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension
import RDConfig


cmimDir=os.path.join(RDConfig.RDBaseDir,'External','cmim-1.0')
libDirs.extend([
         cmimDir,
         rdBitDir,
         ])
libraries.extend(["DataStructs","fastentropy"])

compileArgs.append("-I%s"%(cmimDir))


setup(name="ML.FeatureSelect.rdFeatSelect", version="1.0",
      ext_modules=[Extension("ML.FeatureSelect.rdFeatSelect",
                             ["rdFeatSelect.cpp",
                              ],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_compile_args=compileArgs)])
