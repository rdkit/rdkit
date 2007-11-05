# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension

libDirs.extend([
         os.path.join(RDConfig.RDBaseDir,'Code','ML','InfoTheory'),
         ])
libraries.extend(["InfoTheory"])

setup(name="ML.Data.cQuantize",version="1.0",
  ext_modules=[Extension("ML.Data.cQuantize",["cQuantize.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
			 extra_compile_args=compileArgs)])

