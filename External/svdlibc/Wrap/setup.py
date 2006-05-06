# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension
libDirs.append(os.path.join(RDConfig.RDBaseDir,'External','svdlibc'))
libraries+=["svd"]
compileArgs.append('-I'+extDir)
setup(name="PySVD.cSVD", version="2.0",
      ext_modules=[Extension("PySVD.cSVD",
                             ["PySVD.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
