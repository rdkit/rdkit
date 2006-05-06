# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
  os.path.join(RDConfig.RDBaseDir,'Code','DataStructs'),
  ]

libraries += ["DataStructs","RDBoost",'RDGeneral']

setup(name="DataStructs.cDataStructs",version="2.0",
  ext_modules=[Extension("DataStructs.cDataStructs",
                         ["DataStructs.cpp",
                          "DiscreteValueVect.cpp",
                          "wrap_SparseBV.cpp",
                          "wrap_ExplicitBV.cpp",
                          'wrap_BitOps.cpp',
                          'wrap_Utils.cpp',
                          ],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_link_args=linkArgs,
                         extra_compile_args=compileArgs)])

