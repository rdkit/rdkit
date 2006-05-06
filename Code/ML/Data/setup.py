# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


setup(name="ML.Data.cQuantize",version="1.0",
  ext_modules=[Extension("ML.Data.cQuantize",["cQuantize.cpp"],
			 extra_compile_args=compileArgs)])

