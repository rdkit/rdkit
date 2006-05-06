# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension

setup(name="ML.InfoTheory.cEntropy",version="1.0",
  ext_modules=[Extension("ML.InfoTheory.cEntropy",["cEntropy.cpp"])])

