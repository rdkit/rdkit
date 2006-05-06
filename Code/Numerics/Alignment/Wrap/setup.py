# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'Alignment'),
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         ]

libraries = ["Alignment","RDGeometry"] + libraries

setup(name="Numerics.rdAlignment", version="2.0",
      ext_modules=[Extension("Numerics.rdAlignment",
                             ["rdAlignment.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
