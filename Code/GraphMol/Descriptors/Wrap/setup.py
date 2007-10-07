# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'libs'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'Descriptors'),
         ]

libraries += ["Descriptors","GraphMol","RDGeneral"]

setup(name="Chem.rdMolDescriptors", version="1.0",
      ext_modules=[Extension("Chem.rdMolDescriptors",
                             ["rdMolDescriptors.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
