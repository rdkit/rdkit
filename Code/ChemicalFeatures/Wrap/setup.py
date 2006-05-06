# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         rdLibDir,
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         os.path.join(RDConfig.RDBaseDir,'Code','ChemicalFeatures'),
         ]

libraries += ["ChemicalFeatures","RDGeometry", "RDGeneral",] 


setup(name="Chem.rdChemicalFeatures", version="2.0",
      ext_modules=[Extension("Chem.rdChemicalFeatures",
                             ["rdChemicalFeatures.cpp",
                              "FreeChemicalFeature.cpp"],                             
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])

