# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         rdLibDir,
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','MolChemicalFeatures'),
         vfLibDir,
         ]

libraries += ["MolChemicalFeatures", "SmilesParse", "Substruct", "GraphMol",
              "RDBoost", "RDGeometry", "RDGeneral"] #, vfLibs]

compileArgs.append(vfInc)
libDirs.append(vfLibDir)
libraries.extend(vfLibs)

setup(name="Chem.rdMolChemicalFeatures", version="2.0",
      ext_modules=[Extension("Chem.rdMolChemicalFeatures",
                             ["MolChemicalFeature.cpp",
                              "MolChemicalFeatureFactory.cpp",
                              "ChemicalFeatureUtils.cpp",
                              "rdMolChemicalFeatures.cpp"],
                             
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])

