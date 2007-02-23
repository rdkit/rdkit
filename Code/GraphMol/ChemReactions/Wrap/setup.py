# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'libs'),
         ]

libraries += ["SmilesParse","FileParsers","ChemReactions","Substruct","GraphMol", "RDGeneral"]

compileArgs.append(vfInc)
libDirs.append(vfLibDir)
libraries.extend(vfLibs)


setup(name="Chem.rdChemReactions", version="1.0",
      ext_modules=[Extension("Chem.rdChemReactions",
                             ["rdChemReactions.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
