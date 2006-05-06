# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension


libDirs += [
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'libs'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'MolAlign'),
         os.path.join(RDConfig.RDBaseDir,'Code','Numerics', 'Alignment'),
         os.path.join(RDConfig.RDBaseDir,'Code','Geometry'),
         os.path.join(RDConfig.RDBaseDir,'Code','GraphMol', 'MolTransforms'),
         vfLibDir,
         ]

libraries += ["MolTransforms", "MolAlign", "GraphMol", "Substruct", 
              "Alignment", "RDGeometry", "RDGeneral"]

compileArgs.append(vfInc)
libDirs.append(vfLibDir)
libraries.extend(vfLibs)

setup(name="Chem.rdMolAlign", version="2.0",
      ext_modules=[Extension("Chem.rdMolAlign",
                             ["rdMolAlign.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_link_args=linkArgs,
                             extra_compile_args=compileArgs)])
