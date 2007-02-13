# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension
import RDConfig


#our generic catalog stuff
libDirs.extend([rdCatalogDir,rdLibDir,rdBitDir,rdGeometryDir])
libraries.extend(["MolCatalog",
		  "GraphMol",
		  "Catalogs",
		  "DataStructs",
                  "RDGeometry",
		  ])

#compileArgs.append(vfInc)
#libDirs.append(vfLibDir)
#libraries.extend(vfLibs)


setup(name="Chem.rdMolCatalog",version="1.0",
  ext_modules=[Extension("Chem.rdMolCatalog",["rdMolCatalog.cpp",
					      ],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs)])

