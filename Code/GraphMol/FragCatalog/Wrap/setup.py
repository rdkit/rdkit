# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension
import RDConfig


#our generic catalog stuff
libDirs.extend([rdCatalogDir,rdLibDir,rdBitDir,rdGeometryDir])
libraries.extend(["FragCat",
		  "SmilesParse",
		  "Substruct",
		  "Subgraphs",
		  "FileParsers",
		  "GraphMol",
		  "Catalogs",
		  "DataStructs",
                  "RDGeometry",
		  ])

libDirs.append(lapackpp)
libraries.extend(lapackLibs)
compileArgs.append(lapackInc)

compileArgs.append(vfInc)
libDirs.append(vfLibDir)
libraries.extend(vfLibs)


setup(name="Chem.rdfragcatalog",version="1.0",
  ext_modules=[Extension("Chem.rdfragcatalog",["rdfragcatalog.cpp",
                                               "FragCatalog.cpp",
                                               "FragCatParams.cpp",
                                               "FragCatGenerator.cpp",
                                               "FragFPGenerator.cpp",
                                               ],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs)])

