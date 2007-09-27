# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension

libDirs += [rdLibDir,rdBitDir,rdCatalogDir,rdGeometryDir]
libraries.extend(["SmilesParse","ChemTransforms","Substruct","FileParsers", "Subgraphs",
           "FragCat", "GraphMol", "Catalogs", "Fingerprints","DataStructs",
           "RDBoost","RDGeneral", "RDGeometry",
           ])

compileArgs.append(lapackInc)
libDirs.append(lapackpp)
libraries.extend(lapackLibs)


compileArgs.append(vfInc)
libDirs.append(vfLibDir)
libraries.extend(vfLibs)

setup(name="Chem.rdchem",version="2.0",
  ext_modules=[Extension("Chem.rdchem",["rdchem.cpp",
                                        'Table.cpp','Atom.cpp',
                                        'Bond.cpp','Mol.cpp',
                                        'Conformer.cpp',
                                        'RingInfo.cpp',
                                        'EditableMol.cpp',
                                        ],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs),
               Extension("Chem.rdmolops",['rdmolops.cpp',
                                          'MolOps.cpp',
                                          ],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs),
               Extension("Chem.rdmolfiles",["rdmolfiles.cpp",
                                        # MolSupplier stuff
                                        "SDMolSupplier.cpp",
                                        "TDTMolSupplier.cpp",
                                        "SmilesMolSupplier.cpp",
                                        #Mol Writer stuff
                                        "SmilesWriter.cpp",
                                        "SDWriter.cpp",
                                        "TDTWriter.cpp",
                                        ],
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs),
         ])
      
