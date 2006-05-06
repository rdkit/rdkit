# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from RDBuild import *
from distutils.core import setup,Extension
import RDConfig


libDirs.extend([
         os.path.join(RDConfig.RDBaseDir,'Code','ML','InfoTheory'),
         rdBitDir,
         ])

libraries.extend(["InfoTheory", "DataStructs"])

setup(name="ML.InfoTheory.rdInfoTheory", version="2.0",
      ext_modules=[Extension("ML.InfoTheory.rdInfoTheory",
                             ["InfoBitRanker.cpp",
                              "BitCorrMatGenerator.cpp",
                              "rdInfoTheory.cpp"],
                             library_dirs=libDirs,
                             libraries=libraries,
                             extra_compile_args=compileArgs)])
