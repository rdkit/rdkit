# Run this with:
#  python setup.py install --install-lib=$RDBASE/Python
from distutils.core import setup,Extension
import RDConfig


destDir = RDConfig.RDCodeDir

libDirs=['.']
libraries=['hc']
setup(name="ML.Cluster.Clustering",version="1.0",
      package_dir={'':destDir},
      ext_modules=[Extension("ML.Cluster.Clustering",["Clustering.c"],
                             library_dirs=libDirs,
                             libraries=libraries)])
