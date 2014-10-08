# Run this with:
#  python setup.py install --install-lib=.
from __future__ import print_function
from distutils.core import setup,Extension
import RDConfig

# force the use of g++ please
from distutils import sysconfig
save_init_posix = sysconfig._init_posix
def my_init_posix():
    print('my_init_posix: changing gcc to g++')
    save_init_posix()
    g = sysconfig.get_config_vars()
    g['CC'] = 'g++'
    g['LDSHARED'] = 'g++ -shared'
    g['PY_CFLAGS']= g['PY_CFLAGS'].replace('-O3','')
sysconfig._init_posix = my_init_posix


destDir = RDConfig.RDCodeDir
extDir=RDConfig.RDBaseDir+"/External"

# our stuff
destDir = RDConfig.RDCodeDir
srcDir = RDConfig.RDBaseDir+"/Code"
extDir=RDConfig.RDBaseDir+"/External"

incDirs = [srcDir]

rdbitdir = srcDir + "/DataStructs"
libDirs=[rdbitdir]
libraries=["DataStructs"]

# this is how things are done with BPLv2
boostInc = '-isystem%s'%(extDir+"/boost_1_30_0")

# FIX: there's gotta be a better way of doing this
pyLibDir = '/usr/lib/python2.2/config'
boostLibDir=extDir+"/boost_1_30_0/libs/python/build/bin/libboost_python.so/gcc/debug/runtime-link-dynamic/shared-linkable-true/"
boostLib="boost_python"
libDirs.append(boostLibDir)
libDirs.append(pyLibDir)

libraries.append(boostLib)
libraries.append("python2.2") 

compileArgs=['-ftemplate-depth-150',
	     '-DBOOST_PYTHON_DYNAMIC_LIB',
	     boostInc,
	     ]
setup(name="crossTest",version="1.0",
  ext_modules=[Extension("moduleA",["wrapA.cpp"],
                         include_dirs=incDirs,
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs),
               Extension("moduleB",["moduleB.cpp"],
                         include_dirs=incDirs,
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs),
               Extension("moduleC",["wrapC.cpp"],
                         include_dirs=incDirs,
                         library_dirs=libDirs,
                         libraries=libraries,
                         extra_compile_args=compileArgs)])

