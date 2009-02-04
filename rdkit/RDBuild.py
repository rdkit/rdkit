# $Id$
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" variables and definitions for distutils builds
"""
from distutils.core import setup,Extension
import os
import RDConfig



#-------------------------------------------
# force the use of g++
from distutils import sysconfig
save_init_posix = sysconfig._init_posix
def my_init_posix():
    print 'my_init_posix: changing gcc to g++'
    save_init_posix()
    g = sysconfig._config_vars
    g['CC'] = 'g++'
    g['LDSHARED'] = 'g++ -shared'
sysconfig._init_posix = my_init_posix

# by default, distutils includes -g in the argument list.  This turns
# out to be a complete disaster with boost, so make sure we use
# the proper arguments:
cflags = os.environ.get('RDOPTFLAGS','-O2')
sysconfig.get_config_vars()['OPT']=cflags


#-------------------------------------------

destDir = RDConfig.RDCodeDir
extDir=RDConfig.RDBaseDir+"/External"
boostBase = os.environ['BOOSTBASE']
boostHome = os.environ['BOOSTHOME']
boostInc = '-I%(boostHome)s/include/%(boostBase)s'%locals()


pyName="python"+os.environ['PYTHON_VERSION']
gccVers=os.environ.get('GCCVERSION','gcc')
boostLib="boost_python-gcc%s"%(gccVers)
boostLogLib="boost_log-gcc%s"%(gccVers)

libDirs=["/usr/local/lib"]
libDirs.extend([         
	 os.path.join(os.environ['PYTHON_ROOT'],'lib',pyName,'config'),
	 os.path.join(RDConfig.RDBaseDir,'Code','RDBoost'),
	 os.path.join(RDConfig.RDBaseDir,'Code','RDGeneral'),
	 ])
libraries=["RDBoost","RDGeneral",boostLib,pyName,boostLogLib,'boost_thread-gcc%s-mt'%gccVers]

compileArgs=["-ftemplate-depth-100",
	     '-DBOOST_PYTHON_DYNAMIC_LIB',
	     boostInc,
             "-I%s"%(os.path.join(RDConfig.RDBaseDir,'Code')),
	     ]
linkArgs=['-fPIC']

# MISC RD Stuff:
rdLibDir= os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','libs')
rdBitDir= os.path.join(RDConfig.RDBaseDir,'Code','DataStructs')
rdCatalogDir= os.path.join(RDConfig.RDBaseDir,'Code','Catalogs')
rdGeometryDir= os.path.join(RDConfig.RDBaseDir,'Code','Geometry')


# Lapack++
lapackpp = os.path.join(extDir,"Lapack++")
lapackInc = "-I"+os.path.join(lapackpp,"include")
#lapackLibs=["lapack++", "lapack","lamatrix++","blas++", "blas", 
#            "m", "gfortran"]  
lapackLibs=["lapack++", "lapack","lamatrix++","blas++", "blas", 
            "m"]
if os.environ.get('RDF77LIB',''):
    lapackLibs.append(os.environ['RDF77LIB'])


# vflib
vfBase = os.path.join(extDir,"vflib-2.0")
vfInc = "-I"+os.path.join(vfBase,"include")
vfLibDir = os.path.join(vfBase,"lib")
vfLibs = ["vf"]
