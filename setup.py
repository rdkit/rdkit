#!/usr/bin/env python
#

from distutils.core import setup,Extension
from distutils import sysconfig
import os, sys, re, glob, shutil

version='2009.Q1b1'

module_ext = sysconfig.get_config_var('SO')
if sys.platform=="win32":
     install_base="Lib/site-packages"
else:
     install_base = os.path.join(sysconfig.get_config_var('LIBDEST'),'site-packages')
ext_modules=[]

child_packages = [
               "rdkit.Chem",
               "rdkit.DataManip",
               "rdkit.DataStructs",
               "rdkit.Dbase",
               "rdkit.DistanceGeometry",
               "rdkit.ForceField",
               "rdkit.Geometry",
               "rdkit.Logger",
               "rdkit.ML",
               "rdkit.Numerics",
               "rdkit.SimDivFilters",
               "rdkit.VLib",
               "rdkit.sping",
               "rdkit.utils",
               ]

sos = [(os.path.join(install_base,'rdkit'),['rdkit/rdBase'+module_ext])]

py_packages = ["rdkit"]+child_packages

for pkg in child_packages:
     for root,dirs,files in os.walk(pkg.replace('.','/')):
          if '.svn' in dirs: dirs.remove('.svn')
          if 'test_data' in dirs: dirs.remove('test_data')

          modName=root.replace(os.path.sep,'.')
          if '__init__.py' in files and modName not in py_packages:
               py_packages.append(modName)

          files=[os.path.join(root,file) for file in files if (os.path.splitext(file)[-1]==module_ext or\
                                                                                 'test_data' in root)]
          sos.extend([(os.path.join(install_base,root),files)])

extraBase='share/rdkit'


projects=[]
for root,dirs,files in os.walk('Projects'):
     if '.svn' in dirs: dirs.remove('.svn')

     files=[os.path.join(root,filen) for filen in files]
     projects.append((extraBase+'/'+root,files))

data_files = [(extraBase+'/Data',glob.glob('Data/*.*'))]
data_files.extend([(extraBase,glob.glob('./*.txt'))])
if sys.platform=='win32':
     data_files.extend([(extraBase+'/lib',glob.glob('bin/*.dll'))])
else:
     data_files.extend([(extraBase+'/lib',glob.glob('bin/*'))])
data_files.extend(sos)



documentation = []
for root,dirs,files in os.walk('Docs'):
     if '.svn' in dirs: dirs.remove('.svn')

     files=[os.path.join(root,filen) for filen in files]
     documentation.append((extraBase+'/'+root,files))

if sys.platform=='win32':
     scripts=("postinstall_win32.py",)
else:
     scripts=()


setup(
      name='rdkit',
      version=version,
      description='RDKit Cheminformatics Library',
      long_description="""Data structures, algorithms, and scripts for cheminformatics.""",
      author='Greg Landrum',
      author_email='glandrum@users.sourceforge.net',
      url='http://www.rdkit.org/',
      download_url = 'http://code.google.com/p/rdkit/downloads/list',
      license='BSD',
      platforms=['Windows','Linux','Mac OS-X'],
      scripts=scripts,
      classifiers = ['Development Status :: 5 - Production/Stable',
                     'Environment :: Console',
                     'Intended Audience :: Developers',
                     'Intended Audience :: Science/Research',
                     'Programming Language :: Python',
                     'Programming Language :: C++',
                     'License :: OSI Approved :: BSD License',
                     'Topic :: Scientific/Engineering :: Chemistry',

                     ],
      packages=py_packages,
      ext_modules=ext_modules,
      package_dir={'rdkit':'rdkit'},
      data_files=data_files+documentation+projects,
      )
