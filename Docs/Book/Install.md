# Installation

Below a number of installation recipies is presented, with varying degree of complexity.

## Cross-platform under anaconda python (fastest install)

### Introduction to anaconda

Conda is an open-source, cross-platform, software package manager. It supports the packaging and distribution of software components, and manages their installation inside isolated execution environments. It has several analogies with pip and virtualenv, but it is designed to be more "python-agnostic" and more suitable for the distribution of binary packages and their dependencies.

### How to get conda

The easiest way to get Conda is having it installed as part of the [Anaconda Python distribution](http://docs.continuum.io/anaconda/install.html). A possible (but a bit more complex to use) alternative is provided with the smaller and more self-contained [Miniconda](http://conda.pydata.org/miniconda.html). The conda source code repository is available on [github](https://github.com/conda) and additional documentation is provided by the project [website](http://conda.pydata.org/).

### How to install RDKit with Conda

Creating a new conda environment with the RDKit installed using these  packages requires one single command similar to the following::

  $ conda create -c https://conda.anaconda.org/rdkit -n my-rdkit-env rdkit

Finally, the new environment must be activated, so that the corresponding python interpreter becomes available in the same shell:

  $ source activate my-rdkit-env

If for some reason this does not work, try:

  $ cd [anaconda folder]/bin
  $ source activate my-rdkit-env

Windows users will use a slightly different command:

  C:\> activate my-rdkit-env

### How to build from source with Conda

For more details on building from source with Conda, see the [conda-rdkit repository](https://github.com/rdkit/conda-rdkit)

### Installing and using PostgreSQL and the RDKit PostgreSQL cartridge from a conda environment

Due to the conda python distribution being a different version to the system python, it is easiest to install PostgreSQL and the PostgreSQL python client via conda.

With your envirinment activated, this is done simply by:

    conda install -c https://conda.binstar.org/rdkit rdkit-postgresql

The conda packages PostgreSQL version needs to be initialized by running the initdb command found in [conda folder]/envs/my-rdkit-env/bin

    [conda folder]/envs/my-rdkit-env/bin/initdb -D /folder/where/data/should/be/stored

PostgreSQL can then be run from the terminal with the command:

    [conda folder]/envs/my-rdkit-env/bin/postgres -D /folder/where/data/should/be/stored

For most use cases you will instead need to run PostgreSQL as a daemon, one way to do this is using supervisor. You can find out more and how to install supervisor [here](http://supervisord.org/). The required configuration file will look something like this:

    [program:postgresql]
    command=[conda folder]/envs/my-rdkit-env/bin/postgres -D /folder/where/data/should/be/stored
    user=[your username]
    autorestart=true

Once PostgreSQL is up and running, all of the normal PostgreSQL commands can then be run when your conda environment is activated. Therefore to create a database you can run:

    createdb my_rdkit_db
    psql my_rdkit_db
    # create extension rdkit;

If you are trying to use multiple installations of PostgreSQL in different environments, you will need to setup different pid files, unix sockets and ports by [editing the PostgreSQL config files](https://opensourcedbms.com/dbms/running-multiple-postgresql-9-2-instances-on-one-server-in-centos-6rhel-6fedora/). With the above configurations these files can be found in /folder/where/data/should/be/stored.

## Linux and the Mac

### Installation from repositories

#### Ubuntu 12.04 and later

Thanks to the efforts of the Debichem team, RDKit is available via the Ubuntu repositories. To install:

    sudo apt-get install python-rdkit librdkit1 rdkit-data

#### Fedora, CentOS, and RHEL

Gianluca Sforna creates binary RPMs that can be found here: [http://giallu.fedorapeople.org/rdkit-20XX.XX/](http://giallu.fedorapeople.org/rdkit-20XX.XX/)

#### MacOS

Eddie Cao has produced a homebrew formula that can be used to easily build the RDKit [https://github.com/rdkit/homebrew-rdkit](https://github.com/rdkit/homebrew-rdkit)

### Building from Source

#### Prerequisites

##### Installing prerequisites as packages

###### Ubuntu and other debian-derived systems

Install the following packages using apt-get:

    flex bison build-essential python-numpy cmake python-dev sqlite3 libsqlite3-dev
    libboost-dev libboost-python-dev libboost-regex-dev

###### Fedora, CentOS (5.7+), and RHEL

Install the following packages using yum:

    cmake tk-devel readline-devel zlib-devel bzip2-devel sqlite-devel @development-tools

Packages to install from source (not required on RHEL/CentOS 6.x):

-   python 2.7 : use `./configure CFLAGS=-fPIC --enable-unicode=ucs4 --enable-shared`
-   numpy : do `export LD\_LIBRARY\_PATH="/usr/local/lib"` before `python setup.py install`
-   boost 1.48.0 or later: do `./bootstrap.sh --with-libraries=python,regex; ./b2; ./b2 install`

###### Older versions of CentOS

Here things are more difficult. Check this wiki page for information: https://code.google.com/p/rdkit/wiki/BuildingOnCentOS

##### Installing prerequisites from source

-   Required packages:
  -   cmake. You need version 2.6 (or more recent). http://www.cmake.org if your linux distribution doesn't have an appropriate package.

         > **note**
         >
         > It seems that v2.8 is a better bet than v2.6. It might be worth compiling your own copy of v2.8 even if v2.6 is already installed.

  - The following are required if you are planning on using the Python wrappers
      -   The python headers. This probably means that you need to install the python-dev package (or whatever it's called) for your linux distribution.
      -   sqlite3. You also need the shared libraries. This may require that you install a sqlite3-dev package.
      -   You need to have numpy (http://www.scipy.org/NumPy) installed.
          > **note**
          >
          > for building with XCode4 on the MacOS – there seems to be a problem with the version of numpy that comes with XCode4. Please see below in the (see faq) section for a workaround.
 -   Optional packages
     -   If you would like to install the RDKit InChI support (first available in the Q2 2011 release), follow the instructions in $RDBASE/External/INCHI-API to get a copy of the InChI source and put it in the appropriate place.

###### Installing Boost

If your linux distribution has a boost-devel package including the python and regex libraries, you can use that and save yourself the steps below.

> **note**
>
> if you *do* have a version of the boost libraries pre-installed and you want to use your own version, be careful when you build the code. We've seen at least one example on a Fedora system where cmake compiled using a user-installed version of boost and then linked against the system version. This led to segmentation faults. There is a workaround for this below in the (see FAQ) section.

-   download the boost source distribution from [the boost web site](http://www.boost.org)
-   extract the source somewhere on your machine (e.g. `/usr/local/src/boost_1_45_0`)
-   build the required boost libraries:
   -   `cd $BOOST`
   -   If you want to use the python wrappers: `./bootstrap.sh --with-libraries=python,regex`
   -   If not using the python wrappers: `./bootstrap.sh --with-libraries=regex`
   -   Building on 32 bit systems: `./b2 install`
   -   Building on 64 bit systems: `./b2 address-model=64 cflags=-fPIC cxxflags=-fPIC install`

     If you have any problems with this step, check the boost [installation instructions](http://www.boost.org/more/getting_started/unix-variants.html).

#### Building the RDKit

Fetch the source, here as tar.gz but you could use git as well:

    wget http://downloads.sourceforge.net/project/rdkit/rdkit/QX_20XX/RDKit_20XX_XX_X.tgz

-   Ensure that the prerequisites are installed
-   environment variables:
     -   RDBASE: the root directory of the RDKit distribution (e.g. ~/RDKit)
     -   *Linux:* LD_LIBRARY_PATH: make sure it includes $RDBASE/lib and wherever the boost shared libraries were installed
     -   *Mac:* DYLD_LIBRARY_PATH: make sure it includes $RDBASE/lib and wherever the boost shared libraries were installed
     - The following are required if you are planning on using the Python wrappers:  
       -   PYTHONPATH: make sure it includes $RDBASE
-   Building:
   -   cd to $RDBASE
   -   `mkdir build`
   -   `cd build`
   -   `cmake ..` : See the section below on configuring the build if you need to specify a non-default version of python or if you have boost in a non-standard location
   -   `make` : this builds all libraries, regression tests, and wrappers (by default).
   -   `make install`

See below for a list of FAQ and solutions.

#### Testing the build (optional, but recommended)

-   cd to $RDBASE/build and do `ctest`
-   you're done!

#### Advanced

##### Specifying an alternate Boost installation

You need to tell cmake where to find the boost libraries and header files:

If you have put boost in /opt/local, the cmake invocation would look like:

    cmake -DBOOST_ROOT=/opt/local ..

##### Specifying an alternate Python installation

You need to tell cmake where to find the python library it should link against and the python header files.

Here's a sample command line:

    cmake -D PYTHON_LIBRARY=/usr/lib/python2.5/config/libpython2.5.a -D PYTHON_INCLUDE_DIR=/usr/include/python2.5/ -D PYTHON_EXECUTABLE=/usr/bin/python ..

The `PYTHON_EXECUTABLE` part is optional if the correct python is the first version in your PATH.

##### Disabling the Python wrappers

You can completely disable building of the python wrappers by setting the configuration variable RDK_BUILD_PYTHON_WRAPPERS to nil:

    cmake -D RDK_BUILD_PYTHON_WRAPPERS= ..

#### Recommended extras

-   for structure depiction cairo (for use with Python2) or cairocffi (for use with Python3) and their respective Python bindings are recommended.

##### Building the Java wrappers

*Additional Requirements*

-   SWIG v2.0.x: http://www.swig.org
-   Junit: get a copy of `http://search.maven.org/remotecontent?filepath=junit/junit/4.10/junit-4.10.jar` and put it in the directory `$RDBASE/External/java_lib` (you will need to create the directory) and rename it to junit.jar.

*Building*

-   When you invoke cmake add `-D RDK_BUILD_SWIG_WRAPPERS=ON` to the arguments. For example: `cmake -D RDK_BUILD_SWIG_WRAPPERS=ON ..`
-   Build and install normally using make. The directory `$RDBASE/Code/JavaWrappers/gmwrapper` will contain the three required files: libGraphMolWrap.so (libGraphMolWrap.jnilib on the Mac), org.RDKit.jar, and org.RDKitDoc.jar.

*Using the wrappers*

To use the wrappers, the three files need to be in the same directory, and that should be on your CLASSPATH and in the java.library.path. An example using jython:

    % CLASSPATH=$CLASSPATH:$RDBASE/Code/JavaWrappers/gmwrapper/org.RDKit.jar; jython -Djava.library.path=$RDBASE/Code/JavaWrappers/gmwrapper
    Jython 2.2.1 on java1.6.0_20
    Type "copyright", "credits" or "license" for more information.
    >>> from org.RDKit import *
    >>> from java import lang
    >>> lang.System.loadLibrary('GraphMolWrap')
    >>> m = RWMol.MolFromSmiles('c1ccccc1')
    >>> m.getNumAtoms()
    6L

#### Frequently Encountered Problems

In each case I've replaced specific pieces of the path with `...`.

*Problem:* :

    Linking CXX shared library libSLNParse.so
    /usr/bin/ld: .../libboost_regex.a(cpp_regex_traits.o): relocation R_X86_64_32S against `std::basic_string<char, std::char_traits<char>, std::allocator<char> >::_Rep::_S_empty_rep_storage' can not be used when making a shared object; recompile with -fPIC
    .../libboost_regex.a: could not read symbols: Bad value
    collect2: ld returned 1 exit status
    make[2]: *** [Code/GraphMol/SLNParse/libSLNParse.so] Error 1
    make[1]: *** [Code/GraphMol/SLNParse/CMakeFiles/SLNParse.dir/all] Error 2
    make: *** [all] Error 2

*Solution:*

Add this to the arguments when you call cmake: `-DBoost_USE_STATIC_LIBS=OFF`

More information here: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01119.html

* * * * *

*Problem:* :

    .../Code/GraphMol/Wrap/EditableMol.cpp:114:   instantiated from here
    .../boost/type_traits/detail/cv_traits_impl.hpp:37: internal compiler error: in make_rtl_for_nonlocal_decl, at cp/decl.c:5067

    Please submit a full bug report, with preprocessed source if appropriate. See \<URL:<http://bugzilla.redhat.com/bugzilla>\> for instructions. Preprocessed source stored into /tmp/ccgSaXge.out file, please attach this to your bugreport. make[2]: **\* [Code/GraphMol/Wrap/CMakeFiles/rdchem.dir/EditableMol.cpp.o] Error 1 make[1]:**\* [Code/GraphMol/Wrap/CMakeFiles/rdchem.dir/all] Error 2 make: *\** [all] Error 2

*Solution:*

Add `#define BOOST_PYTHON_NO_PY_SIGNATURES` at the top of `Code/GraphMol/Wrap/EditableMol.cpp`

More information here: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01178.html

* * * * *

*Problem:*

Your system has a version of boost installed in /usr/lib, but you would like to force the RDKit to use a more recent one.

*Solution:*

This can be solved by using cmake version 2.8.3 (or more recent) and providing the `-D Boost_NO_SYSTEM_PATHS=ON` argument:

    cmake -D BOOST_ROOT=/usr/local -D Boost_NO_SYSTEM_PATHS=ON ..

* * * * *

*Problem:*

Building on the Mac with XCode 4

The problem seems to be caused by the version of numpy that is distributed with XCode 4, so you need to build a fresh copy.

*Solution:* Get a copy of numpy and build it like this as root: as root:

    export MACOSX_DEPLOYMENT_TARGET=10.6
    export LDFLAGS="-Wall -undefined dynamic_lookup -bundle -arch x86_64"
    export CFLAGS="-arch x86_64"
    ln -s /usr/bin/gcc /usr/bin/gcc-4.2
    ln -s /usr/bin/g++ /usr/bin/g++-4.2
    python setup.py build
    python setup.py install

Be sure that the new numpy is used in the build:

    PYTHON_NUMPY_INCLUDE_PATH /Library/Python/2.6/site-packages/numpy/core/include

and is at the beginning of the PYTHONPATH:

    export PYTHONPATH="/Library/Python/2.6/site-packages:$PYTHONPATH"

Now it's safe to build boost and the RDKit.

## Windows

### Prerequisites

-   Python 2.7 (from http://www.python.org/)
-   numpy (from http://numpy.scipy.org/ or use `pip install numpy`). Binaries for win64 are available here: http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
-   Pillow: (from https://python-pillow.github.io/> or use `pip install Pillow`)

#### Recommended extras

-   aggdraw: a library for high-quality drawing in Python. Instructions for downloading are here: <http://effbot.org/zone/aggdraw-index.htm> The new (as of May 2008) drawing code has been tested with v1.2a3 of aggdraw. Despite the alpha label, the code is stable and functional.
-   matplotlib: a library for scientific plotting from Python. <http://matplotlib.sourceforge.net/>
-   ipython : a very useful interactive shell (and much more) for Python. <http://ipython.scipy.org/dist/>
> -   win32all: Windows extensions for Python. <http://sourceforge.net/projects/pywin32/>

### Installation of RDKit binaries

-   Get the appropriate windows binary build from: https://sourceforge.net/projects/rdkit/files/rdkit
-   Extract the zip file somewhere without a space in the name, i.e. `c:/`
-   The rest of this will assume that the installation is in `c:/RDKit_2012_12_1`
-   Set the following environment variables:
   -   RDBASE: `c:/RDKit_2012_12_1`
   -   PYTHONPATH: `%RDBASE%` if there is already a PYTHONPATH, put `;%RDBASE%` at the end.
   -   PATH: add `;%RDBASE%/lib` to the end

In Win7 systems, you may run into trouble due to missing DLLs, see one thread from the mailing list: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01632.html You can download the missing DLLs from here: http://www.microsoft.com/en-us/download/details.aspx?id=5555

### Installation from source

#### Extra software to install

-   Microsoft Visual C++ : The Express version has everything necessary and can be downloaded for free (<http://www.microsoft.com/express/download/>). This is a big installation and will take a while. The RDKit has been successfully built with all version of visual c++ since 6.0, so the current version of VC++ (2010 as of this writing) should be fine.
-   cmake : (<http://www.cmake.org/cmake/resources/software.html>) should be installed.
-   boost : It is strongly recommended to download and use a precompiled version of the boost libraries from <http://sourceforge.net/projects/boost/files/boost-binaries/> . When you run the installer, the only binary libraries you need are python, regex, and system. If you want to install boost from source, download a copy from <http://www.boost.org> and follow the instructions in the "Getting Started" section of the documentation. Make sure the libraries and headers are installed to c:boost
-   a git client : *This is only necessary if you are planning on building development versions of the RDKit.* This can be downloaded from <http://git-scm.com/downloads> .
-   Optional packages
   -   If you would like to install the RDKit InChI support, follow the instructions in $RDBASE/External/INCHI-API/README to get a copy of the InChI source and put it in the appropriate place.
   -   If you would like to install the RDKit Avalon toolkit support, follow the instructions in $RDBASE/External/AvalonTool/README to get a copy of the Avalon toolkit source and put it in the appropriate place.

#### Setup and Preparation

This section assumes that python is installed in `c:\Python27`, that the boost libraries have been installed to `c:\boost`, and that you will build the RDKit from a directory named `c:\RDKit`. If any of these conditions is not true, just change the corresponding paths.

-   If you install things in paths that have spaces in their names, be sure to use quotes properly in your environment variable definitions.
-   If you are planning on using a development version of the RDKit: get a copy of the current RDKit source using git. If you're using the command-line client the command is: `git clone  https://github.com/rdkit/rdkit.git c:\RDKit`
-   If you are planning on using a released version of the RDKit : get a copy of the most recent release and extract it into the directory `c:\RDKit`
-   Set the required environment variables (you can set this in cygwin or in windows. If you set them in windows, be sure to restart your cygwin window)
   -   `RDBASE = c:\RDKit`
   -   Make sure `c:\Python27` is in your PATH
   -   Make sure `c:\RDKit\lib` is in your PATH
   -   Make sure `c:\boost\lib` is in your PATH.
   -   Make sure `c:\RDKit is` in your PYTHONPATH

#### Building from the command line (recommended)

-   Create a directory `c:\RDKit\build` and cd into it
-   Run cmake. Here's an example basic command line for 64bit windows that assumes the InChI and Avalon toolkit sources are available (see above): `cmake -DRDK_BUILD_PYTHON_WRAPPERS=0N -DAVALONTOOLS_DIR=c:/avalontoolkit_beta/sourcedistribution -DBOOST_ROOT=c:/boost -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_BUILD_AVALON_SUPPORT=ON -G"Visual Studio 10 Win64" ..`
-   Build the code. Here's an example command line: `c:/Windows/Microsoft.NET/Framework64/v4.0.30319/MSBuild.exe /m:4 /p:Configuration=Release INSTALL.vcxproj`

#### Testing the Build (optional, but recommended)

-   cd to `c:\RDKit\build` and run ctest.
-   you're done!

## License

This document is copyright (C) 2012-2015 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 License. To view a copy of this license, visit <http://creativecommons.org/licenses/by-sa/3.0/> or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple words: “Do whatever you want with it, but please give us some credit.”
