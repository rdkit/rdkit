# Installation

Below a number of installation recipes is presented, with varying degree of complexity.

## Cross-platform under anaconda python (fastest install)

### Introduction to anaconda

Conda is an open-source, cross-platform, software package manager. It supports the packaging and distribution of software components, and manages their installation inside isolated execution environments. It has several analogies with pip and virtualenv, but it is designed to be more "python-agnostic" and more suitable for the distribution of binary packages and their dependencies.

### How to get conda

The easiest way to get Conda is having it installed as part of the [Anaconda Python distribution](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). A possible (but a bit more complex to use) alternative is provided with the smaller and more self-contained [Miniconda](https://conda.io/miniconda.html). The conda source code repository is available on [github](https://github.com/conda) and additional documentation is provided by the project [website](https://conda.io/docs/).

### How to install RDKit with Conda

Creating a new conda environment with the RDKit installed requires one single command similar to the following::

```shellsession
$ conda create -c rdkit -n my-rdkit-env rdkit
```

Finally, the new environment must be activated so that the corresponding python interpreter becomes available in the same shell:

```shellsession
$ conda activate my-rdkit-env
```

If for some reason this does not work, try:

```shellsession
$ cd [anaconda folder]/bin
$ source activate my-rdkit-env
```

Windows users will use a slightly different command:

```
C:\> activate my-rdkit-env
```

#### conda-forge package

A [conda-forge](https://conda-forge.org/#about) RDKit package is also available, which may offer an easier installation for users already using other conda-forge packages. This package can be installed with:

```shellsession
$ conda install -c conda-forge rdkit
```

### How to build from source with Conda

For more details on building from source with Conda, see the [conda-rdkit repository](https://github.com/rdkit/conda-rdkit).

#### macOS 10.12 (Sierra): Python 3 environment

The following commands will create a development environment for macOS Sierra and Python 3. Download
Miniconda3-latest-MacOSX-x86_64.sh from [Conda](http://conda.pydata.org/miniconda.html) and run these
following commands:

```
bash Miniconda3-latest-MacOSX-x86_64.sh
conda install numpy matplotlib
conda install cmake cairo pillow eigen pkg-config
conda install boost-cpp boost py-boost
```

Optionally, add the following packages to your environment as useful development tools.

```
pip install yapf==0.11.1
pip install coverage==3.7.1
```

Then follow the usual build instructions. The `PYTHON_INCLUDE_DIR` must be set in the
cmake command.

```
PYROOT=<path to miniconda3>
cmake -DPYTHON_INCLUDE_DIR=$PYROOT/include/python3.6m  \
  -DRDK_BUILD_AVALON_SUPPORT=ON \
  -DRDK_BUILD_CAIRO_SUPPORT=ON \
  -DRDK_BUILD_INCHI_SUPPORT=ON \
  ..
```

Once `make` and `make install` completed successfully, use the following command to run the tests:

```
RDBASE=$RDBASE DYLD_FALLBACK_LIBRARY_PATH="$RDBASE/lib:$PYROOT/lib" PYTHONPATH=$RDBASE ctest
```

This is required due to the [System Integrity Protection SIP](https://en.wikipedia.org/wiki/System_Integrity_Protection)
introduced in more recent macOS versions.

#### Linux x86_64: Python 3 environment

The following commands will create a development environment for Linux x86_64 and Python 3.

Start by downloading the latest anaconda installer from [Anaconda](https://www.anaconda.com/download/#linux) and install it. Then, install the required packages:

```
bash Anaconda3-5.2.0-x86_64.sh
conda install -y cmake cairo pillow eigen pkg-config
conda install -y boost-cpp boost py-boost
```

Numpy and matplotlib are already part of the base installation of anaconda. Due to the latest boost libraries being currently built with a GLIBC version higher than the default in anaconda, we need to update to a more recent version:

```
conda install -y gxx_linux-64
```

At this point, you should be able to clone the RDKit repository to the desired build location, and start the build. Please consider that it is necessary to indicate the path to the numpy headers for RDKit to find them, since anaconda hides them inside the numpy package:

```
git clone https://github.com/rdkit/rdkit.git
cd rdkit
mkdir build && cd build
cmake -DPy_ENABLE_SHARED=1 \
  -DRDK_INSTALL_INTREE=ON \
  -DRDK_INSTALL_STATIC_LIBS=OFF \
  -DRDK_BUILD_CPP_TESTS=ON \
  -DPYTHON_NUMPY_INCLUDE_PATH="$(python -c 'import numpy ; print(numpy.get_include())')" \
  -DBOOST_ROOT="$CONDA_PREFIX" \
  ..
```

And finally, `make`, `make install` and `ctest`


### Installing and using PostgreSQL and the RDKit PostgreSQL cartridge from a conda environment

Due to the conda python distribution being a different version to the system python, it is easiest to install PostgreSQL and the PostgreSQL python client via conda.

With your environment activated, this is done simply by:

```
conda install -c rdkit rdkit-postgresql
```

The conda packages PostgreSQL version needs to be initialized by running the initdb command found in `[conda folder]/envs/my-rdkit-env/bin`

```
[conda folder]/envs/my-rdkit-env/bin/initdb -D /folder/where/data/should/be/stored
```

PostgreSQL can then be run from the terminal with the command:

```
[conda folder]/envs/my-rdkit-env/bin/postgres -D /folder/where/data/should/be/stored
```

For most use cases you will instead need to run PostgreSQL as a daemon, one way to do this is using supervisor. You can find out more and how to install supervisor [here](http://supervisord.org/). The required configuration file will look something like this:

```
[program:postgresql]
command=[conda folder]/envs/my-rdkit-env/bin/postgres -D /folder/where/data/should/be/stored
user=[your username]
autorestart=true
```

Once PostgreSQL is up and running, all of the normal PostgreSQL commands can then be run when your conda environment is activated. Therefore to create a database you can run:

```
createdb my_rdkit_db
psql my_rdkit_db
# create extension rdkit;
```

If you are trying to use multiple installations of PostgreSQL in different environments, you will need to setup different pid files, unix sockets and ports by [editing the PostgreSQL config files](https://opensourcedbms.com/dbms/running-multiple-postgresql-9-2-instances-on-one-server-in-centos-6rhel-6fedora/). With the above configurations these files can be found in /folder/where/data/should/be/stored.

## Linux and OS X

### Installation from repositories

#### Ubuntu 12.04 and later

Thanks to the efforts of the Debichem team, RDKit is available via the Ubuntu repositories. To install:

```shellsession
$ sudo apt-get install python-rdkit librdkit1 rdkit-data
```

#### Fedora, CentOS, and RHEL

Thanks to Gianluca Sforna's work, binary RPMs for the RDKit are now part of the official Fedora repositories:
https://admin.fedoraproject.org/pkgdb/package/rpms/rdkit/

#### OS X

Eddie Cao has produced a homebrew formula that can be used to easily build the RDKit [https://github.com/rdkit/homebrew-rdkit](https://github.com/rdkit/homebrew-rdkit)

### Building from Source

Starting with the `2018_03` release, the RDKit core C++ code is written in modern C++; for this release that means C++11.
This means that the compilers used to build it cannot be completely ancient. Here are the minimum tested versions:

- g++ v4.8: though note that the SLN parser code cannot be built with v4.8. It will automatically be disabled when this older compiler is used.
- clang v3.9: it may be that older versions of the compiler also work, but we haven't tested them.
- Visual Studio 2015: it may be that older versions of the compiler also work, but we haven't tested them.

#### Installing prerequisites from source

-   Required packages:
    - cmake. You need version 3.1 (or more recent). http://www.cmake.org if your linux distribution doesn't have an appropriate package.
    - The following are required if you are planning on using the Python wrappers
        -   The python headers. This probably means that you need to install the python-dev package (or whatever it's called) for your linux distribution.
        -   sqlite3. You also need the shared libraries. This may require that you install a sqlite3-dev package.
        -   You need to have [numpy](http://www.scipy.org/NumPy) installed.

> **note**
>
> for building with XCode4 on OS X there seems to be a problem with the version of numpy that comes with XCode4. Please see below in the (see faq) section for a workaround.


###### Installing Boost

If your linux distribution has a boost-devel package with a version >= 1.58 including the python and serialization libraries, you can use that and save yourself the steps below.

> **note**
>
> if you *do* have a version of the boost libraries pre-installed and you want to use your own version, be careful when you build the code. We've seen at least one example on a Fedora system where cmake compiled using a user-installed version of boost and then linked against the system version. This led to segmentation faults. There is a workaround for this below in the (see FAQ) section.

-   download the boost source distribution from [the boost web site](http://www.boost.org)
-   extract the source somewhere on your machine (e.g. `/usr/local/src/boost_1_58_0`)
-   build the required boost libraries. The boost site has [detailed instructions](http://www.boost.org/doc/libs/1_58_0/more/getting_started/index.html) for this, but here's an overview:
    -   `cd $BOOST`
    -   If you want to use the python wrappers: `./bootstrap.sh --with-libraries=python,serialization`
    -   If not using the python wrappers: `./bootstrap.sh --with-libraries=serialization`
    -   `./b2 install`

If you have any problems with this step, check the boost [installation instructions](http://www.boost.org/more/getting_started/unix-variants.html).

#### Building the RDKit

Fetch the source, here as tar.gz but you could use git as well:

```shellsession
$ wget https://github.com/rdkit/rdkit/archive/Release_XXXX_XX_X.tar.gz
```

-   Ensure that the prerequisites are installed
-   environment variables:
    -   `RDBASE`: the root directory of the RDKit distribution (e.g. `~/RDKit`)
    -   *Linux:* `LD_LIBRARY_PATH`: make sure it includes `$RDBASE/lib` and wherever the boost shared libraries were installed
    -   *OS X:* `DYLD_LIBRARY_PATH`: make sure it includes `$RDBASE/lib` and wherever the boost shared libraries were installed
    - The following are required if you are planning on using the Python wrappers:
        -   `PYTHONPATH`: make sure it includes `$RDBASE`
-   Building:
    -   `cd $RDBASE`
    -   `mkdir build`
    -   `cd build`
    -   `cmake ..` : See the section below on configuring the build if you need to specify a non-default version of python or if you have boost in a non-standard location
    -   `make` : this builds all libraries, regression tests, and wrappers (by default).
    -   `make install`

See below for a list of FAQ and solutions.

#### Testing the build (optional, but recommended)

-   `cd $RDBASE/build` and do `ctest`
-   you're done!

#### Advanced

##### Specifying install location

You need to turn `RDK_INSTALL_INTRE` off:

```
cmake -DRDK_INSTALL_INTREE=OFF -DCMAKE_INSTALL_PREFIX=/path/as/you/like ..
```

##### Specifying an alternate Boost installation

You need to tell cmake where to find the boost libraries and header files:

If you have put boost in `/opt/local`, the cmake invocation would look like:

```
cmake -DBOOST_ROOT=/opt/local ..
```

*Note* that if you are using your own boost install on a system with a system install, it's normally a good idea to also include the argument `-D Boost_NO_SYSTEM_PATHS=ON` in your cmake command.

##### Specifying an alternate Python installation

If you aren't using the default python installation for your computer, You need to tell cmake where to find the python library it should link against and the python header files.

Here's a sample command line:

```
cmake -D PYTHON_LIBRARY=/usr/lib/python3.6/config/libpython3.6.a -D PYTHON_INCLUDE_DIR=/usr/include/python3.6/ -D PYTHON_EXECUTABLE=/usr/bin/python3 ..
```

The `PYTHON_EXECUTABLE` part is optional if the correct python is the first version in your `PATH`.

##### Disabling the Python wrappers

You can completely disable building of the python wrappers:

```
cmake -DRDK_BUILD_PYTHON_WRAPPERS=OFF ..
```

##### Recommended extras

- You can enable support for generating InChI strings and InChI keys by adding the argument `-DRDK_BUILD_INCHI_SUPPORT=ON` to your cmake command line.
- You can enable support for the Avalon toolkit by adding the argument `-DRDK_BUILD_AVALON_SUPPORT=ON` to your cmake command line.
- If you'd like to be able to generate high-quality PNGs for structure depiction, you should have cairo installed on your system and build the RDKit with cairo support enabled: `-DRDK_BUILD_CAIRO_SUPPORT=ON`
- If you'd like to be able to use the 3D descriptors, you need to have a copy of eigen3 installed. Most operating systems have an appropriate package.

##### Building the Java wrappers

*Additional Requirements*

-   SWIG >v2.0: http://www.swig.org

*Building*

-   When you invoke cmake add `-D RDK_BUILD_SWIG_WRAPPERS=ON` to the arguments. For example: `cmake -D RDK_BUILD_SWIG_WRAPPERS=ON ..`
-   Build and install normally using make. The directory `$RDBASE/Code/JavaWrappers/gmwrapper` will contain the three required files: `libGraphMolWrap.so` (`libGraphMolWrap.jnilib` on OS X), `org.RDKit.jar`, and `org.RDKitDoc.jar`.

*Using the wrappers*

To use the wrappers, the three files need to be in the same directory, and that should be on your `CLASSPATH` and in the `java.library.path`. An example using jython:

```
$ CLASSPATH=$CLASSPATH:$RDBASE/Code/JavaWrappers/gmwrapper/org.RDKit.jar jython -Djava.library.path=$RDBASE/Code/JavaWrappers/gmwrapper
Jython 2.2.1 on java1.6.0_20
Type "copyright", "credits" or "license" for more information.
>>> from org.RDKit import *
>>> from java import lang
>>> lang.System.loadLibrary('GraphMolWrap')
>>> m = RWMol.MolFromSmiles('c1ccccc1')
>>> m.getNumAtoms()
6L
```

##### Optional packages
-   If you would like to install the RDKit InChI support, follow the instructions in `$RDBASE/External/INCHI-API/README`.
-   If you would like to install the RDKit Avalon toolkit support, follow the instructions in `$RDBASE/External/AvalonTool/README`.
-   If you would like to build and install the PostgreSQL cartridge, follow the instructions in `$RDBASE/Code/PgSQL/rdkit/README`.

#### Frequently Encountered Problems

In each case I've replaced specific pieces of the path with `...`.

*Problem:* :

```
Linking CXX shared library libSLNParse.so
/usr/bin/ld: .../libboost_regex.a(cpp_regex_traits.o): relocation R_X86_64_32S against `std::basic_string<char, std::char_traits<char>, std::allocator<char> >::_Rep::_S_empty_rep_storage' can not be used when making a shared object; recompile with -fPIC
.../libboost_regex.a: could not read symbols: Bad value
collect2: ld returned 1 exit status
make[2]: *** [Code/GraphMol/SLNParse/libSLNParse.so] Error 1
make[1]: *** [Code/GraphMol/SLNParse/CMakeFiles/SLNParse.dir/all] Error 2
make: *** [all] Error 2
```

*Solution:*

Add this to the arguments when you call cmake: `-DBoost_USE_STATIC_LIBS=OFF`

More information here: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01119.html

* * * * *

*Problem:* :

```
.../Code/GraphMol/Wrap/EditableMol.cpp:114:   instantiated from here
.../boost/type_traits/detail/cv_traits_impl.hpp:37: internal compiler error: in make_rtl_for_nonlocal_decl, at cp/decl.c:5067

    Please submit a full bug report, with preprocessed source if appropriate. See \<URL:<http://bugzilla.redhat.com/bugzilla>\> for instructions. Preprocessed source stored into /tmp/ccgSaXge.out file, please attach this to your bugreport. make[2]: **\* [Code/GraphMol/Wrap/CMakeFiles/rdchem.dir/EditableMol.cpp.o] Error 1 make[1]:**\* [Code/GraphMol/Wrap/CMakeFiles/rdchem.dir/all] Error 2 make: *\** [all] Error 2
```

*Solution:*

Add `#define BOOST_PYTHON_NO_PY_SIGNATURES` at the top of `Code/GraphMol/Wrap/EditableMol.cpp`

More information here: http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01178.html

* * * * *

*Problem:*

Your system has a version of boost installed in `/usr/lib`, but you would like to force the RDKit to use a more recent one.

*Solution:*

This can be solved by providing the `-D Boost_NO_SYSTEM_PATHS=ON` argument:

```
cmake -D BOOST_ROOT=/usr/local -D Boost_NO_SYSTEM_PATHS=ON ..
```

* * * * *

*Problem:*

Building on OS X with XCode 4

The problem seems to be caused by the version of numpy that is distributed with XCode 4, so you need to build a fresh copy.

*Solution:* Get a copy of numpy and build it like this as root: as root:

```
export MACOSX_DEPLOYMENT_TARGET=10.6
export LDFLAGS="-Wall -undefined dynamic_lookup -bundle -arch x86_64"
export CFLAGS="-arch x86_64"
ln -s /usr/bin/gcc /usr/bin/gcc-4.2
ln -s /usr/bin/g++ /usr/bin/g++-4.2
python setup.py build
python setup.py install
```

Be sure that the new numpy is used in the build:

```
PYTHON_NUMPY_INCLUDE_PATH /Library/Python/2.6/site-packages/numpy/core/include
```

and is at the beginning of the PYTHONPATH:

```
export PYTHONPATH="/Library/Python/2.6/site-packages:$PYTHONPATH"
```

Now it's safe to build boost and the RDKit.

## Windows

### Prerequisites

-   python 3.6+ (from http://www.python.org/)
-   numpy (from http://numpy.scipy.org/ or use `pip install numpy`). Binaries for win64 are available here: http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
-   Pillow: (from https://python-pillow.github.io/> or use `pip install Pillow`)

#### Recommended extras

-   matplotlib: a library for scientific plotting from Python. <http://matplotlib.sourceforge.net/>
-   ipython : a very useful interactive shell (and much more) for Python. <http://ipython.scipy.org/dist/>


### Installation from source

#### Extra software to install

-   Microsoft Visual C++ : The Community version has everything necessary and can be downloaded for free (<https://www.visualstudio.com/vs/community>). This is a big installation and will take a while. The RDKit has been built with Visual Studio 2015 and 2017. More recent versions should be fine.
-   cmake : (<http://www.cmake.org/cmake/resources/software.html>) should be installed.
-   boost : It is strongly recommended to download and use a precompiled version of the boost libraries from <http://sourceforge.net/projects/boost/files/boost-binaries/> . When you run the installer, the only binary libraries you need are python and serialization. If you want to install boost from source, download a copy from <http://www.boost.org> and follow the instructions in the "Getting Started" section of the documentation. Make sure the libraries and headers are installed to C:\boost
-   a git client : *This is only necessary if you are planning on building development versions of the RDKit.* This can be downloaded from <http://git-scm.com/downloads>; git is also included as an optional add-on of Microsoft Visual Studio 2015.

#### Setup and Preparation

This section assumes that python is installed in `C:\Python36 that the boost libraries have been installed to `C:\boost`, and that you will build the RDKit from a directory named `C:\RDKit`. If any of these conditions is not true, just change the corresponding paths.

-   If you install things in paths that have spaces in their names, be sure to use quotes properly in your environment variable definitions.
-   If you are planning on using a development version of the RDKit: get a copy of the current RDKit source using git. If you're using the command-line client the command is: `git clone  https://github.com/rdkit/rdkit.git C:\RDKit`
-   If you are planning on using a released version of the RDKit: get a copy of the most recent release and extract it into the directory `C:\RDKit`
-   Set the required environment variables:
    -   `RDBASE = C:\RDKit`
    -   Make sure `C:\Python36 is in your PATH
    -   Make sure `C:\RDKit\lib` is in your PATH
    -   Make sure `C:\boost\lib` is in your PATH.
    -   Make sure `C:\RDKit is` in your PYTHONPATH

#### Building from the command line (recommended)

-   Create a directory `C:\RDKit\build` and cd into it
-   Run cmake. Here's an example basic command line for 64bit windows that will download the InChI and Avalon toolkit sources from the InChI Trust and SourceForge repositories, respectively, and build the PostgreSQL cartridge for the installed version of PostgreSQL:  
  `cmake -DRDK_BUILD_PYTHON_WRAPPERS=ON -DBOOST_ROOT=C:/boost -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_BUILD_AVALON_SUPPORT=ON -DRDK_BUILD_PGSQL=ON -DPostgreSQL_ROOT="C:\Program Files\PostgreSQL\9.5" -G"Visual Studio 14 2015 Win64" ..`
-   Build the code. Here's an example command line:  
  `C:/Windows/Microsoft.NET/Framework64/v4.0.30319/MSBuild.exe /m:4 /p:Configuration=Release INSTALL.vcxproj`
-   If you have built in PostgreSQL support, you will need to open a shell with administrator privileges, stop the PostgreSQL service, run the `pgsql_install.bat` installation script, then restart the PostgreSQL service (please refer to `%RDBASE%\Code\PgSQL\rdkit\README` for further details):
    -   `"C:\Program Files\PostgreSQL\9.5\bin\pg_ctl.exe" -N "postgresql-9.5" -D "C:\Program Files\PostgreSQL\9.5\data" -w stop`
    -   `C:\RDKit\build\Code\PgSQL\rdkit\pgsql_install.bat`
    -   `"C:\Program Files\PostgreSQL\9.5\bin\pg_ctl.exe" -N "postgresql-9.5" -D "C:\Program Files\PostgreSQL\9.5\data" -w start`
    -   Before restarting the PostgreSQL service, make sure that the Boost libraries the RDKit was built against are in the system PATH, or PostgreSQL will fail to create the `rdkit` extension with a deceptive error message such as:  
      `ERROR: could not load library "C:/Program Files/PostgreSQL/9.5/lib/rdkit.dll": The specified module could not be found.`


#### Testing the Build (optional, but recommended)

-   cd to `C:\RDKit\build` and run ctest. Please note that if you have built in PostgreSQL support, the current logged in user needs to be a PostgreSQL user with database creation and superuser privileges, or the PostgreSQL test will fail. A convenient option to authenticate will be to set the `PGPASSWORD` environment variable to the PostgreSQL password of the current logged in user in the shell from which you are running ctest.
-   You're done!

## License

This document is copyright (C) 2012-2020 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License. To view a copy of this license, visit <http://creativecommons.org/licenses/by-sa/4.0/> or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple words: "Do whatever you want with it, but please give us some credit."
