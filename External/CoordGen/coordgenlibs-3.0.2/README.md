# CoordgenLibs

[![Azure_Build_Status](https://dev.azure.com/patlorton/coordgenlibs/_apis/build/status/schrodinger.coordgenlibs?branchName=master)](https://dev.azure.com/patlorton/coordgenlibs/_build/latest?definitionId=2&branchName=master)

This is **Schrödinger, Inc's** 2D coordinate generation. It was formerly
proprietary code, but is now released under the [BSD
license](https://github.com/schrodinger/coordgenlibs/blob/master/LICENSE). The
emphasis of these algorithms are on quality of 2D coordinates rather than speed
of generation. The algorithm distinguishes itself from many others by doing
well with both macrocycles and metal complexes. It also does extremely well on
typical drug-like small molecules, and has been validated on millions of
compounds.

Schrodinger intends to continue to contribute to this code as it still uses it
inside its products, but will also be happy if others contribute pull-requests
when there are improvements they would like to make. We'll also be happy to
hear bug reports or feature requests from use of this code, though make no
guarantee on our ability to process these.

## Documentation

Examples and documentation will be added/improved over time

## Templates

Coordgen uses templates for some macrocycle systems. The source for the
templates is `templates.mae`. If you're an end user of coordgen, you can add
local templates in a file called `user_templates.mae` in a directory specified
by `CoordgenTemplates::setTemplateDir()`. If you want to update the templates,
add new templates to `templates.mae` and run `mol_generator.py` to generate the
source files.

## Usage example

Code for a sample executable is provided in the `example_dir` directory.
Building the example executable is enabled by default, but can be disabled by
means of the `COORDGEN_BUILD_EXAMPLE` option.

## Automated Testing

Automated testing is still primarily taking place inside Schrodinger's internal
build system, although tests are incrementally being added to the `testing`
directory. Building the tests is enabled by default, but can be disabled by
means of the `COORDGEN_BUILD_TESTS` option.

Memory debugging is, by default, configured to use `valgrind`. It can be run on
the tests by passing `-DCMAKE_BUILD_TYPE=Debug` to cmake, to enable building
the debugging symbols, and then using `ctest -T memcheck` inside the build
directory.

## Building from source

### Requirements

To build coordgen, you will need to have the following installed in your
system:

- **CMake** version 3.2 or later.
- The development files for the **Boost libraries**. At least the **iostreams**
  and **regex** components are required. In case of also building the unit
  tests, the **filesystems** and **unit_test_framework** components will also
  be required.
- A **C++ compiler** supporting the C++11 standard.
- A compiled instance of the **maeparser library** or its source code.

In case **maeparser** is not available on your system, neither as a compiled
library or as source code, if a working `git` executable and an internet
connection are available, the builder can automatically download the source and
build **maeparser** for you.

### Building

1. Create a build directory inside the the one that contains Coordgen, and move
   into it:

	```bash
	mkdir build
	cd build
	```

1. Run `cmake` to configure the build, passing the path to the directory where
   the sources are located (just `..` if you created `build` inside the sources
   directory). At this point, you should add any required flags to the `cmake`
   command. Check the 'Options' section in CMakeLists.txt to see which options
   are available.

	```bash
	cmake .. -Dmaeparser_DIR=/home/schrodinger/maeparser_install -DCMAKE_INSTALL_PREFIX=/home/schrodinger/coordgen_install`
	```

	A few notes on the maeparser dependency:

	- CMake will, by default, search your system's default library paths for
	  the maeparser library. If a `CMAKE_INSTALL_PREFIX` was specified to
	  Coordgen, CMake will also search for maeparser there.

	- If you already built and installed maeparser using the
	  `CMAKE_INSTALL_PREFIX` to set the installation path, you should pass the
	  exact same path to Coordgen with `maeparser_DIR`.

	- If CMake cannot find a compiled library for maeparser, it will attempt to
	  download the source code from GitHub and build it. The release to be
	  downloaded if the library is not found can be set using the
	  `-DMAEPARSER_VERSION` flag. The sources will be stored in a directory
	  named like `maeparser-{MAEPARSER_VERSION}` under the coordgen sources.

	- If `maeparser_DIR` was passed to CMake, and the library was not found,
	  CMake will **NOT** download the sources from GitHub (since we expected to
	  find a compiled library).

	- If a copy of maeparser's source is found under the proper path, it be
	  used, instead of being downloaded again.

	- If you want to use Coordgen in a CMake project that also depends on
	  maeparser, set up the maeparser first, as Coordgen will be able to find
	  and use it, without searching for further libraries or compiling it again
	  from the source code.

1. Build and install:

	```bash
	make -j install
	```
