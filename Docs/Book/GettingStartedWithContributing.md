# Getting Started with Contributing to the RDKit

We are really happy to have new contributors help make the RDKit better. This document aims to provide an in-depth overview of what and how to make useful contributions to the RDKit community. 

## What and How to Contribute

### If you are not a developer

Really helpful things you can do that don't involve serious coding:

- Being an active part of the community and asking/answering questions on Github
- Submitting high-quality bug reports via [GitHub issues](https://github.com/rdkit/rdkit/issues)
- Contributing documentation
    - Fixing pieces you find to be unclear or wrong
    - Adding references where available 
    - Adding new sections or examples
    - In-code documentation: comments explaining what the code does or references
    - This Getting Started guide! 
- Contributing tutorials
- Ideas for improvements (with the appropriate tags!) 
- Write blog posts (either your own or for the [RDKit blog](https://greglandrum.github.io/rdkit-blog/))

#### GitHub Issues vs Discussions

### If you are a developer

Submitting new features and or bug fixes is certainly nice. For ideas on where to get started have a look at the [GitHub issues](https://github.com/rdkit/rdkit/issues). We try to remember to tag issues that are relatively small and self-contained as [good for beginners](https://github.com/rdkit/rdkit/labels/good%20for%20beginners); take a look at that list and see if there's anything that looks interesting. If you have questions, ask by adding a comment to the issue.

It’s generally a good idea to check with Greg or one of the maintainers before adding significant new functionality.

Code contributions can be made as follows:
- Core C++ functionality + relevant wrappers
- Python scripts
- JavaScript
- Java
- KNIME nodes

More details on what code contributions (including documentation) should look like are given in the following sections. 

### If you have 5 minutes to spare 
- Cleaning up GitHub issues that have been resolved or that should be under discussions
- Answering questions on the [GitHub Discussions](https://github.com/rdkit/rdkit/discussions) board or [Stack Overflow](https://stackoverflow.com/questions/tagged/rdkit)

### How to Submit Code/Docs Contributions

Contributions are made to the RDKit Codebase via GitHub pull requests. A summary of how to do this is given here. For a more in-depth overview see [this fantastic set of slides](https://github.com/rdkit/UGM_2016/blob/master/Presentations/Landrum_Schneider_GitHub_Git_and_RDKit.pdf) from Greg and Nadine. 



## Contributing to the RDKit Docs 

## Contributing to the Code - Python 

This guide already assumes that you are familiar with [the rdkit basics](https://www.rdkit.org/docs/GettingStartedInPython.html), version control with [Git](https://git-scm.com/) and GitHub features such as commits and pull requests. If not, there is a very [neat presentation from UGM 2016](https://github.com/rdkit/UGM_2016/blob/master/Presentations/Landrum_Schneider_GitHub_Git_and_RDKit.pdf) that you should definitely look at before you start.



### Setting up your Environment

Wim summarizes and links whatever is here: https://github.com/rdkit/rdkit/issues/3052

Instructions how to setup IDE are here: https://github.com/rdkit/rdkit/issues/3052 . Extended instruction on how to set up your environment for Python development are here in a post on the [RDKit blog](https://greglandrum.github.io/rdkit-blog/posts/2020-03-30-setting-up-an-environment.html). At the bottom, it includes a recipe for cloning a local fork on which you can then run and test the updated code.



### Coding Style

RDKit does not at the moment has a dedicated style guide or a linter that you should use. However, there are a few rules that you might want to adhere to so that the code and API design stays coherent:

- use **camel case** for **ModuleNames**, "ClassNames", "MethodNames"

## Contributing to the Code - C++

### Adding a new unit test

New unit tests should use the Catch2 test framework ([https://github.com/catchorg/Catch2](https://github.com/catchorg/Catch2)). Test files are called “catch\_something.cpp”, and are placed directly in the code folders. If you contribute a small change / bug-fix, you don't need a new file, but just add a TEST\_CASE. Example:

```c++
TEST_CASE("test basic Mol features", "[ROMol]") {
  SECTION("basics") {
    {
      // shortcut to create a molecule from SMILES
      auto m = "CCO"_smiles;
      // check that the molecule was created successfully
      REQUIRE(m);
      // test some functionality
      int nAtoms = m->getNumAtoms();
      CHECK(nAtoms == 3);
    }
  }
}
```

### Adding a new test file

For larger features, you can create a separate test file “catch\_myfeature.cpp” in the respective folder:

```c++
#include "RDGeneral/test.h"
#include <catch2 /catch_all.hpp>
#include <GraphMol /SmilesParse/SmilesParse.h>

using namespace RDKit;

TEST_CASE("test basic Mol features", "[ROMol]") {
  SECTION("basics") {
    {
    // shortcut to create a molecule from SMILES
    auto m = "CCO"_smiles;
    // check that the molecule was created successfully
    REQUIRE(m);
    // test some functionality
    int nAtoms = m->getNumAtoms();
    CHECK(nAtoms == 3);
    }
  }
}
```

Additionally, adapt the CMakeLists.txt of the respective folder. You will have to adapt filenames and add needed libraries (in this example, FileParsers contains the functionality to parse SMILES.

```
rdkit_catch_test(dummyTest catch_dummy.cpp LINK_LIBRARIES FileParsers)
```

## Miscellaneous Contributions

### Java

### JavaScript

### Tutorials

### KNIME nodes
