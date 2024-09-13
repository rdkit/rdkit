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
- Ideas for improvements 
- Write blog posts (either your own or for the [RDKit blog](https://greglandrum.github.io/rdkit-blog/))

#### GitHub Issues vs Discussions


### If you are a developer

It’s generally a good idea to check with Greg or one of the maintainers before adding significant new functionality.

● Contribute interesting scripts/libraries for the Contrib folder 


### If you have 5 minutes to spare 
cleaning up issues
answering questions

## Contributing to the RDKit Docs 

## Contributing to the Code - Python 

This guide already assumes that you are familiar with [the rdkit basics](https://www.rdkit.org/docs/GettingStartedInPython.html), version control with [Git](https://git-scm.com/) and GitHub features such as commits and pull requests. If not, there is a very [neat presentation from UGM 2016](https://github.com/rdkit/UGM_2016/blob/master/Presentations/Landrum_Schneider_GitHub_Git_and_RDKit.pdf) that you should definitely look at before you start.



### Setting up your Environment

Wim summarizes and links whatever is here: https://github.com/rdkit/rdkit/issues/3052

Instructions how to setup IDE here: https://github.com/rdkit/rdkit/issues/3052 . Extended instruction on how to set up your environment for Python development are here: https://greglandrum.github.io/rdkit-blog/posts/2020-03-30-setting-up-an-environment.html



### Coding Style

RDKit does not at the moment has a dedicated style guide or linter setup that you should use for Python. However, there are a few rules that you might want to adhere to so that the code and API design stays coherent (sidenote: it is currently not):

- Use `CamelCase` for `PackageNames`, `ModuleNames`, `ClassNames`, `MethodNames`, `FunctionNames`, but note the different style for `functionAttributeNames` and `methodAttributeNames`
 - These recommendations are made based on the most commonly found patterns in the `Chem` module. However, you may see small inconsistencies here and then...
   - ![style](style.png "Abnormal Style")
- Use `snake_case` names for local variables inside methods and functions. However, `justlowercase` is also fine.
- Use multi-line Python strings to document modules and methods, i.e.:

 ```python
 def xy():
    """ This is a documentation string.
    
    It is continued here...
    """
    pass
 ```
- Use double quotes (`"`) for string literals.

## Contributing to the Code - C++

## Miscellaneous Contributions

### Java

### JavaScript

### Tutorials

### KNIME nodes
