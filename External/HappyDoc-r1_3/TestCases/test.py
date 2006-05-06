#!/usr/bin/python
#
# Here is the module documentation
#
"""Simple test file for HappyDoc.

  This module contains various test cases for testing the HappyDoc
  documentation generator.

  """

import string
import CommandLineApp
import os, types
import one.two
from sys import path
from token import *
from a.b import c, d
from CommandLineApp import TestApp
    

class DefaultClassInst:
    "Class used for default parameter with a class instance."
    pass

class DefaultClassInstWithParams:
    "Class used for default parameter with a class instance taking parameters."
    
    def __init__(self, *args):
        "Initialize DefaultClassInstWithParams instance."
        pass

class DottedBaseClass(A.B):
    "Class to test subclassing from a base class with dots in the name."
    pass

class MultipleBaseClasses(DefaultClassInst, DefaultClassInstWithParams,
                          CommandLineApp.CommandLineApp):
    "Class testing multiple inheritence."
    pass

# Here is a module level variable definition.
foo=1

def example_function_with_args(arg1, arg2,
                               arg3withDefault='hi there',
                               arg3aWithDefault="'hi again'",
                               arg3bWithDefault='"hi there again"',
                               arg4DefaultInt=101,
                               arg5DefaultTuple=(1,2),
                               arg6DefaultList=[3,4],
                               arg7DefaultNone=None,
                               arg8DefaultName=foo,
                               arg9DefaultInstance=DefaultClassInst(),
                               arg10DefaultInstanceWithParams= \
                               DefaultClassInstWithParams(1, 2,
                                                          ('tuple', 'param'),
                                                          ['list', 'param']
                                                          ),
                               stringArgWithHTML='<h1>Hi, Dick & Jane!</h1>',
                               ):
    "This is an example function for testing purposes."
    if one:
        raise IOError('RAISE_class')
    else:
        raise 'RAISE_blah2'
    for i in range(1, 10):
        raise 'RAISE_loop'
    raise 'RAISE_main_level'
    return None

def example_function_without_args():
    "This example function has no arguments."
    pass

def example_function_with_varargs(*args):
    "This example function takes a variable number of arguments."
    pass

def example_function_with_kwargs(**kw):
    """
    This example function takes variable keyword arguments.
    """
    pass

class BaseClass1:
    "First base class."

    #
    # classMember definition example
    #
    classMember=1

    #
    # Documentation for method_of_bc1
    #
    def method_of_bc1(self):
        pass

    def method2_of_bc1(self):
        #
        # Documentation for method2_of_bc1 after name
        #
        pass
    

class SubClass1(BaseClass1):
    "First subclass of BaseClass1."
    pass

class BaseClass2:
    "Second base class."
    pass

class SubClass2(BaseClass1, BaseClass2):
    "Second subclass of BaseClass1 and BaseClass2."
    
    def anotherMethod(self):
        "A method not defined in the base classes."
        pass
    
#
# Class documentation for SubClass3.
#
#   This class documentation section includes a lot of text,
#   in several paragraphs.
#
#   See, here is another paragraph.
#
#   We even have an example:
#
#      a -> b -> c -> a
#
#   And a bulleted list:
#
#      * line one
#
#      * line two
#
#      * line three
#
class SubClass3(SubClass2, MultipleBaseClasses):
    pass

#
# Random comment inserted in text
#

class SubClass4(SubClass1, SubClass3):
    #
    # Class documentation for SubClass4 after the name
    #
    pass



#
# External docs for five are skipped because
# of the blank line following the comment block
#

class five:
    #
    # internal docs for five
    #
    pass


class OverRecursion:
    pass

class OverRecursion(OverRecursion):
    pass

