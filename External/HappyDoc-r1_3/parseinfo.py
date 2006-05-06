#!/usr/bin/env python
#
# Time-stamp: <01/01/21 20:02:20 dhellmann>
#
# COPYRIGHT
#
#   Permission to use, copy, modify, and distribute this software and
#   its documentation for any purpose and without fee is hereby
#   granted, provided that the above copyright notice appear in all
#   copies and that both that copyright notice and this permission
#   notice appear in supporting documentation, and that the name of
#   Doug Hellmann not be used in advertising or publicity pertaining
#   to distribution of the software without specific, written prior
#   permission.
#
# DISCLAIMER
# 
#   DOUG HELLMANN DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
#   SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
#   FITNESS, IN NO EVENT SHALL DOUG HELLMANN BE LIABLE FOR ANY
#   SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#   WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN
#   AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
#   ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
#   THIS SOFTWARE.
# 

"""Extract information from a Python code parse tree.

  This module is based on the Demos/parser/example.py module
  distributed with the Python source distribution.

Original comments

  Simple code to extract class & function docstrings from a module.

  This code is used as an example in the library reference manual in
  the section on using the parser module.  Refer to the manual for a
  thorough discussion of the operation of this code.

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: parseinfo.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sun, 12-Mar-2000 11:20:10 EST',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.2 $',
    'date':'$Date: 2001/03/21 22:31:06 $',
    'locker':'$Locker:  $',
    }

#
# Import system modules
#
import parser
import symbol
import token
import types
import pprint
import sys
import string
import re
import os
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

#
# Import Local modules
#
import StructuredText
from hdpathmanagement import applyPrefixToPath

#
# Module
#

def extractComments(
    text,
    extractRe=re.compile('^((?P<blankline>\s*$)|(?P<namedobj>\s*(?P<nametype>(class|def))\s+(?P<name>[0-9A-Za-z_]+))|(?P<commentline>\s*#+(?P<comment>.*)))').search,
    ):
    """Given a block of Python source, extract the comments.

    The comment text is associated with nearby named objects
    (functions, methods, classes, etc.).  This function returns
    a dictionary of names and the associated comment text.

    """
    dbg=None

    comment_info = {}

    comment_text = ''
    current_name = None

    f = StringIO(text)
    line = f.readline()
    while line:
        if dbg and (dbg >= 2):
            print '=>%s' % string.rstrip(line)
        match_obj = extractRe(line)
        if match_obj:
            #
            # Documentation before named object
            #
            match_dict = match_obj.groupdict()
            comment = match_dict['comment']
            name = match_dict['name']
            blankline = match_dict['blankline']

            if match_dict['commentline'] and not comment:
                comment = ' '

            if comment:
                # Append new text to running buffer.
                if dbg:
                    print 'PARSEINFO: Adding comment text.'
                comment_text = '%s%s\n' % (comment_text, comment,)

            elif name and comment_text:
                if current_name:
                    # Hit a new name, store the comment_text buffer
                    # for the current_name and switch current_name.
                    if dbg:
                        print 'PARSEINFO: Storing comment for %s, switching to %s' % \
                              (current_name, name)
                    comment_info[current_name] = comment_text
                    comment_text = ''
                    current_name = name
                else:
                    # Hit a new name with existing comment_text,
                    # store the comment along with that name.
                    if dbg:
                        print 'PARSEINFO: Storing comment for %s' % name
                    comment_info[name] = comment_text
                    comment_text = ''
                    current_name = None

            elif name:
                # Recognized new name definition.
                if dbg:
                    print 'PARSEINFO: New name %s' % name
                current_name = name

            elif blankline:
                # Reset when a blank line separates comment from
                # named stuff.
                if dbg:
                    print 'PARSEINFO: blank line'
                if comment_text and current_name:
                    if not comment_info.get(current_name, None):
                        if dbg:
                            print 'PARSEINFO: Storing comment after name %s' \
                                  % current_name
                        comment_info[current_name] = comment_text
                else:
                    if dbg:
                        if comment_text:
                            print 'PARSEINFO: Discarding comment "%s"' % comment_text
                current_name = None
                comment_text = ''

            elif current_name and comment_text:
                # Store all comment text for the current_name.
                if dbg:
                    print 'PARSEINFO: Storing comment for %s' % current_name
                comment_info[current_name] = comment_text
                comment_text = ''
                current_name = None

        else:
            if dbg:
                print 'PARSEINFO: Not matched (%s)' % string.strip(line)
            current_name = None
            comment_text = ''

        line = f.readline()

    f.close()

    if current_name and comment_text:
        # Final storage to make sure we have everything.
        if dbg:
            print 'PARSEINFO: Final storage of comment for %s' % current_name
        comment_info[current_name] = comment_text

    return comment_info


def getDocs(fileName, includeComments=1):
    """Retrieve information from the parse tree of a source file.

    Parameters
    
      fileName --
        Name of the file to read Python source code from.

      includeComments=1 --
        Flag to indicate whether comments should be parsed for
        cases where __doc__ strings are not available.
        
    """
    f = open(fileName,'rU')
    #
    # Read file and add an extra newline to fix problem
    # reported with files containing only a single docstring
    # line.
    #
    source = '%s\n' % f.read()
    f.close()
    basename = os.path.basename(os.path.splitext(fileName)[0])
    ast = parser.suite(source)
    tup = parser.ast2tuple(ast)
    if includeComments:
        comment_info = extractComments(source)
    else:
        comment_info = {}
    mod_info = ModuleInfo(tree=tup,
                          name=basename,
                          fileName=fileName,
                          commentInfo=comment_info)
    return mod_info




class SuiteInfoBase:
    """Base class for information gathering classes.

    Default implementation assumes that the user is interested
    in learning about functions and classes defined within the
    parse tree used for initialization.  This makes implementation
    of MethodInfo easy.  Other derived classes add behavior to
    find other information.
    """
    
    #_docstring = ''
    _docstring_summary = None
    _name = ''

    def __init__(self, name = '', parent = None, tree = None, commentInfo={}):
        """Initialize the info extractor.

        Parameters:

            name -- name of this object
        
            parent -- parent object (e.g. Module for Class)
            
            tree -- parse tree from which to extract information

            commentInfo -- comments extracted from source file where
            this object was found
            
        """
        self._class_info = {}
        self._function_info = {}
        self._namespaces = ( self._class_info, self._function_info )
        self._name = name
        self._parent = parent
        self._comment_info = commentInfo
        self._comments = commentInfo.get(self._name, '')
        if tree:
            self._extractInfo(tree)
        return

    def __getitem__(self, itemName):
        for ns in self._namespaces:
            if ns.has_key(itemName):
                return ns[itemName]
        raise KeyError('Unrecognized name: "%s"' % itemName, itemName)

    _summary_pattern = re.compile(r'^\s*([^\n]+)\n')
    def _extractSummary(self, text):
        "Extract a summary text from a larger body."
        text = string.strip(text)
        #
        # Remove surrounding quotes, if present.
        #
        while text and (text[0] in ('"', "'")):
            text = text[1:]
        while text and (text[-1] in ('"', "'")):
            text = text[:-1]
        #
        # Pull out the first line, and return it if
        # we can find it.  Otherwise, return the whole
        # string since that means that the whole thing
        # is just one line.
        #
        matchObj = self._summary_pattern.search(text)
        if matchObj:
            return string.strip(matchObj.group(0))
        else:
            return text
    
    def _extractInfo(self, tree):
        "Pull information out of the parse tree."
        # extract docstring
        if len(tree) == 2:
            found, vars = match(DOCSTRING_STMT_PATTERN[1], tree[1])
        elif len(tree)==4:
            found, vars = match(DOCSTRING_STMT_PATTERN, tree[3])
        else:
            found = False

        if found:
            self._docstring = eval(vars['docstring'])
        else:
            self._docstring = ''
        # discover inner definitions
        for node in tree[1:]:
            found, vars = match(COMPOUND_STMT_PATTERN, node)
            if found:
                cstmt = vars['compound']
                if cstmt[0] == symbol.funcdef:
                    name = cstmt[2][1]
                    self._function_info[name] = FunctionInfo(
                        tree=cstmt,
                        parent=self,
                        commentInfo=self._comment_info,
                        )
                    #pprint.pprint(cstmt)
                elif cstmt[0] == symbol.classdef:
                    name = cstmt[2][1]
                    self._class_info[name] = ClassInfo(
                        tree=cstmt,
                        parent=self,
                        commentInfo=self._comment_info
                        )
        return
    
    def getDocString(self):
        "Return any __doc__ string value found for the object."
        dstring = '%s\n\n%s' % (self._docstring, self._comments)
        #print 'DOC STRING for %s is ' % self._name, dstring
        return dstring

    def getSummary(self):
        "Return a summary of the __doc__ string for this object."
        if self._docstring_summary is None:
            self._docstring_summary = \
                                    self._extractSummary(self.getDocString())
        return self._docstring_summary

    def getName(self):
        "Return the name of the object."
        return self._name

    def __str__(self):
        return '%s(%s) ' % (repr(self), self.getName())

    def getFullyQualifiedName(self, filenamePrefix=''):
        "Return a complete, unique, name representing this object."
        #print 'BBB getFullyQualifiedName(%s)' % self
        if not self._parent:
            #print 'BBB (no parent)'
            name = self.getQualifiedName(filenamePrefix)
        else:
            #print 'BBB (with parent)'
            parent_name = self._parent.getFullyQualifiedName(filenamePrefix)
            parent_base, parent_ext = os.path.splitext( parent_name )
            name = '%s_%s%s' % ( parent_base, self.getName(), parent_ext )
        #print 'BBB =>%s' % name
        return name

    def getQualifiedName(self, filenamePrefix='', transTable=string.maketrans('/', '_')):
        #print 'YYY getQualifiedName(%s)' % self.getName()
        if not self._parent:
            #
            # Start with the filename for this object
            #
            if filenamePrefix:
                name = applyPrefixToPath(self._filename, filenamePrefix)
            else:
                name = self._filename
            #
            # Remove preceding slashes to make name relative
            #
            while name[0] in './':
                name = name[1:]
        else:
            name = '%s_%s' \
                   % (os.path.basename(self._parent.getQualifiedName(filenamePrefix)),
                      self.getName())
        #print 'YYY\t=>%s' % name
        return name
        
        
    def getClassNames(self):
        "Return the names of classes defined within the module."
        return self._class_info.keys()

    def getClassInfo(self, name):
        "Return a ClassInfo object for the class by name."
        return self._class_info[name]




        
class SuiteFuncInfo:
    #  Mixin class providing access to function names and info.

    def getFunctionNames(self):
        return self._function_info.keys()

    def getFunctionInfo(self, name):
        return self._function_info[name]





    

    
    
from types import ListType, TupleType

def joinCodeSnippets(first, second, separator):
    """Join two code snippets into one string.

    Use some general code content rules to try to make the
    resulting snippet look nice.
    """
    
    if string.strip(second) in ('.',):
        sep_to_be_used = ''
        
    elif second and ( second[0] in ('.', ',', '(',) ):
        if second[0] == '(' and first and first[-1] == ',':
            sep_to_be_used = separator
        else:
            sep_to_be_used = ''
        
    elif (not first) or (first and first[-1] in ('.',)) or (first in ('-',)):
        sep_to_be_used = ''
        
    elif ( (first and ( first[-1] in ('(', '[', '{') )) and
           (second and ( second[-1] in (')', ']', '}') ))
           ):
        sep_to_be_used = ''
        
    elif first and ( first[-1] in ('(', '[', '{') ):
        sep_to_be_used = separator
        
    else:
        sep_to_be_used = separator
                
    text = '%s%s%s' % (first, sep_to_be_used, second)
    return text



def parseTreeToString(tree, separator=' '):
    """Convert a parse tree to a string which would have parsed in that way.

    Given a parse tree, walk it to determine the original string
    which would have been parsed to produce that tree.
    """
    #pprint 'STRINGING: ',
    #pprint.pprint(tree)
    text = ''
    if tree and type(tree) in (types.TupleType, types.ListType):
        if type(tree[0]) in (types.TupleType, types.ListType):
            tree_parts = tree
        else:
            tree_parts = tree[1:]
        sub_parts = map( lambda x, s=separator: parseTreeToString(x, s),
                         tree_parts)
        for one_part in sub_parts:
            text = joinCodeSnippets(text, one_part, separator)
    else:
        text = str(tree)
    return text




def match(pattern, data, vars=None, dbg=0):
    """Match `data' to `pattern', with variable extraction.

    pattern --
        Pattern to match against, possibly containing variables.

    data --
        Data to be checked and against which variables are extracted.

    vars --
        Dictionary of variables which have already been found.  If not
        provided, an empty dictionary is created.

    The `pattern' value may contain variables of the form ['varname'] which
    are allowed to match anything.  The value that is matched is returned as
    part of a dictionary which maps 'varname' to the matched value.  'varname'
    is not required to be a string object, but using strings makes patterns
    and the code which uses them more readable.

    This function returns two values: a boolean indicating whether a match
    was found and a dictionary mapping variable names to their associated
    values.
    """
    if dbg:
        print 'PATTERN: ',
        pprint.pprint(pattern)
        print 'DATA:',
        pprint.pprint(data)
    if vars is None:
        vars = {}
    if type(pattern) is ListType:       # 'variables' are ['varname']
        if dbg:
            print 'storing "%s" for variable "%s"' % (data, pattern[0])
        vars[pattern[0]] = data
        return 1, vars
    if type(pattern) is not TupleType:
        if dbg:
            print 'end recursion'
        return (pattern == data), vars
    if len(data) != len(pattern):
        if dbg:
            print 'shortcut, length does not match'
        return 0, vars
    for pattern, data in map(None, pattern, data):
        if dbg:
            print 'recursing'
        same, vars = match(pattern, data, vars, dbg=dbg)
        if not same:
            break
    return same, vars





def lenientMatch(pattern, data, vars=None, dbg=0):
    """Match `data' to `pattern', with variable extraction.

    pattern --
        Pattern to match against, possibly containing variables.

    data --
        Data to be checked and against which variables are extracted.

    vars --
        Dictionary of variables which have already been found.  If not
        provided, an empty dictionary is created.

    The `pattern' value may contain variables of the form ['varname'] which
    are allowed to match anything.  The value that is matched is returned as
    part of a dictionary which maps 'varname' to the matched value.  'varname'
    is not required to be a string object, but using strings makes patterns
    and the code which uses them more readable.

    This function is based on the match() function, but is more lenient.
    The pattern does not have to completely describe the tree.  Instead,
    it can be the 'top' portion of the tree.  Everything must match down
    to the leaves of the pattern.  At that point, the matching stops.  If
    a match was found at all, the return values indicate a match.
    
    This function returns two values: a boolean indicating whether a match
    was found and a dictionary mapping variable names to their associated
    values.
    """
    if dbg:
        print 'PATTERN : ',
        pprint.pprint(pattern)
        print 'DATA    :',
        #pprint.pprint(data)
        print data
    if vars is None:
        vars = {}
    if type(pattern) is ListType:       # 'variables' are ['varname']
        if dbg:
            print 'storing "%s" for variable "%s"' % (data, pattern[0])
        vars[pattern[0]] = data
        return 1, vars
    if type(pattern) is not TupleType:
        if dbg:
            print 'end recursion'
        return (pattern == data), vars
    found_match = 0
    if pattern and data:
        for pattern, data in map(None, pattern, data):
            if dbg:
                print 'recursing'
            same, vars = lenientMatch(pattern, data, vars, dbg=dbg)
            if not same:
                break
            else:
                found_match = same
    return found_match, vars




#  This pattern identifies compound statements, allowing them to be readily
#  differentiated from simple statements.
#
COMPOUND_STMT_PATTERN = (
    symbol.stmt,
    (symbol.compound_stmt, ['compound'])
    )



# This pattern matches the name of an item which appears in a
# testlist sequence.  This can be used for finding the
# base classes of a class, or the parameters to a function or
# method.
#
# BASE_CLASS_NAME_PATTERN = (
#     symbol.test,
#     (symbol.and_test,
#      (symbol.not_test,
#       (symbol.comparison,
#        (symbol.expr,
#         (symbol.xor_expr,
#          (symbol.and_expr,
#           (symbol.shift_expr,
#            (symbol.arith_expr,
#             (symbol.term,
#              (symbol.factor,
#               (symbol.power,
#                (symbol.atom, 
#                 (token.NAME, ['name'])
#                 )))))))))))))
BASE_CLASS_NAME_PATTERN = (
    symbol.test,
    (symbol.and_test,
     (symbol.not_test,
      (symbol.comparison,
       (symbol.expr,
        (symbol.xor_expr,
         (symbol.and_expr,
          (symbol.shift_expr,
           (symbol.arith_expr,
            (symbol.term,
             (symbol.factor, ['power'])
             ))))))))))



#  This pattern will match a 'stmt' node which *might* represent a docstring;
#  docstrings require that the statement which provides the docstring be the
#  first statement in the class or function, which this pattern does not check.
#
DOCSTRING_STMT_PATTERN = (
    symbol.stmt,
    (symbol.simple_stmt,
     (symbol.small_stmt,
      (symbol.expr_stmt,
       (symbol.testlist,
        (symbol.test,
         (symbol.and_test,
          (symbol.not_test,
           (symbol.comparison,
            (symbol.expr,
             (symbol.xor_expr,
              (symbol.and_expr,
               (symbol.shift_expr,
                (symbol.arith_expr,
                 (symbol.term,
                  (symbol.factor,
                   (symbol.power,
                    (symbol.atom,
                     (token.STRING, ['docstring'])
                     )))))))))))))))),
     (token.NEWLINE, '')
     ))



class ImportInfo:
    """Collects info about imports for a module.
    """

    def __init__(self):
        self._straight_imports = []
        self._named_imports = {}
        return

    def addImport(self, moduleName, symbolName=None):
        """Add information about an import statement to the saved info.

        Parameters

          moduleName -- The name of the module involved in the import.
          For example, in 'from X import Y', X is the moduleName and
          in 'import A.B', A.B is the moduleName.

          symbolName -- The name of the symbol being imported.  For
          example, in 'from X import Y', Y is the symbolName.
        """
        dbg=0
        
        if symbolName:
            if dbg: print '\nIMPORT SYMBOL %s from MODULE %s' % (symbolName, moduleName)
            name_list = self._named_imports.get(moduleName, [])
            if symbolName not in name_list:
                if dbg: print '\t*** added'
                name_list.append(symbolName)
                self._named_imports[moduleName] = name_list
                
        else:
            if dbg: print '\nIMPORT MODULE: %s' % moduleName
            if moduleName not in self._straight_imports:
                if dbg: print '\t*** added'
                self._straight_imports.append(moduleName)
        if dbg:
            print 'STRAIGHT: ',
            pprint.pprint(self._straight_imports)
            print 'NAMED: ',
            pprint.pprint(self._named_imports)
            print 'CURRENT IMPORTS: ', self.items()
        return
    
    def importedSymbols(self, moduleName):
        if self._named_imports.has_key(moduleName):
            return self._named_imports[moduleName]
        else:
            raise ValueError('No symbols imported for module', moduleName)
        return

    def __str__(self):
        return '(%s)' % string.join( map(str, self.items()),
                                     '\n'
                                     )
    
    def items(self):
        """Returns a sequence of tuples containing module names and the
        symbols imported from them.
        """
        all_names = self._straight_imports[:]
        for name in self._named_imports.keys():
            if name not in all_names:
                all_names.append(name)
        all_names.sort()
        all_items = []
        for name in all_names:
            if name in self._straight_imports:
                all_items.append( (name, None) )
            if self._named_imports.has_key(name):
                all_items.append( (name, self._named_imports[name]) )
        return all_items


class ModuleInfo(SuiteInfoBase, SuiteFuncInfo):
    """Information gatherer for source code modules.

    Extract information about a source module from
    its parse tree.
    """
    
    def __init__(self, parent = None,
                 tree = None,
                 name = "<string>",
                 fileName = None,
                 commentInfo = {}
                 ):
        """Initialize the info extractor.

        Parameters:

            tree -- parse tree from which to extract information

            name -- name of the module

            fileName -- name of the file containing the module

            commentInfo -- comments extracted from the file

        """
        if os.path.basename(fileName) == '__init__.py':
            self._filename = os.path.dirname(fileName)
        else:
            self._filename = fileName
        SuiteInfoBase.__init__(self, name=name, parent=parent, tree=tree,
                               commentInfo=commentInfo)
        if tree:
            #
            # Look for doc string
            #
            found, vars = match(DOCSTRING_STMT_PATTERN, tree[1])
            if found:
                self._docstring = vars["docstring"]
            #
            # Look for imported modules
            #
            self._import_info = self._extractImportedModules(tree)
        return

    def getFileName(self):
        "Returns the name of the file containing the module."
        return self._filename
    
    def getImportData(self):
        "Returns a list of which symbols are imported."
        return self._import_info.items()

    def _extractImportedModules(self, tree):
        """Returns modules imported by code in tree.
        
        Scan the parse tree for import statements
        and return the names of all modules imported.
        """
        dbg=0
        IMPORT_STMT_WITH_LIST_PATTERN =(
            symbol.stmt,
            (symbol.simple_stmt,
             (symbol.small_stmt,
              ['import_stmt']
              ),
             (token.NEWLINE, '')
             )
            )
        imported_modules = ImportInfo()
        for subtree in tree[1:]:
            if dbg: print '\nNEW IMPORT SUBTREE'
            found, vars = match(IMPORT_STMT_WITH_LIST_PATTERN, subtree)
            if dbg: print 'vars: ',
            if dbg: pprint.pprint(vars)
            if found:
                # vars['import_stmt'] should be an import statement
                # in one of several different forms
                import_stmt = vars['import_stmt']
                if import_stmt[0] == symbol.import_stmt:
                    first = import_stmt[1]
                    if (first[0] == token.NAME) and (first[1] == 'import'):
                        for import_module in import_stmt[2:]:
                            if dbg: print 'import_module: ',
                            if dbg: pprint.pprint(import_module)

                            if import_module[0] == symbol.dotted_name:
                                # Get the tuples with the name
                                module_name_parts = import_module[1:]
                                # Get the strings in the 2nd part of each tuple
                                module_name_parts = map(lambda x: x[1],
                                                        module_name_parts)
                                # Combine the strings into the name
                                module_name = string.join(module_name_parts,
                                                          '')
                                if dbg: print 'ADDING module_name=%s' % module_name
                                imported_modules.addImport(module_name)
                            elif import_module[0] == symbol.dotted_as_name:
                                module_name_parts = import_module[1][1:]
                                module_name_parts = map(lambda x: x[1],
                                                        module_name_parts)
                                # Combine the strings into the name
                                module_name = string.join(module_name_parts,'')
                                if dbg: print 'ADDING2 module_name=%s' % module_name
                                imported_modules.addImport(module_name)
                    elif (first[0] == token.NAME) and (first[1] == 'from'):
                        if dbg: print 'from x import y'
                        module_name = parseTreeToString(import_stmt[2])
                        #imported_modules[module_name] = imported_modules.get(
                        #    module_name, [])
                        try:
                            symbol_list = imported_modules.importedSymbols(module_name)
                        except ValueError:
                            symbol_list = []
                        names = import_stmt[4:]
                        if dbg: print 'NAMES: ', names
                        for n in names:
                            if n[0] == token.NAME:
                                #symbol_list.append(n[1])
                                imported_modules.addImport(module_name, n[1])
                            elif n[0] == symbol.import_as_name:
                                if n[1][0] == token.NAME:
                                    if len(n)==2:
                                        imported_modules.addImport(module_name, n[1][1])
                                    else:
                                        name_parts = n[1:]
                                        name_parts = map(lambda x: x[1],name_parts)
                                        extras = string.join(name_parts,' ')
                                        imported_modules.addImport(module_name,extras)
                            elif n[0] == token.STAR:
                                #symbol_list.append('*')
                                imported_modules.addImport(module_name, '*')
                    if dbg:
                        for part in import_stmt[1:]:
                            pprint.pprint(part)
            if dbg: print 'ITERATION IMPORTS: ', imported_modules
        if dbg: print 'FINAL IMPORTS: ', imported_modules, '\n'
        return imported_modules




    
class ClassInfo(SuiteInfoBase):
    "Gather information about a Python class from its parse tree."
    
    def __init__(self, parent = None, tree = None, commentInfo = {}):
        """Initialize the info extractor.

        Parameters:

            parent -- parent object for this class (e.g. Module)

            tree -- parse tree from which to extract information

            commentInfo -- comments extracted from the source file
            for this class

        """
        SuiteInfoBase.__init__(self, name=tree[2][1], parent=parent,
                               tree=(tree and tree[-1] or None),
                               commentInfo=commentInfo)
        self._base_class_info = self._extractBaseClasses(tree)
        #print self._base_class_info
        self._class_member_info = self._extractClassMembers(tree)
        #print self._class_member_info
        return

    def _extractBaseClasses(self, tree):
        "Returns a list of all base classes from which this class is derived."
        #print
        #pprint.pprint(tree)
        base_class_names = []
        for subtree in tree[1:]:
            #pprint.pprint(subtree)
            if subtree[0] == symbol.testlist:
                for test in subtree[1:]:
                    found, vars = lenientMatch(BASE_CLASS_NAME_PATTERN, test)
                    #pprint.pprint(vars)
                    if found and vars.has_key('power'):
                        #base_class_names.append(vars['name'])
                        name = parseTreeToString(vars['power'])
                        base_class_names.append(name)
        return base_class_names

    def _extractClassMembers(self, tree):
        """Returns a list of all variable assignments
        in the class member context."""
        CLASS_MEMBER_STMT_PATTERN = (
            symbol.stmt,
            (symbol.simple_stmt,
             (symbol.small_stmt,
              (symbol.expr_stmt,
               (symbol.testlist,
                (symbol.test,
                 (symbol.and_test,
                  (symbol.not_test,
                   (symbol.comparison,
                    (symbol.expr,
                     (symbol.xor_expr,
                      (symbol.and_expr,
                       (symbol.shift_expr,
                        (symbol.arith_expr,
                         (symbol.term,
                          (symbol.factor,
                           (symbol.power,
                            (symbol.atom,
                             (token.NAME, ['member_name']),
                             ))))))))))))))))))
        #print 'CLASS_MEMBER_STMT_PATTERN: ',
        #pprint.pprint(CLASS_MEMBER_STMT_PATTERN)
        #print
        #print 'TREE: ',
        #pprint.pprint(tree)
        #print
        #
        # Find the suite defining the class
        #
        for subtree in tree[1:]:
            #print 'SUBTREE: ',
            #pprint.pprint(subtree)
            #print
            if subtree[0] == symbol.suite:
                search_in = subtree
                break

        #print 'SEARCH_IN: ',
        #pprint.pprint(search_in)
        #print

        class_members = []
        
        for subtree in search_in[1:]:
            found, vars = lenientMatch(CLASS_MEMBER_STMT_PATTERN,
                                       subtree,
                                       dbg=0)
            if found and vars:
                #print vars
                class_members.append(vars['member_name'])
            #sys.exit(1)

        #print class_members
        return class_members

    def getMethodNames(self):
        "Returns a list of the names of methods defined for this class."
        return self._function_info.keys()

    def getMethodInfo(self, name):
        "Returns a FunctionInfo object for the method 'name', if it exists."
        return self._function_info[name]

    def getBaseClassNames(self):
        "Returns a list of the names of the base classes for this class."
        return self._base_class_info

    def getExceptionNames(self):
        "Returns a list of the names of all exceptions raised by this class."
        exception_names = []
        return exception_names


    
class FunctionInfo(SuiteInfoBase, SuiteFuncInfo):
    "Gather information about a function or method definition."
    
    def __init__(self, parent=None, tree = None, commentInfo={}):
        """Initialize the info extractor.

        Parameters:

            parent -- parent object for this object (e.g. Module or Function)

            tree -- parse tree from which to extract information

            commentInfo -- comments extracted from the source file holding
            this module

        """
        SuiteInfoBase.__init__(self,
                               name = tree[2][1],
                               parent=parent,
                               tree=(tree and tree[-1] or None),
                               commentInfo=commentInfo)
        parameter_data = self._extractFunctionParameters(tree)
        self._constructParameterInfo(parameter_data)
        self._exception_info = self._extractThrownExceptions(tree)
        #if self._exception_info:
        #    print 'EXCEPTIONS: ',
        #    pprint.pprint(self._exception_info)
        return

    ##
    ## EXCEPTIONS
    ##
    
    EXCEPTION_BY_NAME_PATTERN = (
        symbol.raise_stmt,
        (token.NAME,),
        (symbol.test,
         (symbol.and_test,
          (symbol.not_test,
           (symbol.comparison,
            (symbol.expr,
             (symbol.xor_expr,
              (symbol.and_expr,
               (symbol.shift_expr,
                (symbol.arith_expr,
                 (symbol.term,
                  (symbol.factor, ['exception'])
                  ))))))))))
        )
    
    EXCEPTION_STRING_PATTERN = (
        symbol.raise_stmt,
        (token.NAME,),
        (symbol.test,
         (symbol.and_test,
          (symbol.not_test,
           (symbol.comparison,
            (symbol.expr,
             (symbol.xor_expr,
              (symbol.and_expr,
               (symbol.shift_expr,
                (symbol.arith_expr,
                 (symbol.term,
                  (symbol.factor,
                   (symbol.power,
                    (symbol.atom,
                     (token.STRING, ['exception'])
                     )))))))))))))
        )

    def getExceptionNames(self):
        "Return a list of the names of any exceptions raised by the function."
        return self._exception_info.keys()

    def getExceptionInfo(self, exceptionName):
        """Returns a type value for an exception.

        The return value will be one of (token.NAME, token.STRING)
        indicating whether the exception was thrown as a string
        or a named object.
        """
        return self._exception_info[exceptionName]

    def _extractThrownExceptions(self, tree):
        "Return a dictionary of exception->exception_type values."
        thrown_exceptions = {}
        if not tree:
            return thrown_exceptions
        if type(tree) in (types.ListType, types.TupleType):
            if tree[0] == symbol.raise_stmt:
                #print 'found raise: ', parseTreeToString(tree)
                #print 'parsing...'

                #
                # Initialize
                #
                exception_name = None
                exception_type = None
                
                if not exception_name:
                    # Look for throwing a string
                    found, vars = lenientMatch(
                        self.EXCEPTION_STRING_PATTERN,
                        tree,
                        #dbg=1
                        )
                    if found:
                        try:
                            exception_name = vars['exception']
                            exception_type = token.STRING
                            #print 'GOT EXCEPTION: ', exception_name
                        except KeyError:
                            pass

                if not exception_name:
                    # Look for throwing a name.
                    found, vars = lenientMatch(
                        self.EXCEPTION_BY_NAME_PATTERN,
                        tree,
                        #dbg=1
                        )
                    if found:
                        #print 'FOUND EXCEPTION: ', vars
                        # Threw a named thing, record the name.
                        try:
                            exception_name = parseTreeToString(vars['exception'])
                            exception_type = token.NAME
                            #print 'GOT EXCEPTION: ', exception_name
                        except KeyError:
                            pass
                        
                if not exception_name:
                    #print 'NO NAME,',
                    try:
                        slice=tree[2:]
                        if slice:
                            #print 'using slice of 2:1=', slice
                            exception_name = parseTreeToString(slice)
                    except:
                        #print 'using whole tree=', tree
                        execption_name = parseTreeToString(tree)

                if exception_name:
                    #print 'STORING REFERENCE'
                    thrown_exceptions[exception_name] = exception_type
            else:
                for subtree in tree[1:]:
                    thrown_exceptions.update(
                        self._extractThrownExceptions(subtree)
                        )
        #print 'EXCEPTIONS: ', thrown_exceptions.keys()
        return thrown_exceptions

    ##
    ## PARAMETERS
    ##
    
    # This pattern matches the name of an item which appears in a
    # testlist sequence.  This can be used for finding the
    # base classes of a class, or the parameters to a function or
    # method.
    #
    PARAMETER_DEFAULT_PATTERN = (
        symbol.test,
        (symbol.and_test,
         (symbol.not_test,
          (symbol.comparison,
           (symbol.expr,
            (symbol.xor_expr,
             (symbol.and_expr,
              (symbol.shift_expr,
               (symbol.arith_expr, ['term'], ['trailer'], ['trailer_bits'])
               ))))))))
    
    oldPARAMETER_DEFAULT_PATTERN = (
        symbol.test,
        (symbol.and_test,
         (symbol.not_test,
          (symbol.comparison,
           (symbol.expr,
            (symbol.xor_expr,
             (symbol.and_expr,
              (symbol.shift_expr,
               (symbol.arith_expr,
                (symbol.term,
                 (symbol.factor,
                  (symbol.power, ['atom'])
                    )))))))))))
    PARAMETER_DEFAULT_WITH_TRAILER_PATTERN = (
        symbol.test,
        (symbol.and_test,
         (symbol.not_test,
          (symbol.comparison,
           (symbol.expr,
            (symbol.xor_expr,
             (symbol.and_expr,
              (symbol.shift_expr,
               (symbol.arith_expr,
                (symbol.term,
                 (symbol.factor,
                  (symbol.power, ['atom'], ['trailer'], ['trailer_bits'])
                  ))))))))))
        )
    PARAMETER_ARITH_DEFAULT_WITH_TRAILER_PATTERN = (
        symbol.test,
        (symbol.and_test,
         (symbol.not_test,
          (symbol.comparison,
           (symbol.expr,
            (symbol.xor_expr,
             (symbol.and_expr,
              (symbol.shift_expr,
               (symbol.arith_expr,
                #(symbol.term,
                # (symbol.factor,
                #  (symbol.power, ['atom'], ['trailer'], ['trailer_bits'])
                #  ))
                ['expression'], ['trailer'], ['trailer_bits']
                ))))))))
        )

    def _constructParameterInfo(self, parameterData):
        """Construct storable parameter data from a parameter list.
        
        Given the sequence of tuples extracted as a parameter list,
        store the names (in order) in self._parameter_names and the
        information about the parameter in self._parameter_info where
        the keys are the parameter name and the info is a tuple containing:

        (default_specified, default_value, default_value_type)

        Where:

            default_specified -- boolean indicating whether a default value
                                 was specified

            default_value -- the default value given, if any

            default_value_type -- the type of the default value (token.STRING,
                                  token.NAME, None). A type of None means
                                  unknown.
            
        """
        parameter_info = {}
        parameter_names = []
        for (param, default_specified,
             default_value, default_value_type) in parameterData:
            parameter_names.append(param)
            parameter_info[param] = ( default_specified,
                                      default_value,
                                      default_value_type,
                                      )
        self._parameter_names = tuple(parameter_names)
        self._parameter_info = parameter_info
        return

    def getParameterNames(self):
        """Returns a list of the names of all
        parameters to the function, in order."""
        return self._parameter_names

    def getParameterInfo(self, paramName):
        """Returns the info record for a parameter.

        The returned tuple consists of:

        (default_specified, default_value, default_value_type)

        Where:

            default_specified -- boolean indicating whether a default value
                                 was specified

            default_value -- the default value given, if any

            default_value_type -- the type of the default value (token.STRING,
                                  token.NAME, None). A type of None means
                                  unknown.
                                  
        """
        return self._parameter_info[paramName]
    
    def _extractFunctionParameters(self, tree):
        "Extract information about a function's parameters."
        dbg=0
        if dbg: print
        if dbg: print self._name
        function_parameters = []
        parameters = tree[3]
        if dbg: pprint.pprint(parameters)
        if parameters[1][0] != token.LPAR:
            raise 'Unrecognized parse result %s in %s' % (parameters[1],
                                                          parameters)
        if parameters[2][0] == token.RPAR:
            # No parameters: def func()
            return function_parameters
        if parameters[2][0] != symbol.varargslist:
            raise 'Unrecognized parse result %s in %s' % (parameters[2],
                                                          parameters)
        #
        # Move down the parse tree and process the argument list
        #
        parameters = parameters[2]
        if dbg: pprint.pprint(parameters)
        found_varargs = 0 # are we looking at a variable argument parameter?
        found_kwargs = 0  # are we looking at keyword argument parameter?
        name = None # what is the name of the parameter?
        found_default_value = None # did we see a default value for the param?
        default_value = None # what is the default value?
        default_value_type = None # what is the type of the default value?
        
        for parameter in parameters[1:]:

            # Shortcut cases
            if parameter[0] == token.COMMA:
                continue
            
            if parameter[0] == token.STAR:
                # Start variable argument definition
                found_varargs = 1

            if parameter[0] == token.DOUBLESTAR:
                # Start keyword argument definition
                found_kwargs = 1

            if (parameter[0] in (token.NAME, symbol.fpdef)) and name:
                # We've already found a name,
                # handle adding the previous
                # def to a list.
                function_parameters.append( (name,
                                             found_default_value,
                                             default_value,
                                             default_value_type) )
                name = found_default_value = None
                default_value = None
                default_value_type = None

            if parameter[0] == token.NAME:
                #
                # (Possibly fix and)
                # remember the new name
                #
                name = parameter[1]
                if found_varargs:
                    name = '*%s' % name
                elif found_kwargs:
                    name = '**%s' % name
                continue
                    
            if parameter[0] == symbol.fpdef:
                # Here we've made the concious decision
                # to include 'self' in the parameter list,
                # even if this is a method.  Renderers
                # will know (by context) whether we are
                # a method, and at that point can decide to
                # leave out the first parameter in the
                # paramter list.  This safeguards us from
                # [a] having to know whether we are a method
                # and [b] having to know whether the author
                # of the code used 'self' as the name of
                # 'self'.
                #
                name = parameter[1][1]
                continue

            if parameter[0] == token.EQUAL:
                #
                # Default value for the current parameter
                # coming up...
                #
                found_default_value = 1
                continue
            
            if parameter[0] == symbol.test:
                #
                # This is a parameter definition.
                #

                # Look for ARITH_EXPR parameter
                found, vars = lenientMatch(
                    #self.PARAMETER_DEFAULT_WITH_TRAILER_PATTERN,
                    self.PARAMETER_DEFAULT_PATTERN,
                    parameter,
                    #dbg=1
                    )
                
                if found:

                    if dbg: print 'FOUND: %s:' % name,
                    if dbg: pprint.pprint(vars)

                    if vars.has_key('term'):
                        default_value, default_value_type = \
                                       self._reconstructValueFromAtom(
                            vars['term'],
                            [
                            vars['trailer'],
                            vars['trailer_bits'],
                            ]
                            )

                else:
                    print 'UNRECOGNIZED:',
                    pprint.pprint(parameter)
                        
        ##
        ## <end> for parameter in parameters[1:]:
        ##
            
        if name:
            # Handle the last parameter
            #
            function_parameters.append( (name,
                                         found_default_value,
                                         default_value,
                                         default_value_type) )
        if dbg: print 'FOUND PARAMETERS: ',
        if dbg: pprint.pprint(function_parameters)
        return function_parameters

    
    def _reconstructValueFromAtom(self, atom, trailer=[]):
        """Convert an atom portion of a parse tree into the value.

        If the atom represents a string, number or name
        """
        dbg=0
        if dbg: print '\nRECONSTRUCTING VALUE FROM ATOM:',
        if dbg: pprint.pprint(atom)
        if trailer and dbg:
            print 'AND TRAILER:',
            pprint.pprint(trailer)
        if len(atom) == 2:
            if atom[1][0] == token.STRING:
                if dbg: print '\tSTRING'
                value = atom[1][1]
                value = value[1:-1]
                value_type = token.STRING
            elif atom[1][0] == token.NUMBER:
                if dbg: print '\tNUMBER'
                value = atom[1][1]
                value = eval(value)
                value_type = token.NUMBER
            elif atom[1][0] == token.NAME:
                if dbg: print '\tNAME'
                value = atom[1][1]
                value_type = token.NAME
                if value == 'None':
                    value = eval(value)
                else:
                    if trailer and filter(lambda x: x, trailer):
                        if dbg: print '\t\tVALUE: ', value
                        if dbg: print '\t\tTRAILER: ',
                        if dbg: pprint.pprint(trailer)
                        if dbg: print '\t\tSTRING TRAILER: "%s"' % parseTreeToString(trailer)
                        trailer_string = ''
                        for trailer_part in trailer:
                            if not trailer_part: continue
                            part_string = parseTreeToString(trailer_part)
                            trailer_string = joinCodeSnippets(trailer_string,
                                                              part_string,
                                                              ' ')
                        value = joinCodeSnippets(value, trailer_string, ' ')
                        value_type = None
            elif atom[1][0] == symbol.factor:
                if dbg: print '\tFACTOR'
                value = parseTreeToString(atom)
                value_type = None
                if trailer and filter(lambda x: x, trailer):
                    if dbg: print '\t\tVALUE: ', value
                    if dbg: print '\t\tTRAILER: ',
                    if dbg: pprint.pprint(trailer)
                    if dbg: print '\t\tSTRING TRAILER: "%s"' % parseTreeToString(trailer)
                    trailer_string = ''
                    for trailer_part in trailer:
                        if not trailer_part: continue
                        part_string = parseTreeToString(trailer_part)
                        trailer_string = joinCodeSnippets(trailer_string,
                                                          part_string,
                                                          ' ')
                    value = joinCodeSnippets(value, trailer_string, ' ')
                    value_type = None
            else:
                if dbg: print 'UNHANDLED SINGLETON: ', atom[1][0], ':',
                if dbg: pprint.pprint(atom)
                value = parseTreeToString(atom)
                value_type = None
        elif atom[1][0] in (token.LPAR, token.LSQB, token.LBRACE):
            # Convert the sequence back into a string
            # since that's good enough for representing
            # it in documentation.
            #
            if dbg: print '\tPARSE TREE TO STRING'
            value = parseTreeToString(atom)
            if dbg: print '\t', value
            value_type = None
        else:
            if dbg: print 'UNHANDLED MULTIPLE: ',
            if dbg: pprint.pprint(atom)
            value = parseTreeToString(atom)
            value_type = None
        if dbg: print '\tRETURNING: (%s, %s)' % (value, value_type)
        return value, value_type
    
