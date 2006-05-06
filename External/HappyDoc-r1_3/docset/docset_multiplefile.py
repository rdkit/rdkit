#!/usr/bin/env python
#
# Time-stamp: <01/02/03 11:58:20 dhellmann>
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


"""Control class for writing out a set of documentation.

$Id: docset_multiplefile.py,v 1.1 2001/03/21 20:59:59 glandrum Exp $
"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: docset_multiplefile.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sun, 26-Mar-2000 11:19:54 EST',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.1 $',
    'date':'$Date: 2001/03/21 20:59:59 $',
    'locker':'$Locker:  $',
    }

#
# Import system modules
#
import UserList
import sys
import os
import string
import re
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO
import token
import symbol
import pprint
import glob
import parser

#
# Import Local modules
#
import optiontools

#
# Module
#

def entryPoint():
    "Return info about this module to the dynamic loader."
    return { 'name':'docset',
             'factory':DocSet,
             }

#
# Do not use multiple lines in the class docstring, since it shows up as
# part of the generated help message for the app and multiple lines throws
# off the formatting.
#
class DocSet(UserList.UserList):
    """Documentation is written to multiple files.

        Parameters

            includeComments -- Boolean.  False means to skip the
                               comment parsing step in the parser.
                               Default is True.
            
            includePrivateNames -- Boolean.  False means to ignore
                                   names beginning with _.  Default
                                   is True.

            prewrittenFileBasenames -- Base names (no extensions) of
                                       StructuredText files which are
                                       to be converted to the output
                                       format and included in the
                                       docset.

            statusMessageFunc -- function which will print a status
                                 message for the user

            title -- the title of the documentation set

            useRecursion -- Recurse into subdirectories looking for
                            subdirectories and files within them.

    """

    def __init__(self,
                 
                 formatterFactory,
                 parserFunc,
                 inputModuleNames,
                 
                 appName='HappyDoc',
                 appHome='',
                 appVersion='',
                 
                 author='',
                 baseDirectory=None,
                 descriptionFilename=None,
                 formatterParameters={},
                 ignoreDirFunc=None,
                 includeComments=1,
                 includePrivateNames=1,
                 prewrittenFileBasenames=( 'LICENSE',
                                           'CHANGES',
                                           'ANNOUNCE',
                                           'README',
                                           ),
                 statusMessageFunc=None,
                 title='HappyDoc Generated Documentation',
                 useRecursion=1,
                 
                 **extraNamedParameters
                 ):
        """Initialize the documentation set.

        Parameters

            formatterFactory -- a callable object which creates the type
                                of formatter class to be used for formatting
                                the documentation set.  The object will be
                                called and passed the DocSet along with other
                                configuration values for the formatter.

            parserFunc -- Parser function which returns the info for a module.

            inputModuleNames -- List of modules or directories to be documented.
            
            baseDirectory -- the name of the root directory where output
                             file(s) should be written

            descriptionFilename -- File describing the docset as a whole.

            formatterParameters -- other named configuration values to be passed
                                   to the formatter when it is created through the
                                   factory.  Any unrecognized values will be
                                   quietly ignored.

            ignoreDirFunc -- Function which returns true if the directory should
                             be ignored.

            *For others, see class documentation.*
                                   
        """        
        #
        # Store parameters
        #
        self._formatter_factory = formatterFactory
        self._parser_func = parserFunc
        self._input_module_names = inputModuleNames
        self._title = title
        self._base_directory = baseDirectory
        self._status_message_func = statusMessageFunc
        self._description_filename = descriptionFilename
        self._formatter_configuration = formatterParameters
        self._include_private_names = optiontools.getBooleanArgumentValue(includePrivateNames)
        self._include_comments = optiontools.getBooleanArgumentValue(includeComments)
        self._use_recursion = useRecursion
        self._ignore_dir_name = ignoreDirFunc
        self._prewritten_file_basenames = prewrittenFileBasenames
        self._prewritten_files = []
        
        #
        # Initialize this class
        #
        self._open_handles = []
        self._all_modules = {}
        self._all_classes = {}
        #
        # Initialize base class
        #
        UserList.UserList.__init__(self)
        self.statusMessage('Initializing documentation set %s...' % \
                           self._title)
        #
        # Handle unrecognized named parameters.
        #
        for extra_param, extra_value in extraNamedParameters.items():
            self.statusMessage(
                'WARNING: Parameter "%s" (%s) unrecognized by docset %s.' % \
                (extra_param, extra_value, self.__class__.__name__)
                )
        #
        # Process the modules specified.
        #
        self.processFiles(inputModuleNames)
        return
    
    def getFileInfo(self, fileName):
        "Parse the file and return the parseinfo instance for it."
        #self.statusMessage('Getting info for %s' % fileName)
        return self._parser_func(fileName, self._include_comments)

    def lookForPrewrittenFiles(self, dirName):
        """Look for prewritten StructuredText files in 'dirName'.
        """
        files = []
        for basename in self._prewritten_file_basenames:
            search_pat = os.path.join( dirName, '%s*' % basename )
            self.statusMessage('Looking for external documentation in %s' % search_pat,
                               2)
            found = glob.glob( search_pat )
            for f in found:
                #
                # Eliminate any backup files, etc. created by
                # text editors.
                #
                if f[-1] in '*~#':
                    continue
                #
                # Record this name
                #
                self.statusMessage(
                    '  found external documentation file\n    %s' \
                    % f, 2)
                files.append(f)
        return files
    
    def processFiles(self,
                     fileNames,
                     moduleFileName=re.compile(r'^.*\.py$').match,
                     ):
        "Process a list of files."
        for file_name in fileNames:
            if ( os.path.isdir(file_name)
                 and
                 (not self._ignore_dir_name(file_name))
                 and
                 (self._use_recursion >= 0)
                 ):
                #
                # Find pre-written files within the directory
                #
                self._prewritten_files = self._prewritten_files + \
                                         self.lookForPrewrittenFiles(file_name)
                #
                # Find modules and directories within to
                # recurse.
                #
                dir_contents = glob.glob('%s/*' % file_name)
                self.statusMessage('Recursing into %s' % file_name)
                if not self._use_recursion:
                    self._use_recursion = -1
                self.processFiles(dir_contents)
            elif moduleFileName(file_name):
                try:
                    file_info = self.getFileInfo(file_name)
                except parser.ParserError, msg:
                    self.statusMessage('ERROR: %s' \
                                       % string.join(msg.args, ' '))
                else:
                    self.append(file_info)
        return

    def getBaseDirectory(self):
        "Returns the base directory for this documentation set."
        return self._base_directory

    def getClassInfo(self, className):
        "Returns class info if have it, None otherwise."
        return self._all_classes.get(className, None)

    def statusMessage(self, message='', verboseLevel=1):
        "Print a status message for the user."
        if self._status_message_func:
            self._status_message_func(message, verboseLevel)
        return

    def getStructuredTextFile(self, filename):
        """Get the contents of a prewritten StructuredText file.
        """
        file_contents = ''
        try:
            file_contents = open(filename, 'rt').read()
        except IOError, msg:
            self.statusMessage( 'Could not open %s (%s)' % (filename, str(msg)))
        return file_contents
        

    def write(self):
        "Write the documentation set to the output."
        self.statusMessage('Beginning to write...')
        #
        # Create the formatter
        #
        self._formatter = apply( self._formatter_factory,
                                 ( self, ),
                                 self._formatter_configuration )
        #
        # Get the name of and open the docset root file
        #
        self._root_name = self._formatter.getFullOutputNameForObject(None)
        self._root_node = self.openOutput( self._root_name, self._title, '' )
        #
        # Get the description for the docset from the
        # file specified.
        #
        if (len(self._input_module_names) == 1) and os.path.isdir(self._input_module_names[0]):
            description_dir = self._input_module_names[0]
        else:
            description_dir = '.'
        self._writeDescription( description_dir, self._root_node )
        #
        # Write the output
        #
        self._writeModules()
        self._writePrewrittenFiles()
        self._writeTOC()
        self._writeIndex()
        #
        # Close things up
        #
        self.close()
        return

    def _writeDescription(self, directoryName, output):
        """Write the contents of the description file in 'directoryName' to 'output'.
        """
        if self._description_filename and (self._description_filename != '-'):
            description = None
            description_filename = os.path.join( directoryName,
                                                 self._description_filename )
            description = self.getStructuredTextFile( description_filename )
            if description:
                self._formatter.writeText( description, output, 0 )
        return
    
    def _writeTOC(self):
        "Output the TOC."
        self.statusMessage()
        self.statusMessage('Writing table of contents...')
        #
        # Open a new section and list
        #
        self._formatter.comment('BEGIN: TOC', self._root_node)
        self._formatter.pushSectionLevel(self._root_node)
        self._formatter.sectionHeader(self._root_node, 'Modules')
        #
        # Write module references to the list
        #
        #module_names = self._all_modules.keys()
        #module_names.sort()
        items = self._all_modules.items()
        module_ids = map(lambda x: (os.path.dirname(x[1].getFileName()),
                                    x[1].getFileName(),
                                    x[0]), items)
        module_ids.sort()
        modules_by_dir = {}
        for dir_name, file_name, module_name in module_ids:
            l = modules_by_dir.get(dir_name, [])
            l.append((file_name, module_name))
            modules_by_dir[dir_name] = l
        dirs = modules_by_dir.keys()
        dirs.sort()

        #
        # Figure out the name to the root node minus the base
        # directory.
        #
        output_reduced_name = self._root_node.name[len(self._base_directory)+1:]

        for dir_name in dirs:
            relative_dir_name = dir_name
            while relative_dir_name and (relative_dir_name[0] in './'):
                relative_dir_name = relative_dir_name[1:]
            #
            # Open the section for this directory
            #
            module_set = modules_by_dir[dir_name]
            extra_files = self._directories_with_extra_docs.get(relative_dir_name, [])
            if dir_name and (dir_name[-1] != '/'):
                dir_name = dir_name + '/'
            self._formatter.descriptiveListHeader(self._root_node,
                                                  dir_name)
            #
            # Write the list of extra files
            # which are in the directory being processed
            #
            for extra in extra_files:
                full_extra_name = os.path.join(dir_name, extra)
                while full_extra_name and (full_extra_name[0] in './'):
                    full_extra_name = full_extra_name[1:]
                #print 'FULL EXTRA NAME = "%s"' % full_extra_name
                #print 'EXTRA = "%s"' % extra
                self.statusMessage('\tAdding reference to %s to TOC' % extra)
                self._formatter.descriptiveListItem(
                    self._root_node,
                    self._formatter.getReference( full_extra_name,
                                                  output_reduced_name),
                    '')
            #
            # Write the list of modules
            # which are in the directory being processed
            #
            for file_name, module_name in module_set:
                module = self._all_modules[module_name]
                self.statusMessage('\tAdding %s to TOC' % module_name)
                self._formatter.descriptiveListItem(
                    self._root_node,
                    self._formatter.getReference(module,
                                                 output_reduced_name),
                    module.getSummary()
                    )
            self._formatter.descriptiveListFooter(self._root_node)
        #
        # Close the section
        #
        self._formatter.sectionFooter(self._root_node)
        self._formatter.popSectionLevel(self._root_node)
        self._formatter.comment('END: TOC', self._root_node)
        return

    def _writeIndex(self):
        "Output the index."
        self.statusMessage()
        self.statusMessage('IMPLEMENT Writing index...')
        return

    def _writePrewrittenFiles(self):
        """Convert the format of the discovered pre-written files.

        Convert the format of the discovered pre-written files.
        and write out new files as part of the docset.
        """
        directories_with_extra_docs = {}
        for filename in self._prewritten_files:
            body = self.getStructuredTextFile(filename)
            if not body:
                continue
            #
            # Convert the file using the formatter.
            #
            output_filename = self._formatter.getFullOutputNameForFile(
                filename, usePrefix=1 )
            short_output_filename = output_filename[len(self._base_directory)+1:]
            self.statusMessage('\tRewriting %s\n\t       to %s' \
                               % (filename, short_output_filename)
                               )
            output = self.openOutput( output_filename, self._title, '' )
            body = self.getStructuredTextFile( filename )
            self._formatter.writeText( body, output )
            self.closeOutput( output )
            #
            # Generate a reference to the file
            #
            relative_output_name = self._formatter.getOutputNameForFile( filename )
            dir, base = os.path.split(relative_output_name)
            file_list = directories_with_extra_docs.get( dir, [] )
            file_list.append( base )
            directories_with_extra_docs[dir] = file_list
        self._directories_with_extra_docs = directories_with_extra_docs
        return
    
    def _writeModules(self):
        "Output documentation for all modules."
        self.statusMessage()
        self.statusMessage('Writing module documentation...')
        module_names = self._all_modules.keys()
        module_names.sort()
        for module_name in module_names:
            self._writeModule( module_name )
        return

    def _filterNames(self, nameList):
        if not self._include_private_names:
            nameList = filter(lambda x: ( (x[0] != '_') or (x[:2] == '__') ),
                              nameList)
        return nameList

    def _describeClassInModuleNode(self, output, class_output_file_name, class_info):
        ref = self._formatter.getReference(class_info, class_output_file_name)
        self._formatter.descriptiveListItem(output, ref, class_info.getSummary())
        return
    
    def _writeModule(self, module_name):
        "Output the documentation for the module named."
        module = self._all_modules[module_name]
        output_name = self._formatter.getFullOutputNameForObject(module)
        output = self.openOutput(output_name,
                                 'Module: %s' % module_name,
                                 module.getFileName())
        self._formatter.comment('BEGIN: Module %s' % module_name, output)
        #
        # Write the doc string
        #
        self._formatter.writeText( module.getDocString(), output )
        #
        # Start the indented section
        #
        self._formatter.pushSectionLevel(output)
        #
        # Get some pre-formatted text we're going to reuse a lot
        #
        from_keyword = self._formatter.formatKeyword('from')
        import_keyword = self._formatter.formatKeyword('import')
        #
        # List the dependant modules
        #
        imported_modules = module.getImportData()
        if imported_modules:
            self._formatter.sectionHeader(output, 'Imported modules')
            self._formatter.listHeader(output, None, allowMultiColumn=0)
            output_reduced_name = output.name[len(self._base_directory)+1:]
            for name, symbols in imported_modules:
                if self._all_modules.has_key(name):
                    i_module = self._all_modules[name]
                    ref = self._formatter.getReference(
                        i_module,
                        output_reduced_name
                        )
                else:
                    ref = self._formatter.getPythonReference( name )
                    i_module = None
                if symbols:
                    if i_module:
                        #
                        # Process the list of imported names and
                        # generate references to the ones which are
                        # recognized by our parser.
                        #
                        import_list = []
                        for i_name in symbols:
                            try:
                                i_info = i_module.getClassInfo(i_name)
                                info_source = 'class'
                            except KeyError:
                                try:
                                    i_info = i_module.getFunctionInfo(i_name)
                                    info_source = 'function'
                                except KeyError:
                                    # do not have this name?
                                    import_list.append(i_name)
                                    continue

                            #
                            # Now that we have determined whether we have
                            # any information about the name, work out how
                            # to write the knowledge out.
                            #
                            output_reduced_name = output.name[len(
                                self._base_directory)+1:]
                            
                            if info_source == 'class':
                                #
                                # Get a reference to the file documenting the
                                # class.
                                #
                                i_ref = self._formatter.getReference(
                                    i_info,
                                    output_reduced_name,
                                    )
                                import_list.append(i_ref)
                                
                            elif info_source == 'function':
                                #
                                # Get a reference to the file documenting the
                                # module containing the function, and include
                                # the reference *into* the file to get right at
                                # the function.
                                #
                                i_ref = self._formatter.getNamedReference(
                                    i_module,
                                    i_name,
                                    output_reduced_name,
                                    )
                                import_list.append(i_ref)
                                
                            else:
                                #
                                # We do not know how to convert the name to a
                                # meaningful reference, so just put the name
                                # in the list to be output.
                                #
                                #print 'UNKNOWN TYPE: ', i_name
                                import_list.append(i_name)
                    else:
                        import_list = symbols
                    self._formatter.listItem( output,
                                               '%s %s %s %s' % \
                                               ( from_keyword,
                                                 ref,
                                                 import_keyword,
                                                 string.join(import_list, ', ')
                                                 )
                                               )
                else:
                    self._formatter.listItem( output,
                                               '%s %s' % (import_keyword, ref)
                                               )
            self._formatter.listFooter(output)
            self._formatter.sectionFooter(output)
        #
        # Write the info for the functions in this module
        #
        function_names = self._filterNames(module.getFunctionNames())
        if function_names:
            self._formatter.comment('BEGIN: Functions of module %s' % module_name, output)
            self._formatter.sectionHeader( output, 'Functions' )
            function_names.sort()
            #
            # TOC list
            #
            self._formatter.listHeader( output )
            for function_name in function_names:
                self._formatter.listItem(
                    output,
                    self._formatter.getInternalReference(
                    module.getFunctionInfo(function_name)
                    )
                    )
            self._formatter.listFooter( output )
            #
            # Function descriptions
            #
            for function_name in function_names:
                self._writeFunction(function_name,
                                    module.getFunctionInfo,
                                    output)
            self._formatter.sectionFooter(output)
            self._formatter.comment('END: Functions of module %s' % module_name, output)
        #
        # Write the info for the classes in this module
        #
        class_names = self._filterNames(module.getClassNames())
        if class_names:
            self._formatter.comment('BEGIN: Classes of module %s' % module_name, output)
            self._formatter.sectionHeader( output, 'Classes' )
            self._formatter.descriptiveListHeader( output, None )
            class_names.sort()
            output_reduced_name = output.name[len(self._base_directory)+1:]
            for class_name in class_names:
                c = module.getClassInfo(class_name)
                class_output_name = self._formatter.getFullOutputNameForObject(c)
                #print 'XXX writing docs for %s to %s' % (c, class_output_name)
                self._describeClassInModuleNode(output, class_output_name , c)
                class_output = self.openOutput(class_output_name,
                                               'Class: %s' % class_name,
                                               module.getFileName())
                self._writeClass( module, class_name, class_output )
                self.closeOutput(class_output)
            self._formatter.descriptiveListFooter(output)
            self._formatter.sectionFooter( output )
            self._formatter.comment('END: Classes of module %s' % module_name, output)
        #
        # Finish that indented level.
        #
        self._formatter.sectionFooter(output)
        self._formatter.popSectionLevel(output)
        self._formatter.comment('END: Module %s' % module_name, output)
        #
        # Close the output file
        #
        self.closeOutput(output)
        return


    
    def _writeBaseclassNames(self, parent, classInfo, output, indent=0):
        "Output the base class hierarchy for the given class."
        base_classes = classInfo.getBaseClassNames()
        if base_classes:
            if indent: self._formatter.indent(output)
            self._formatter.listHeader(output, None, allowMultiColumn=0)
            output_reduced_name = output.name[len(self._base_directory)+1:]
            for name in base_classes:
                try:
                    child = parent.getClassInfo(name)
                except KeyError:
                    self._formatter.listItem( output, name )
                else:
                    self._formatter.listItem(
                        output,
                        self._formatter.getReference(child,
                                                     output_reduced_name)
                        )
                    if name != classInfo.getName():
                        self._writeBaseclassNames(parent, child, output, 1)
            self._formatter.listFooter(output)
            if indent: self._formatter.dedent(output)
        return

    def _writeClass(self, parent, class_name, output):
        "Output the documentation for the class in the parent object."
        class_info = parent.getClassInfo(class_name)
        self._formatter.comment('BEGIN: Class %s' % class_name, output)
        self._formatter.writeText( class_info.getDocString(), output)
        #
        # Base class hierarchy
        #
        base_classes = class_info.getBaseClassNames()
        if base_classes:
            self._formatter.pushSectionLevel(output)
            self._formatter.sectionHeader(output, 'Base Classes')
            self._writeBaseclassNames(parent, class_info, output)
            self._formatter.sectionFooter(output)
            self._formatter.popSectionLevel(output)
        #
        # Start the indented section
        #
        self._formatter.pushSectionLevel(output)
        #
        # Write the info for the methods of this class
        #
        method_names = self._filterNames(class_info.getMethodNames())
        if method_names:
            self._formatter.sectionHeader( output, 'Methods' )
            method_names.sort()
            #
            # TOC list
            #
            self._formatter.listHeader( output )
            for method_name in method_names:
                self._formatter.listItem(
                    output,
                    self._formatter.getInternalReference(
                    class_info.getMethodInfo(method_name)
                    )
                    )
            self._formatter.listFooter( output )
            for method_name in method_names:
                self._writeFunction(method_name,
                                    class_info.getMethodInfo,
                                    output)
            self._formatter.sectionFooter(output)
        #
        # Finish that indented level.
        #
        self._formatter.sectionFooter(output)
        self._formatter.popSectionLevel(output)
        self._formatter.comment('END: Class %s' % class_name, output)
        return

    def _writeFunction(self, function_name, getInfo, output):
        "Output the documentation for the function in the parent object."
        function = getInfo(function_name)
        #
        # Header
        #
        self._formatter.itemHeader( output, function )
        #
        # Function signature
        #
        self._formatter.writeFunctionSignature( output, function )
        #
        # Docstring
        #
        self._formatter.writeText( function.getDocString(), output )
        #
        # Exceptions
        #
        exception_names = function.getExceptionNames()
        if exception_names:
            self._formatter.pushSectionLevel(output)
            self._formatter.sectionHeader(output, 'Exceptions')
            self._formatter.writeExceptionListForFunction(output, function, None)
            self._formatter.sectionFooter(output)
            self._formatter.popSectionLevel(output)
        return

    def close(self):
        "Close the open documentation set."
        for f in self._open_handles:
            try:
                self.closeOutput(f)
            except:
                pass
        return

    def append(self, infoObject):
        self._all_modules[ infoObject.getName() ] = infoObject
        for c in infoObject.getClassNames():
            self._all_classes[ c ] = infoObject.getClassInfo(c)
        UserList.UserList.append(self, infoObject)
        return

    def openOutput(self, name, title, subtitle):
        """
        Using this method to open output destinations
        means they will automatically be closed.
        """
        self.statusMessage(
            '\tDocumenting : "%s"' % title)
        f = self._formatter.openOutput(name, title, subtitle)
        self._open_handles.append(f)
        return f
    
    def closeOutput(self, output):
        "Close the output handle."
        self._formatter.closeOutput(output)
        return

    

    
