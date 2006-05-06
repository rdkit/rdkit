#!/usr/bin/env python
#
# Time-stamp: <01/02/03 12:11:42 dhellmann>
#
# COPYRIGHT
#
#   Permission to use, copy, modify, and distribute this software and
#   its documentation for any purpose and without fee is hereby
#   granted, provided that the above copyright notice appear in all
#   copies and that both that copyright notice and this permission
#   notice appear in supporting documentation, and that the name of Doug
#   Hellmann not be used in advertising or publicity pertaining to
#   distribution of the software without specific, written prior
#   permission.
# 
# DISCLAIMER
#
#   DOUG HELLMANN DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
#   INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN
#   NO EVENT SHALL DOUG HELLMANN BE LIABLE FOR ANY SPECIAL, INDIRECT OR
#   CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
#   OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
#   NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
#   CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 


"""Define a class to handle pluggable module loading.

$Id: pluginloader.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: pluginloader.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 19:23:48 EDT',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.1 $',
    'date':'$Date: 2001/03/21 20:59:27 $',
    }

#
# Import system modules
#
import sys
import glob
import os
import pprint
import UserDict

#
# Import Local modules
#


#
# Module
#

#
# Where is the application installed?
#
APP_HOME_DIR=os.path.dirname( sys.argv[0] )
sys.path.append(APP_HOME_DIR)

class PluginLoader(UserDict.UserDict):
    """A class to handle pluggable module loading."""
    
    def __init__(self, module):
        """Create a PluginLoader.

        Parameters

          module -- The imported module, from which we get name, path,
          etc. to find the sub-modules.
          
        """
        #
        # Initialize the base class
        #
        UserDict.UserDict.__init__(self)
        #
        # Name for these plugins
        #
        self.plugin_set_name = module.__name__
        #
        # Where the plugins will be
        #
        self.plugin_dir = module.__path__[0]
        #
        # List of plugins
        #
        _module_list = glob.glob( os.path.join( self.plugin_dir,
                                                '%s*py' % self.plugin_set_name,
                                                )
                                  )
        _module_list.sort()
        #
        # Load the modules
        #
        for _module_name in _module_list:
            _module_name = os.path.basename(_module_name)[:-3]
            _import_name = '%s.%s' % ( self.plugin_set_name,
                                       _module_name )
            try:
                _module = __import__( _import_name )
            except:
                sys.stderr.write('ERROR: Could not import %s' % _import_name)
                continue
            _module = getattr(_module, _module_name)
            try:
                info = _module.entryPoint()
            except AttributeError:
                pass
            else:
                self.addEntryPoint(info)
        return

    def addEntryPoint(self, infoDict):
        """Add an entry point into a module to our lookup tables.

        This method must be implemented by the subclass.
        """
        raise 'Not implemented for %s' % self.__class__.__name__

    
# _docset_module_list = glob.glob( os.path.join(DOCSET_DIR, 'docset*.py') )
# _docset_module_list.sort()
# for _docset_module_name in _docset_module_list:
#     #print 'FILE NAME:', _docset_module_name
#     _docset_module_name = os.path.basename(_docset_module_name)[:-3]
#     #print 'MODULE NAME:', _docset_module_name
#     _docset_module = __import__('docset.%s' % _docset_module_name,
#                                 {}, # globals
#                                 {}, # locals
#                                 [], # from list
#                                 )
#     #print 'MODULE:', _docset_module,
#     #pprint.pprint(_docset_module.__dict__)
#     _docset_module = getattr(_docset_module, _docset_module_name)
#     #print 'MODULE 2:', _docset_module,
#     #pprint.pprint(_docset_module.__dict__)
#     info = _docset_module.getInfo()
#     #print info
#     sys.exit(0)
