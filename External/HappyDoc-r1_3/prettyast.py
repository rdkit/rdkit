#!/usr/bin/env python
#
# Time-stamp: <00/07/04 13:26:36 dhellmann>
#
# Copyright 2000 Doug Hellmann
#
#                         All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that copyright notice and this permission
# notice appear in supporting documentation, and that the name of Doug
# Hellmann not be used in advertising or publicity pertaining to
# distribution of the software without specific, written prior
# permission.
#
# DOUG HELLMANN DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN
# NO EVENT SHALL DOUG HELLMANN BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

"Pretty print the AST for a .py file."

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: prettyast.py,v $',
    'rcs_id'       : '$Id: prettyast.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $',
    'creator'      : 'Doug Hellmann <doughellmann@bigfoot.com>',
    'project'      : 'HappyDoc',
    'created'      : 'Thu, 09-Mar-2000 13:41:44 EST',

    #
    #  Current Information
    #
    'author'       : '$Author: glandrum $',
    'version'      : '$Revision: 1.1 $',
    'date'         : '$Date: 2001/03/21 20:59:27 $',
    'locker'       : '$Locker:  $',
}

#
# Import system modules
#
import pprint
import os
import sys
import parser
import token, symbol
import types

#
# Import Local modules
#
from parseinfo import getDocs
from CommandLineApp import CommandLineApp

#
# Module
#

#
# Build a dictionary with all of the token and symbol values
#
nameDict = {}
nameDict.update(token.tok_name)
nameDict.update(symbol.sym_name)

def astListFixNames(astList):
    global nameDict
    if astList:
        name = nameDict.get(astList[0], astList[0])
        #print 'Replacing "%s" with %s' % (astList[0], name)
        astList[0] = name
    if (len(astList) == 2) and (type(astList[1]) == types.ListType):
        astListFixNames(astList[1])
    elif len(astList) > 2:
        for subList in astList[1:]:
            astListFixNames(subList)
    else:
        # do not recurse on non-lists
        pass
    return

class PrettyAST(CommandLineApp):

    output_directory = '.'
    
    def appInit(self):
        self.statusMessage('%s version %s' \
                           % (self.__class__.__name__, __rcs_info__['version']))
        return

    def prettyPrint(self, fileName):
        "Dump parse tree in readable format."
        source = open(fileName).read()
        basename = os.path.basename(os.path.splitext(fileName)[0])
        ast = parser.suite(source)
        l = parser.ast2list(ast)
        astListFixNames(l)
        pprint.pprint(l)
        return
    
    def main(self, *args):
        for file_name in args:
            #
            # Open the file
            #
            try:
                text = open(file_name).read()
            except IOError:
                self.errorMessage('Could not read %s' % file_name)
                continue
            else:
                self.statusMessage('Processing %s' % file_name, 2)
            self.prettyPrint(file_name)
        return


if __name__ == '__main__':
    try:
        PrettyAST().run()
    except PrettyAST.HelpRequested:
        pass

    
