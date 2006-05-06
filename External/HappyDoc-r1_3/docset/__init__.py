#!/usr/bin/env python
#
# Time-stamp: <01/01/13 18:07:30 dhellmann>
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


"""Docset collection initialization module.

$Id: __init__.py,v 1.1 2001/03/21 20:59:59 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: __init__.py,v $',
    'creator':'Doug Hellmann <doughellmann@bigfoot.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 19:06:39 EDT',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.1 $',
    'date':'$Date: 2001/03/21 20:59:59 $',
    }

#
# Import system modules
#
import os
import sys
import glob
import pprint


#
# Import Local modules
#
import pluginloader

#
# Module
#
import docset

class DocSetLoader(pluginloader.PluginLoader):
    "Load pluggable docset types."

    def __init__(self):
        pluginloader.PluginLoader.__init__(self, docset)
        return

    def addEntryPoint(self, infoDict):
        "Add the information about a docset to our lookup table."
        name = infoDict['name']
        factory = infoDict['factory']
        self[name] = factory
        return

plugins = DocSetLoader()

