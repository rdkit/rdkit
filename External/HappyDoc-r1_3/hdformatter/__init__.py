#!/usr/bin/env python
#
# Time-stamp: <00/12/28 10:31:55 dhellmann>
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


"""Formatter collection initialization module.

$Id: __init__.py,v 1.1 2001/03/21 21:58:54 glandrum Exp $


"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name':'$RCSfile: __init__.py,v $',
    'creator':'Doug Hellmann <doughellmann@home.com>',
    'project':'HappyDoc',
    'created':'Sat, 03-Jun-2000 19:06:39 EDT',
    #
    #  Current Information
    #
    'author':'$Author: glandrum $',
    'version':'$Revision: 1.1 $',
    'date':'$Date: 2001/03/21 21:58:54 $',
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
import hdformatter

class FormatterLoader(pluginloader.PluginLoader):
    "Load pluggable formatter types."

    def __init__(self):
        pluginloader.PluginLoader.__init__(self, hdformatter)
        return

    def addEntryPoint(self, infoDict):
        "Add the information about a docset to our lookup table."
        name = infoDict['name']
        factory = infoDict['factory']
        self[name] = factory
        return

plugins = FormatterLoader()

