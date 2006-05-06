#!/usr/bin/env python
#
# $Id: hdpathmanagement.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $
#
# Time-stamp: <01/02/03 14:48:01 dhellmann>
#
# Copyright 2001 Doug Hellmann.
#
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

"""Provide a common set of path management functions.

"""

__rcs_info__ = {
    #
    #  Creation Information
    #
    'module_name'  : '$RCSfile: hdpathmanagement.py,v $',
    'rcs_id'       : '$Id: hdpathmanagement.py,v 1.1 2001/03/21 20:59:27 glandrum Exp $',
    'creator'      : 'Doug Hellmann <DougHellmann@bigfoot.com>',
    'project'      : 'UNSPECIFIED',
    'created'      : 'Sat, 03-Feb-2001 12:49:56 EST',

    #
    #  Current Information
    #
    'author'       : '$Author: glandrum $',
    'version'      : '$Revision: 1.1 $',
    'date'         : '$Date: 2001/03/21 20:59:27 $',
}

#
# Import system modules
#
import string
import os
import sys

#
# Import Local modules
#


#
# Module
#


def rmkdir(path):
    "Create a directory and all of its children."
    if not path:
        return
    parts = os.path.split(path)
    if len(parts) > 1:
        parent, child = parts
        if not (os.path.isdir(parent) or os.path.islink(parent)):
            rmkdir(parent)
    if not (os.path.isdir(path) or os.path.islink(path)):
        os.mkdir(path)
    return


def applyPrefixToPath(path, prefix):
    "Add the prefix value to every part of a given path."
    parts = string.split( path, os.sep )
    parts = map(lambda x, p=prefix: '%s%s' % (p, x), parts)
    name = apply(os.path.join, parts)
    return name

def removePrefix(path, prefix):
    "Remove prefix from the beginning of path, if present."
    common_prefix = os.path.commonprefix( (os.path.dirname(path),
                                           prefix)
                                          )
    if common_prefix == prefix:
        path = path[len(common_prefix):]
    while path and (path[0] == os.sep):
        path = path[1:]
    return path
    

def computeRelativeHTMLLink(fromName, toName, baseDirectory):
    """Compute the relative link between fromName and toName.

    Parameters

      'fromName' -- Named output node from which to compute the link.

      'toName' -- Named output node to which link should point.

      'baseDirectory' -- Name of the base directory in which both
      fromName and toName will reside.
      
    Both fromName and toName are strings refering to URL targets.
    This method computes the relative positions of the two nodes
    and returns a string which, if used as the HREF in a link in
    fromName will point directly to toName.
    """
    dbg=0
    if dbg: print '\nFROM: ', fromName
    if dbg: print 'TO  : ', toName
    if dbg: print 'BASE: ', baseDirectory

    #
    # Normalize directory names
    #
    fromName = os.path.normpath(fromName)
    toName = os.path.normpath(toName)
    if dbg: print 'FROM NORMALIZED : ', fromName
    if dbg: print 'TO NORMALIZED   : ', toName
    
    #
    # Remove the base directory prefix from both
    #
    fromName = removePrefix(fromName, baseDirectory)
    toName = removePrefix(toName, baseDirectory)
    if dbg: print 'FROM - PREFIX : ', fromName
    if dbg: print 'TO   - PREFIX : ', toName
    
    if fromName == toName:
        if dbg: print '\tsame name'
        relative_link = os.path.basename(toName)
    else:
        common_prefix = os.path.commonprefix( (os.path.dirname(fromName),
                                               os.path.dirname(toName))
                                              )
        #if not common_prefix:
        #    raise ValueError('No common prefix to names.', fromName, toName)
        from_name_no_prefix = fromName[len(common_prefix):]
        while from_name_no_prefix and (from_name_no_prefix[0] == os.sep):
            from_name_no_prefix = from_name_no_prefix[1:]
            if dbg: print '\tfrom, no prefix:', from_name_no_prefix
        subdir_path = os.path.dirname(from_name_no_prefix)
        parts = string.split(subdir_path, os.sep)
        if dbg: print '\tparts:', parts
        if parts and parts[0]:
            levels = len(string.split(subdir_path, os.sep))
        else:
            levels = 0
        up_levels = (os.pardir + os.sep) * levels
        to_name_no_prefix = toName[len(common_prefix):]
        if to_name_no_prefix and (to_name_no_prefix[0] == os.sep):
            to_name_no_prefix = to_name_no_prefix[1:]
        relative_link = "%s%s" % (up_levels, to_name_no_prefix)
    if dbg: print 'LINK: ', relative_link, '\n'
    return relative_link
