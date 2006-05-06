# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
import time,re,random,types
from xml.dom import minidom

def outlItemToDOM(item,headers=None,doc=None,parentNode=None,clearIDs=0):
  if headers is None:
    hdr = item.listView().header()
    headers = []
    for i in range(hdr.count()):
      headers.append(str(hdr.label(i)))
  if doc is None:
    doc = minidom.Document()
  node = doc.createElement("outline")
  for key in item.attribs().keys():
    if not clearIDs or key != 'id':
      val = str(item.attribs()[key])
      if val != "":
        node.setAttribute(key,val)

  childItem = item.firstChild()
  while childItem:
    childNode = outlItemToDOM(childItem,headers,doc,node)
    node.appendChild(childNode)
    childItem = childItem.nextSibling()

  return node

def findHelper(node,field,val,substring=0,caseSensitive=1):
  attrs = node.attribs()
  if attrs.has_key(field):
    v = attrs[field]
    if v == val:
      # easiest case
      return 1
    elif type(val) in types.StringTypes:
      if caseSensitive:
        val = val.upper()
        v = v.upper()
      if not substring:
        return val==v
      else:
        if v.find(val) >= 0:
          return 1
  return 0
