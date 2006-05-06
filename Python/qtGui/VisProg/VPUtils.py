#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qtGui import qtUtils
import re
whitespace = re.compile('\s+')

def quickCDXMLParse(fName):
  """ quickly extracts some information from a CDXML file

    returns (nReactants,nProducts)
  """
  from xml.dom import minidom
  try:
    inF = open(fName,'r')
  except:
    qtUtils.error('cannot open reaction file %s for reading'%(fName))
    return 0,0
  dom = minidom.parseString(inF.read())
  steps = dom.getElementsByTagName('step')
  if len(steps) < 1:
    qtUtils.error('no steps in reaction')
    return None,None
  if len(steps) > 1:
    qtUtils.warning('%d steps in reaction, extras ignored'%len(steps))
  step = steps[0]

  attrStr = step.getAttribute('ReactionStepReactants')
  splits = whitespace.split(attrStr)
  nReact = len(splits) - splits.count('')

  attrStr = step.getAttribute('ReactionStepProducts')
  splits = whitespace.split(attrStr)
  nProd = len(splits) - splits.count('')

  return nReact,nProd

