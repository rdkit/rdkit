# $Id$
#
#  Copyright (C) 2002-2006  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
from qt import *
from icons import icons
import re,time,math
from xml.dom import minidom
try:
  from PIL import Image,ImageDraw
except ImportError:
  hasPil=0
else:
  hasPil=1
  
from StringIO import StringIO

_today=time.time()
_secsPerDay=60*60*24
_secsPerWeek=_secsPerDay*7

def TextifyDate(date,smart=1,cutoff=_secsPerWeek):
  """ date: a time tuple
  """
  try:
    delta = time.mktime(date) - _today
  except (OverflowError,TypeError):
    return ''
  nDays = int(math.ceil(delta/_secsPerDay))

  if not smart or nDays<-1:
    res = time.strftime('%d %B %Y',date)
  else:
    if delta < cutoff:
      if nDays == 0:
        res = 'Today'
      elif nDays == 1:
        res = "Tomorrow"
      elif nDays == -1:
        res = "Yesterday"
      else:
        res = '%d days'%(nDays)
    else:
      res = time.strftime('%d %B %Y',date)
  return res    
  
def GenCompletionPixmap(percent,size=(20,20),bgColor=(255,255,255),
                       doneColor=(255,100,100),undoneColor=(200,200,200),
                       outlineColor=(100,100,255)):
  pm = QPixmap()
  png = None
  if hasPil:
    img = GenCompletionImage(percent,size,bgColor,doneColor,undoneColor,outlineColor)
    sio = StringIO()
    try:
      img.save(sio,format='png')
      png = sio.getvalue()
    except:
      pass

  if png is None:
    idx = int(float(percent)/10)
    png = icons.pngs[idx]
  pm.loadFromData(png)
    
    

  return pm

def GenCompletionImage(percent,size=(20,20),bgColor=(255,255,255),
                       doneColor=(255,100,100),undoneColor=(200,200,200),
                       outlineColor=(100,100,255)):
  """ returns a PIL image

    percent can be either a fraction (0->1) or a percentage (0->100)
    
  """
  assert percent >= 0 and percent <= 100,"bad percent value: %f"%(percent)

  #if percent > 1:
  percent = float(percent)/100.

  img = Image.new('RGB',size,bgColor)
  startP = -90
  endP = -90 + percent*360
  if endP == startP:
    endP = startP+.1
  draw = ImageDraw.Draw(img)
  box = (0,0,int(size[0])-1,int(size[1])-1)
  draw.pieslice(box,int(startP),int(endP),
                fill=doneColor,outline=outlineColor)
  draw.pieslice(box,int(endP),int(startP),
                fill=undoneColor,outline=outlineColor)
  return img


def DelimitedToOPML(data,delimiter=','):
  """ converts a delimited text document (like from brainforest) to
    OPML
  """
  doc = minidom.Document()
  opml = doc.createElement('opml')
  opml.setAttribute('version','1.0')
  doc.appendChild(opml)
  head = doc.createElement('head')
  modDate='%s GMT'%time.asctime(time.gmtime())
  n = doc.createElement('dateModified')
  txt = doc.createTextNode(modDate)
  n.appendChild(txt)
  head.appendChild(n)
  opml.appendChild(head)
  body = doc.createElement('body')
  opml.appendChild(body)
  splitD = data.split('\r\n')
  if len(splitD) == 1:
    splitD = data.split('\n')

  # we're gonna be cute with this dictionary
  nodeDict = {}

  firstLine = splitD[0]
  headNode = doc.createElement('outline')
  headNode.setAttribute('text',firstLine)
  body.appendChild(headNode)
  nodeDict[''] = headNode

  lineNum = 1
  while lineNum < len(splitD):
    line = splitD[lineNum]
    if len(line):
      newNode = doc.createElement('outline')

      splitLine = line.split(delimiter)
      i = 0
      while len(splitLine[i]) == 0:
        i += 1
      startLev = i
      line = delimiter.join(splitLine[startLev:])

      # okay, we've cleaved off the opening delimiters, find the parent
      number = line.split(' ')[0]
      splitNum = number.split('.')[:-1]
      parentId = '.'.join(splitNum[:-1])
      ourId = '.'.join(splitNum)
      parent = nodeDict[parentId]

      # so far so good... get the text, priority and % done
      entry = re.findall(r'\ *.*\.\ *(.*)\ +\(p(.*),\ *(.*)\%',line)[0]
      text = entry[0]
      priority = entry[1]
      percentComplete = entry[2]

      if len(text.strip()) == 0:
        text = '-- No Name --'
      newNode.setAttribute('text',text)
      newNode.setAttribute('priority',priority)
      newNode.setAttribute('percentDone',percentComplete)

      entry = re.findall(r'\(.*\%\,\ due\ (.*)\s*$',line)
      if len(entry):
        due = entry[0]
        newNode.setAttribute('due',due)

      lineNum += 1
      entry = re.findall(r'<Note: (.*?)\>?\s*$',splitD[lineNum])
      if len(entry):
        note = entry[0]
        if not len(re.findall(r'\>\s*$',splitD[lineNum])):
          lineNum += 1
          while lineNum < len(splitD) and not len(re.findall('\>\s*$',splitD[lineNum])):
            note = '%s\n%s'%(note,splitD[lineNum])
            lineNum += 1
          toAdd = re.findall('(.*)\>\s*$',splitD[lineNum])
          note = '%s\n%s'%(note,toAdd[0])
        lineNum += 1
        newNode.setAttribute('note',note)
      nodeDict[ourId] = newNode
      parent.appendChild(newNode)
    else:
      lineNum += 1
  return doc.toxml()

BrainForestFormat=0
HumanFormat1=1
def _delimRecurse(node,header,level,delimiter,format=BrainForestFormat):
  res = []
  count = 1
  while node:
    if node.nodeType == node.ELEMENT_NODE and node.nodeName == 'outline':
      try:
        pri = node.getAttribute('priority')
      except:
        pri = 1
      txt = node.getAttribute('text')
      percentDone = node.getAttribute('percentDone')
      if format == BrainForestFormat:
        line = '%s%s%d. %s (p%s, %s%%)'%(delimiter*level,header,count,txt,pri,percentDone)
      elif format==HumanFormat1:
        if str(percentDone)=='':
          percentDone = "0"
        if str(pri) == '':  
          line = '%s%s%d. %s (%s%%)'%(delimiter*level,header,count,txt,percentDone)
        else:
          line = '%s%s%d. %s (Pri: %s, %s%%)'%(delimiter*level,header,count,
                                               txt,pri,percentDone)
      else:
        raise ValueError,'bad format'
        
      res.append(line)

      child = node.firstChild
      if child:
        below = _delimRecurse(child,'%s%d.'%(header,count),level+1,delimiter,
                              format=format)
        res.extend(below)
      count += 1
    node = node.nextSibling
  return res

def OPMLToDelimited(data,delimiter=',',header='from OPML',format=BrainForestFormat):
  """ converts OPML document to a delimited text document
  (like from brainforest).

  """
  res = []
  dom = minidom.parseString(data)
  root = dom.getElementsByTagName('body')[0]
  child = root.firstChild
  nodes = _delimRecurse(child,'',0,delimiter,format=format)
  res.extend(nodes)

  return '\n'.join(res)

if __name__ == '__main__':
  import sys
  fName = sys.argv[1]
  data =open(fName,'r').read()

  txt = OPMLToDelimited(data,format=HumanFormat1,delimiter="  ")
  print txt
