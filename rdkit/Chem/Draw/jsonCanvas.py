# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import re
import json  
from canvasbase import CanvasBase

class Canvas(CanvasBase):

#------------------------------------------------------------------------------

  def __init__(self, size=None):
    self.json=''
    self.objects = []
    self.size=size
    self.scale=1

#------------------------------------------------------------------------------

  def flush(self):
    self.json = json.dumps(self.objects)

#------------------------------------------------------------------------------
    
  def _addCanvasLine(self,p1,p2,color,**kwargs):
    line = dict()
    line['type'] = 'path'
    if kwargs.get('dash',(0,0)) == (0,0):
        line['path'] = "M%s,%sL%s,%s" % (p1[0], p1[1], p2[0], p2[1])
    else:
        dash = kwargs['dash']
        pts = self._getLinePoints(p1,p2,dash)     
        currDash = 0
        dashOn = True
        path = ''
        while currDash<(len(pts)-1):
          if dashOn:
            p1 = pts[currDash]
            p2 = pts[currDash+1]
            path += "M%s,%sL%s,%s" % (p1[0], p1[1], p2[0], p2[1])
          currDash+=1
          dashOn = not dashOn
        line['path'] = path  
        
    line['stroke'] = 'rgb' + str(color)
    line['stroke-width'] = self.scale*kwargs.get('linewidth',1)    
    return line  

#------------------------------------------------------------------------------

  def addCanvasLine(self,p1,p2,color=(0,0,0),color2=None,**kwargs):
    color = tuple(map(lambda x: x*255, color))
    color2 = tuple(map(lambda x: x*255, color2))
    if color2 and color2!=color:
      mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
      line1 = self._addCanvasLine(p1,mp,color,**kwargs)
      line2 = self._addCanvasLine(mp,p2,color2,**kwargs)
      self.objects.append(line1)
      self.objects.append(line2)
    else:
      line = self._addCanvasLine(p1,p2,color,**kwargs)
      self.objects.append(line)

#------------------------------------------------------------------------------    

  def addCanvasText(self,text,pos,font,color=(0,0,0),**kwargs):
    color = tuple(map(lambda x: x*255, color))
    bgColor = tuple(map(lambda x: x*255,kwargs.get('bgColor',(1,1,1))))
    tex = dict()
    tex['type'] = 'text'
    tex['x'] = pos[0]
    tex['y'] = pos[1]
    text = re.sub(r'\<.+?\>','',text)
    tex['text'] = text
    tex['font-size'] = font.size
    tex['stroke'] = 'rgb' + str(color)
    if font.weight=='bold':
        tex['font-weight'] = "bold"
    backgroundRect = dict()
    backgroundRect['type'] = 'rect'
    backgroundRect['x'] = pos[0] - (font.size * len(text)) / 2
    backgroundRect['y'] = pos[1] - font.size / 2
    backgroundRect['width'] = font.size * len(text)
    backgroundRect['height'] = font.size
    backgroundRect['fill'] = 'rgb' + str(bgColor)
    backgroundRect['stroke'] = 'rgb' + str(bgColor)
    self.objects.append(backgroundRect)    
    self.objects.append(tex)

#------------------------------------------------------------------------------
    
  def addCanvasPolygon(self,ps,color=(0,0,0),fill=True,stroke=False,**kwargs):
    if not fill and not stroke: return
    color = tuple(map(lambda x: x*255, color))
    polygon = dict()
    polygon['type'] = 'path'
    path = ''
    path += "M%s,%s" % (ps[0][0],ps[0][1])
    for p in ps[1:]:
      path += "L%s,%s" % (p[0],p[1])
    path += 'Z'  
    polygon['path'] = path 
    
    if fill:
        polygon['fill'] = 'rgb' + str(color)
    if stroke:
        polygon['stroke'] = 'rgb' + str(color)     
    self.objects.append(polygon)

#------------------------------------------------------------------------------      

  def addCanvasDashedWedge(self,p1,p2,p3,dash=(2,2),color=(0,0,0),
                           color2=None,**kwargs):
    color = tuple(map(lambda x: x*255, color))                      
    wedge = dict()
    wedge['stroke-width'] = self.scale*kwargs.get('linewidth',1)
    wedge['stroke'] = 'rgb' + str(color)
    wedge['type'] = 'path'
                     
    dash = (3,3)
    pts1 = self._getLinePoints(p1,p2,dash)
    pts2 = self._getLinePoints(p1,p3,dash)

    if len(pts2)<len(pts1): pts2,pts1=pts1,pts2

    path = ''
    for i in range(len(pts1)):
      path += "M%s,%s" % (pts1[i][0],pts1[i][1])
      path += "L%s,%s" % (pts2[i][0],pts2[i][1])

    wedge['path'] = path
    self.objects.append(wedge) 

#------------------------------------------------------------------------------

  def addCircle(self,center,radius,color=(0,0,0),fill=True,stroke=False,alpha=1.0,
                **kwargs):
    if not fill and not stroke: return
    color = tuple(map(lambda x: x*255, color))
    circle = dict()
    circle['type'] = 'circle'
    circle['cx'] = center[0]
    circle['cy'] = center[1]
    circle['r'] = radius
    if fill:
        circle['fill'] = 'rgb' + str(color)
    if stroke:
        circle['stroke'] = 'rgb' + str(color)
    if alpha:
        circle['opacity'] = alpha      
    self.objects.append(circle)
    
#------------------------------------------------------------------------------    
