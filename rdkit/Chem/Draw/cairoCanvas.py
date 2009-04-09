#
#
#  Copyright (C) 2008 Greg Landrum
#  Copyright (C) 2009 Uwe Hoffmann
#
#   @@ All Rights Reserved  @@
#
import cairo 
import math
import rdkit.RDConfig
import os
import array

#missing up to now
#faceMap={'sans':os.path.join(
#    rdkit.RDConfig.RDCodeDir,'Chem','Draw','FreeSans.ttf')}

class Canvas(object):
  def __init__(self,
               image=None,  # PIL image
               size=None,
               ctx=None,
               imageType=None, # determines file type
               fileName=None,  # if set determines output file name
               ):
    """
    Canvas can be used in four modes:
    1) using the supplied PIL image
    2) using the supplied cairo context ctx
    3) writing to a file fileName with image type imageType
    4) creating a cairo surface and context within the constructor
    """
    self.image=None
    self.imageType=imageType
    if image is not None:
      imgd = image.tostring()
      a = array.array('B',imgd)
      stride = cairo.ImageSurface.format_stride_for_width (cairo.FORMAT_ARGB32,
                                                           image.size[0])
      if stride != image.size[0] * 4:
        raise Exception ,"invalid stride"
      surface = cairo.ImageSurface.create_for_data (
        a, cairo.FORMAT_ARGB32,
        image.size[0], image.size[1], stride)
      ctx = cairo.Context(surface)
      size=image.size[0], image.size[1]
      self.image=image
    elif size is not None:
      ##cairo.HAS_GLITZ_SURFACE
      #cairo.HAS_IMAGE_SURFACE
      ##cairo.HAS_QUARTZ_SURFACE
      ##cairo.HAS_USER_FONT
      #cairo.HAS_WIN32_SURFACE
      #cairo.HAS_XCB_SURFACE
      ##cairo.HAS_XLIB_SURFACE
      #cairo.HAS_PNG_FUNCTIONS
      #print imageType
      if cairo.HAS_PDF_SURFACE and imageType == "pdf":
        surface = cairo.PDFSurface (fileName, size[0], size[1])
      elif cairo.HAS_SVG_SURFACE and imageType == "svg":
        surface = cairo.SVGSurface (fileName, size[0], size[1])
      elif cairo.HAS_PS_SURFACE and imageType == "ps":
        surface = cairo.PSSurface (fileName, size[0], size[1])
      #elif cairo.HAS_IMAGE_SURFACE and imageType == "png":
      elif imageType == "png":
        surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, size[0], size[1])
      else:
        raise ValueError, "Unrecognized file type. Valid choices are pdf, svg, ps, and png"
      ctx = cairo.Context(surface)
      ctx.set_source_rgb(255,255,255)
      ctx.paint()
    else:
      surface=ctx.get_target()
      try:
        size=surface.get_width(),surface.get_height() 
      except AttributeError:
        size=None
    self.ctx=ctx
    self.size=size
    self.surface=surface
    self.fileName=fileName


  def flush(self):
    """temporary interface, must be splitted to different methods,
       (elif self.imageType == "png":)"""
    if self.fileName and self.imageType=='png':
      self.surface.write_to_png(self.fileName)
    elif self.image is not None:
      self.image.fromstring(self.surface.get_data(),
                            "raw","RGBA",0,1)
      self.surface.finish()
    elif self.imageType == "png":
      buffer=self.surface.get_data()
      return buffer


def convertColor(color):
  color = (int(color[0]*255),int(color[1]*255),int(color[2]*255))
  return color

def _getLinePoints(p1,p2,dash):
  x1,y1=p1
  x2,y2=p2
  dx = x2-x1
  dy = y2-y1
  lineLen = math.sqrt(dx*dx+dy*dy)
  theta = math.atan2(dy,dx)
  cosT = math.cos(theta)
  sinT = math.sin(theta)

  pos = (x1,y1)
  pts = [pos]
  dist = 0
  currDash = 0
  while dist < lineLen:
    currL = dash[currDash%len(dash)]
    if(dist+currL > lineLen): currL = lineLen-dist
    endP = (pos[0] + currL*cosT, pos[1] + currL*sinT)
    pts.append(endP)
    pos = endP
    dist += currL
    currDash += 1
  return pts

def _doLine(canvas,p1,p2,**kwargs):
  if kwargs.get('dashes',(0,0)) == (0,0):
    canvas.ctx.move_to(p1[0],p1[1])
    canvas.ctx.line_to(p2[0],p2[1])
  else:
    # the antialiasing makes the dashes appear too small
    dash = [x*4 for x in kwargs['dashes']]
    pts = _getLinePoints(p1,p2,dash)

    currDash = 0
    dashOn = True
    while currDash<(len(pts)-1):
      if dashOn:
        p1 = pts[currDash]
        p2 = pts[currDash+1]
        canvas.ctx.move_to(p1[0],p1[1])
        canvas.ctx.line_to(p2[0],p2[1])
      currDash+=1
      dashOn = not dashOn

def addCanvasLine(canvas,p1,p2,color=(0,0,0),color2=None,**kwargs):
  canvas.ctx.set_line_width(kwargs.get('linewidth',1))
  if color2 and color2!=color:
    mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
    color = convertColor(color)
    canvas.ctx.set_source_rgb(*color)
    _doLine(canvas,p1,mp,**kwargs)
    canvas.ctx.stroke()
    color2 = convertColor(color2)
    canvas.ctx.set_source_rgb(*color2)
    _doLine(canvas,mp,p2,**kwargs)
    canvas.ctx.stroke()
  else:
    color = convertColor(color)
    canvas.ctx.set_source_rgb(*color)
    _doLine(canvas,p1,p2,**kwargs)
    canvas.ctx.stroke()

def addCanvasText(canvas,text,pos,font,color=(0,0,0),**kwargs):
  color = convertColor(color)
  canvas.ctx.select_font_face("Georgia",
                              cairo.FONT_SLANT_NORMAL,
                              cairo.FONT_WEIGHT_BOLD)
  canvas.ctx.set_font_size(font.size)
  w,h=canvas.ctx.text_extents(text)[2:4]
  bw,bh=w*1.8,h*1.4
  dPos = pos[0]-bw/2.,pos[1]-bh/2.
  bgColor=kwargs.get('bgColor',(1,1,1))
  bgColor = convertColor(bgColor)
  canvas.ctx.set_source_rgb(*bgColor)
  canvas.ctx.rectangle(dPos[0],dPos[1],bw,bh)
  canvas.ctx.fill()
  dPos = pos[0]-w/2.,pos[1]+h/2.
  canvas.ctx.set_source_rgb(*color)
  canvas.ctx.move_to(*dPos)
  canvas.ctx.show_text(text)

def addCanvasPolygon(canvas,ps,color=(0,0,0),**kwargs):
  dps = []
  color = convertColor(color)
  canvas.ctx.set_source_rgb(*color)
  canvas.ctx.move_to(ps[0][0],ps[0][1])
  for p in ps[1:]:
    canvas.ctx.line_to(p[0],p[1])
  canvas.ctx.close_path()
  canvas.ctx.fill()

def addCanvasDashedWedge(canvas,p1,p2,p3,dash=(2,2),color=(0,0,0),
                         color2=None,**kwargs):
  canvas.ctx.set_line_width(kwargs.get('linewidth',1))
  color = convertColor(color)
  canvas.ctx.set_source_rgb(*color)
  dash = (3,3)
  pts1 = _getLinePoints(p1,p2,dash)
  pts2 = _getLinePoints(p1,p3,dash)

  if len(pts2)<len(pts1): pts2,pts1=pts1,pts2

  for i in range(len(pts1)):
    canvas.ctx.move_to(pts1[i][0],pts1[i][1])
    canvas.ctx.line_to(pts2[i][0],pts2[i][1])
  canvas.ctx.stroke()
