# $Id$
#
#  Copyright (C) 2008 Greg Landrum
#  Copyright (C) 2009 Uwe Hoffmann
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import cairo 
import math
import rdkit.RDConfig
import os
import array

from canvasbase import CanvasBase

class Canvas(CanvasBase):
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
      imgd = image.tostring("raw","BGRA")
      a = array.array('B',imgd)
      stride=image.size[0]*4
      surface = cairo.ImageSurface.create_for_data (
        a, cairo.FORMAT_ARGB32,
        image.size[0], image.size[1], stride)
      ctx = cairo.Context(surface)
      size=image.size[0], image.size[1]
      self.image=image
    elif size is not None:
      if cairo.HAS_PDF_SURFACE and imageType == "pdf":
        surface = cairo.PDFSurface (fileName, size[0], size[1])
      elif cairo.HAS_SVG_SURFACE and imageType == "svg":
        surface = cairo.SVGSurface (fileName, size[0], size[1])
      elif cairo.HAS_PS_SURFACE and imageType == "ps":
        surface = cairo.PSSurface (fileName, size[0], size[1])
      elif imageType == "png":
        surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, size[0], size[1])
      else:
        raise ValueError, "Unrecognized file type. Valid choices are pdf, svg, ps, and png"
      ctx = cairo.Context(surface)
      ctx.set_source_rgb(1,1,1)
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
    """
    if self.fileName and self.imageType=='png':
      self.surface.write_to_png(self.fileName)
    elif self.image is not None:
      # on linux at least it seems like the PIL images are BGRA, not RGBA:
      self.image.fromstring(self.surface.get_data(),
                            "raw","BGRA",0,1)
      self.surface.finish()
    elif self.imageType == "png":
      buffer=self.surface.get_data()
      return buffer

  def _doLine(self, p1, p2, **kwargs):
    if kwargs.get('dash',(0,0)) == (0,0):
      self.ctx.move_to(p1[0],p1[1])
      self.ctx.line_to(p2[0],p2[1])
    else:
      # the antialiasing makes the dashes appear too small
      dash = [x*4 for x in kwargs['dash']]
      pts = self._getLinePoints(p1,p2,dash)

      currDash = 0
      dashOn = True
      while currDash<(len(pts)-1):
        if dashOn:
          p1 = pts[currDash]
          p2 = pts[currDash+1]
          self.ctx.move_to(p1[0],p1[1])
          self.ctx.line_to(p2[0],p2[1])
        currDash+=1
        dashOn = not dashOn

  def addCanvasLine(self,p1,p2,color=(0,0,0),color2=None,**kwargs):
    self.ctx.set_line_width(kwargs.get('linewidth',1))
    if color2 and color2!=color:
      mp = (p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.
      self.ctx.set_source_rgb(*color)
      self._doLine(p1,mp,**kwargs)
      self.ctx.stroke()
      self.ctx.set_source_rgb(*color2)
      self._doLine(mp,p2,**kwargs)
      self.ctx.stroke()
    else:
      self.ctx.set_source_rgb(*color)
      self._doLine(p1,p2,**kwargs)
      self.ctx.stroke()

  def addCanvasText(self,text,pos,font,color=(0,0,0),**kwargs):
    self.ctx.select_font_face(font.face,
                                cairo.FONT_SLANT_NORMAL,
                                cairo.FONT_WEIGHT_BOLD)
    self.ctx.set_font_size(font.size)
    w,h=self.ctx.text_extents(text)[2:4]
    bw,bh=w*1.8,h*1.4
    dPos = pos[0]-bw/2.,pos[1]-bh/2.
    bgColor=kwargs.get('bgColor',(1,1,1))
    self.ctx.set_source_rgb(*bgColor)
    self.ctx.rectangle(dPos[0],dPos[1],bw,bh)
    self.ctx.fill()
    dPos = pos[0]-w/2.,pos[1]+h/2.
    self.ctx.set_source_rgb(*color)
    self.ctx.move_to(*dPos)
    self.ctx.show_text(text)

  def addCanvasPolygon(self,ps,color=(0,0,0),fill=True,stroke=False,**kwargs):
    if not fill and not stroke: return
    dps = []
    self.ctx.set_source_rgb(*color)
    self.ctx.move_to(ps[0][0],ps[0][1])
    for p in ps[1:]:
      self.ctx.line_to(p[0],p[1])
    self.ctx.close_path()
    if stroke:
      if fill:
        self.ctx.stroke_preserve()
      else:
        self.ctx.stroke()
    if fill:
      self.ctx.fill()

  def addCanvasDashedWedge(self,p1,p2,p3,dash=(2,2),color=(0,0,0),
                           color2=None,**kwargs):
    self.ctx.set_line_width(kwargs.get('linewidth',1))
    self.ctx.set_source_rgb(*color)
    dash = (3,3)
    pts1 = self._getLinePoints(p1,p2,dash)
    pts2 = self._getLinePoints(p1,p3,dash)

    if len(pts2)<len(pts1): pts2,pts1=pts1,pts2

    for i in range(len(pts1)):
      self.ctx.move_to(pts1[i][0],pts1[i][1])
      self.ctx.line_to(pts2[i][0],pts2[i][1])
    self.ctx.stroke()
