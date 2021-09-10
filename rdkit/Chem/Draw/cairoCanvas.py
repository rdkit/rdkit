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
import array
import math
import os
import re

from PIL import Image

from rdkit.Chem.Draw.canvasbase import CanvasBase

have_cairocffi = False
# for Python3, import cairocffi preferably
try:
  import cairocffi as cairo
except ImportError:
  import cairo
else:
  have_cairocffi = True
have_pango = False
if 'RDK_NOPANGO' not in os.environ:
  if have_cairocffi:
    import cffi
    import platform
    ffi = cffi.FFI()
    ffi.cdef('''
        /* GLib */
        typedef void* gpointer;
        typedef void cairo_t;
        typedef void PangoFontDescription;
        void g_object_unref (gpointer object);

        /* Pango and PangoCairo */
        #define PANGO_SCALE 1024
        typedef ... PangoLayout;
        typedef enum {
            PANGO_ALIGN_LEFT,
            PANGO_ALIGN_CENTER,
            PANGO_ALIGN_RIGHT
        } PangoAlignment;
        typedef struct PangoRectangle {
          int x;
          int y;
          int width;
          int height;
        } PangoRectangle;
        PangoLayout *pango_cairo_create_layout (cairo_t *cr);
        void pango_cairo_update_layout (cairo_t *cr, PangoLayout *layout);
        void pango_cairo_show_layout (cairo_t *cr, PangoLayout *layout);
        void pango_layout_set_alignment (
            PangoLayout *layout, PangoAlignment alignment);
        void pango_layout_set_markup (
            PangoLayout *layout, const char *text, int length);
        void pango_layout_get_pixel_extents (PangoLayout *layout,
            PangoRectangle *ink_rect, PangoRectangle *logical_rect);
        PangoFontDescription *pango_font_description_new (void);
        void pango_font_description_free (PangoFontDescription *desc);
        void pango_font_description_set_family (PangoFontDescription *desc,
            const char *family);
        void pango_font_description_set_size (PangoFontDescription *desc,
            int size);
        void pango_layout_set_font_description (PangoLayout *layout,
            const PangoFontDescription *desc);
    ''')
    if platform.system() == 'Windows':
      defaultLibs = {
        'pango_default_lib': 'libpango-1.0-0.dll',
        'pangocairo_default_lib': 'libpangocairo-1.0-0.dll',
        'gobject_default_lib': 'libgobject-2.0-0.dll'
      }
    else:
      defaultLibs = {
        'pango_default_lib': 'pango-1.0',
        'pangocairo_default_lib': 'pangocairo-1.0',
        'gobject_default_lib': 'gobject-2.0'
      }
    import ctypes.util
    for libType in ['pango', 'pangocairo', 'gobject']:
      envVar = 'RDK_' + libType.upper() + '_LIB'
      envVarSet = False
      if envVar in os.environ:
        envVarSet = True
        libName = os.environ[envVar]
      else:
        libName = defaultLibs[libType + '_default_lib']
      libPath = ctypes.util.find_library(libName)
      exec(libType + ' = None')
      importError = False
      if libPath:
        try:
          exec(libType + ' = ffi.dlopen("' + libPath.replace('\\', '\\\\') + '")')
        except:
          if envVarSet:
            importError = True
          else:
            pass
      else:
        importError = True
      if importError:
        raise ImportError(envVar + ' set to ' + libName + ' but ' + libType.upper() +
                          ' library cannot be loaded.')
    have_pango = (pango and pangocairo and gobject)
  else:
    for libType in ['pango', 'pangocairo']:
      try:
        exec('import ' + libType)
      except ImportError:
        exec(libType + ' = None')
    have_pango = (pango and pangocairo)

if (not hasattr(cairo.ImageSurface, 'get_data')
    and not hasattr(cairo.ImageSurface, 'get_data_as_rgba')):
  raise ImportError('cairo version too old')

scriptPattern = re.compile(r'\<.+?\>')


class Canvas(CanvasBase):

  def __init__(
      self,
      image=None,  # PIL image
      size=None,
      ctx=None,
      imageType=None,  # determines file type
      fileName=None,  # if set determines output file name
  ):
    """
        Canvas can be used in four modes:
        1) using the supplied PIL image
        2) using the supplied cairo context ctx
        3) writing to a file fileName with image type imageType
        4) creating a cairo surface and context within the constructor
        """
    self.image = None
    self.imageType = imageType
    if image is not None:
      try:
        imgd = image.tobytes("raw", "BGRA")
      except SystemError:
        r, g, b, a = image.split()
        mrg = Image.merge("RGBA", (b, g, r, a))
        imgd = mrg.tobytes("raw", "RGBA")

      a = array.array('B', imgd)
      stride = image.size[0] * 4
      surface = cairo.ImageSurface.create_for_data(a, cairo.FORMAT_ARGB32, image.size[0],
                                                   image.size[1], stride)
      ctx = cairo.Context(surface)
      size = image.size[0], image.size[1]
      self.image = image
    elif ctx is None and size is not None:
      if hasattr(cairo, "PDFSurface") and imageType == "pdf":
        surface = cairo.PDFSurface(fileName, size[0], size[1])
      elif hasattr(cairo, "SVGSurface") and imageType == "svg":
        surface = cairo.SVGSurface(fileName, size[0], size[1])
      elif hasattr(cairo, "PSSurface") and imageType == "ps":
        surface = cairo.PSSurface(fileName, size[0], size[1])
      elif imageType == "png":
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, size[0], size[1])
      else:
        raise ValueError("Unrecognized file type. Valid choices are pdf, svg, ps, and png")
      ctx = cairo.Context(surface)
      ctx.set_source_rgb(1, 1, 1)
      ctx.paint()
    else:
      surface = ctx.get_target()
      if size is None:
        try:
          size = surface.get_width(), surface.get_height()
        except AttributeError:
          size = None
    self.ctx = ctx
    self.size = size
    self.surface = surface
    self.fileName = fileName

  def flush(self):
    """temporary interface, must be splitted to different methods,
        """
    if self.fileName and self.imageType == 'png':
      self.surface.write_to_png(self.fileName)
    elif self.image is not None:
      # on linux at least it seems like the PIL images are BGRA, not RGBA:
      if hasattr(self.surface, 'get_data'):
        self.image.frombytes(bytes(self.surface.get_data()), "raw", "BGRA", 0, 1)
      else:
        self.image.frombytes(bytes(self.surface.get_data_as_rgba()), "raw", "RGBA", 0, 1)
      self.surface.finish()
    elif self.imageType == "png":
      if hasattr(self.surface, 'get_data'):
        buffer = self.surface.get_data()
      else:
        buffer = self.surface.get_data_as_rgba()
      return buffer

  def _doLine(self, p1, p2, **kwargs):
    if kwargs.get('dash', (0, 0)) == (0, 0):
      self.ctx.move_to(p1[0], p1[1])
      self.ctx.line_to(p2[0], p2[1])
    else:
      dash = kwargs['dash']
      pts = self._getLinePoints(p1, p2, dash)

      currDash = 0
      dashOn = True
      while currDash < (len(pts) - 1):
        if dashOn:
          p1 = pts[currDash]
          p2 = pts[currDash + 1]
          self.ctx.move_to(p1[0], p1[1])
          self.ctx.line_to(p2[0], p2[1])
        currDash += 1
        dashOn = not dashOn

  def addCanvasLine(self, p1, p2, color=(0, 0, 0), color2=None, **kwargs):
    self.ctx.set_line_width(kwargs.get('linewidth', 1))
    if color2 and color2 != color:
      mp = (p1[0] + p2[0]) / 2., (p1[1] + p2[1]) / 2.
      self.ctx.set_source_rgb(*color)
      self._doLine(p1, mp, **kwargs)
      self.ctx.stroke()
      self.ctx.set_source_rgb(*color2)
      self._doLine(mp, p2, **kwargs)
      self.ctx.stroke()
    else:
      self.ctx.set_source_rgb(*color)
      self._doLine(p1, p2, **kwargs)
      self.ctx.stroke()

  def _addCanvasText1(self, text, pos, font, color=(0, 0, 0), **kwargs):
    if font.weight == 'bold':
      weight = cairo.FONT_WEIGHT_BOLD
    else:
      weight = cairo.FONT_WEIGHT_NORMAL
    self.ctx.select_font_face(font.face, cairo.FONT_SLANT_NORMAL, weight)
    text = scriptPattern.sub('', text)
    self.ctx.set_font_size(font.size)
    w, h = self.ctx.text_extents(text)[2:4]
    bw, bh = w + h * 0.4, h * 1.4
    offset = w * pos[2]
    dPos = pos[0] - w / 2. + offset, pos[1] + h / 2.
    self.ctx.set_source_rgb(*color)
    self.ctx.move_to(*dPos)
    self.ctx.show_text(text)

    if 0:
      self.ctx.move_to(dPos[0], dPos[1])
      self.ctx.line_to(dPos[0] + bw, dPos[1])
      self.ctx.line_to(dPos[0] + bw, dPos[1] - bh)
      self.ctx.line_to(dPos[0], dPos[1] - bh)
      self.ctx.line_to(dPos[0], dPos[1])
      self.ctx.close_path()
      self.ctx.stroke()

    return (bw, bh, offset)

  def _addCanvasText2(self, text, pos, font, color=(0, 0, 0), **kwargs):
    if font.weight == 'bold':
      weight = cairo.FONT_WEIGHT_BOLD
    else:
      weight = cairo.FONT_WEIGHT_NORMAL
    self.ctx.select_font_face(font.face, cairo.FONT_SLANT_NORMAL, weight)
    orientation = kwargs.get('orientation', 'E')

    plainText = scriptPattern.sub('', text)

    # for whatever reason, the font size using pango is larger
    # than that w/ default cairo (at least for me)
    pangoCoeff = 0.8

    if have_cairocffi:
      measureLout = pangocairo.pango_cairo_create_layout(self.ctx._pointer)
      pango.pango_layout_set_alignment(measureLout, pango.PANGO_ALIGN_LEFT)
      pango.pango_layout_set_markup(measureLout, plainText.encode('latin1'), -1)
      lout = pangocairo.pango_cairo_create_layout(self.ctx._pointer)
      pango.pango_layout_set_alignment(lout, pango.PANGO_ALIGN_LEFT)
      pango.pango_layout_set_markup(lout, text.encode('latin1'), -1)
      fnt = pango.pango_font_description_new()
      pango.pango_font_description_set_family(fnt, font.face.encode('latin1'))
      pango.pango_font_description_set_size(fnt,
                                            int(round(font.size * pango.PANGO_SCALE * pangoCoeff)))
      pango.pango_layout_set_font_description(lout, fnt)
      pango.pango_layout_set_font_description(measureLout, fnt)
      pango.pango_font_description_free(fnt)
    else:
      cctx = pangocairo.CairoContext(self.ctx)
      measureLout = cctx.create_layout()
      measureLout.set_alignment(pango.ALIGN_LEFT)
      measureLout.set_markup(plainText)
      lout = cctx.create_layout()
      lout.set_alignment(pango.ALIGN_LEFT)
      lout.set_markup(text)
      fnt = pango.FontDescription('%s %d' % (font.face, font.size * pangoCoeff))
      lout.set_font_description(fnt)
      measureLout.set_font_description(fnt)

    # this is a bit kludgy, but empirically we end up with too much
    # vertical padding if we use the text box with super and subscripts
    # for the measurement.
    if have_cairocffi:
      iext = ffi.new('PangoRectangle *')
      lext = ffi.new('PangoRectangle *')
      iext2 = ffi.new('PangoRectangle *')
      lext2 = ffi.new('PangoRectangle *')
      pango.pango_layout_get_pixel_extents(measureLout, iext, lext)
      pango.pango_layout_get_pixel_extents(lout, iext2, lext2)
      w = lext2.width - lext2.x
      h = lext.height - lext.y
    else:
      iext, lext = measureLout.get_pixel_extents()
      iext2, lext2 = lout.get_pixel_extents()
      w = lext2[2] - lext2[0]
      h = lext[3] - lext[1]
    pad = [h * .2, h * .3]
    # another empirical correction: labels draw at the bottom
    # of bonds have too much vertical padding
    if orientation == 'S':
      pad[1] *= 0.5
    bw, bh = w + pad[0], h + pad[1]
    offset = w * pos[2]
    if 0:
      if orientation == 'W':
        dPos = pos[0] - w + offset, pos[1] - h / 2.
      elif orientation == 'E':
        dPos = pos[0] - w / 2 + offset, pos[1] - h / 2.
      else:
        dPos = pos[0] - w / 2 + offset, pos[1] - h / 2.
      self.ctx.move_to(dPos[0], dPos[1])
    else:
      dPos = pos[0] - w / 2. + offset, pos[1] - h / 2.
      self.ctx.move_to(dPos[0], dPos[1])

    self.ctx.set_source_rgb(*color)
    if have_cairocffi:
      pangocairo.pango_cairo_update_layout(self.ctx._pointer, lout)
      pangocairo.pango_cairo_show_layout(self.ctx._pointer, lout)
      gobject.g_object_unref(lout)
      gobject.g_object_unref(measureLout)
    else:
      cctx.update_layout(lout)
      cctx.show_layout(lout)

    if 0:
      self.ctx.move_to(dPos[0], dPos[1])
      self.ctx.line_to(dPos[0] + bw, dPos[1])
      self.ctx.line_to(dPos[0] + bw, dPos[1] + bh)
      self.ctx.line_to(dPos[0], dPos[1] + bh)
      self.ctx.line_to(dPos[0], dPos[1])
      self.ctx.close_path()
      self.ctx.stroke()

    return (bw, bh, offset)

  def addCanvasText(self, text, pos, font, color=(0, 0, 0), **kwargs):
    if have_pango:
      textSize = self._addCanvasText2(text, pos, font, color, **kwargs)
    else:
      textSize = self._addCanvasText1(text, pos, font, color, **kwargs)
    return textSize

  def addCanvasPolygon(self, ps, color=(0, 0, 0), fill=True, stroke=False, **kwargs):
    if not fill and not stroke:
      return
    self.ctx.set_source_rgb(*color)
    self.ctx.move_to(ps[0][0], ps[0][1])
    for p in ps[1:]:
      self.ctx.line_to(p[0], p[1])
    self.ctx.close_path()
    if stroke:
      if fill:
        self.ctx.stroke_preserve()
      else:
        self.ctx.stroke()
    if fill:
      self.ctx.fill()

  def addCanvasDashedWedge(self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs):
    self.ctx.set_line_width(kwargs.get('linewidth', 1))
    self.ctx.set_source_rgb(*color)
    dash = (3, 3)
    pts1 = self._getLinePoints(p1, p2, dash)
    pts2 = self._getLinePoints(p1, p3, dash)

    if len(pts2) < len(pts1):
      pts2, pts1 = pts1, pts2

    for i in range(len(pts1)):
      self.ctx.move_to(pts1[i][0], pts1[i][1])
      self.ctx.line_to(pts2[i][0], pts2[i][1])
    self.ctx.stroke()

  def addCircle(self, center, radius, color=(0, 0, 0), fill=True, stroke=False, alpha=1.0,
                **kwargs):
    if not fill and not stroke:
      return
    self.ctx.set_source_rgba(color[0], color[1], color[2], alpha)
    self.ctx.arc(center[0], center[1], radius, 0, 2. * math.pi)
    self.ctx.close_path()
    if stroke:
      if fill:
        self.ctx.stroke_preserve()
      else:
        self.ctx.stroke()
    if fill:
      self.ctx.fill()
