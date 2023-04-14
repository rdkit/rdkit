# pidWX.py
# base on piddleWX.py -- a wxPython backend for PIDDLE
# Copyright (c) 2000  Paul & Kevin Jacobs
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2, or (at your option)
#   any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
"""WXCanvas

This class implements a PIDDLE Canvas object that draws using wxPython, the
wxWindows Python language bindings into a GUI window.
"""
# Status as of version 1.0:
#
# -  All test cases work correctly!  (except for the symbol font)
#
# -  The full interactive canvas API is supported, but not extensively
#    tested.  However, the default status bar object displays interactive
#    events nicely.
#
# -  The drawing code has been separated into a PiddleWxDc class that
#    enables any wxPython application to use Piddle methods as a
#    DeviceContext.
#
# History:
#
#    0.1     First release.  Most PIDDLE functionality is there.
#    0.5     PIL support and code factoring added by Jeffrey Kunce.
#    1.0     All PIDDLE functionality supported.  Rotated test, multiline text,
#            and pil image support added.
#
# Known Bugs:
#
# -  Due to limitations in wxWindows, some fonts will not display when rotated;
#    As the default GUI font is not true-type, it cannot be rotated, and
#    PiddleWX will automatically re-assign the offending font to a type face
#    that can.
#
#
# Extensions:
#
# -  The canvas init function can take up to two additional arguments, both
#    of which are on by default.  The first is a boolean which enables the
#    interactive mode when true.  The second is a boolean that will add a
#    status bar to the canvas when true.
#
# -  The canvas can handle a left button released mouse event.  It's
#    handler is called onClickUp and is similar to the other event handlers.
#    handlers.  The standard onClick event registers left-click button down
#    events.
#
# -  The canvas can switch between interactive and non-interactive modes using
#    the SetInteractive function
#
# -  The canvas also supports an onLeaveWindow event, to trap when the mouse
#    leaves the canvas window, but due to a wxWindows bug it may not work.
#
# Limitations:
#
# -  There is no symbol font support at the moment.  It should be possible
#    under MS Windows.  Not sure about UNIX (wxGTK).
#
# -  The default status bar cursor co-ordinates display may not behave
#    properly due to several bugs.
#
# TODO:
#
# -  Make the status bar more user-configurable.
#
# -  Support printing & generation of Postscript output (just because it is
#    easy to do).

from pidWxDc import PiddleWxDc
from wxPython.wx import *

__version__ = "1.0"
__date__ = "February 6, 2000"
__author__ = "Paul & Kevin Jacobs"


class _WXCanvasDefaultStatusBar(wxStatusBar):
  """This status bar displays clear and quit buttons, as well as the
   location of the mouse and if the left button is down."""

  def __init__(self, canvas, pos, size):

    # These are terribly wrong here in the init under wxGTK but they firm up after a while
    # I think they need to be set in OnSize event

    #wxStatusBar.__init__(self, canvas.window, -1) # known to work w/ wxGTK

    wxStatusBar.__init__(self, canvas.window, -1, pos, size)  # original
    self.parentwindow = canvas.window
    self.parentwindow.SetStatusBar(self)  # added because it seems to be necessary

    self.SetFieldsCount(3)

    self.sizeChanged = false
    self.text = ""
    self.old_text = ""

    self.canvas = canvas

    qID = NewId()
    cID = NewId()
    sID = NewId()
    bID = NewId()

    # under linux this field value is wrong at this point in constructor so use OnSize event as well
    # Note: under windows the field sizes are correct at this point

    self.quitButton = wxButton(self, qID, "Quit")
    self.clearButton = wxButton(self, cID, "Clear")
    self.click = wxCheckBox(self, bID,
                            "Click")  # This checkbox is down when the left button is pressed

    self.Reposition()  # set up sizes for layout

    EVT_BUTTON(self, qID, canvas._OnQuit)
    EVT_BUTTON(self, cID, canvas._OnClear)
    EVT_PAINT(self, self.repaint)
    EVT_SIZE(self, self.OnSize)
    EVT_IDLE(self, self.OnIdle)

  def repaint(self, event):
    dc = wxPaintDC(self)
    self.draw(dc)

  def redraw(self):
    dc = wxClientDC(self)
    self.draw(dc)

  def draw(self, dc):
    #
    field = self.GetFieldRect(1)
    extents = dc.GetTextExtent(self.old_text)
    dc.SetPen(wxTRANSPARENT_PEN)
    dc.SetBrush(dc.GetBackground())
    field = self.GetFieldRect(2)
    dc.DrawRectangle(field.x, field.y, field.width, field.height)
    dc.SetFont(wxFont(9, wxDEFAULT, wxNORMAL, wxNORMAL))
    extents = dc.GetTextExtent(self.text)
    dc.DrawText(self.text, field.x + field.width / 2 - extents[0] / 2, field.y)
    self.old_text = self.text

  def SetStatusText(self, s):
    self.text = s
    self.redraw()

  def OnOver(self, x, y):
    self.text = repr(x) + "," + repr(y)
    self.redraw()

  def OnClick(self, x, y):
    self.text = repr(x) + "," + repr(y)
    self.click.SetValue(true)
    self.redraw()

  def OnLeaveWindow(self):
    self.click.SetValue(false)
    self.text = ""
    self.redraw()

  def OnClickUp(self, x, y):
    self.click.SetValue(false)

  def OnSize(self, evt):
    self.Reposition()

    # Set a flag so the idle time handler will also do the repositioning.
    # It is done this way to get around a buglet where GetFieldRect is not
    # accurate during the EVT_SIZE resulting from a frame maximize.
    self.sizeChanged = true

  def OnIdle(self, evt):
    if self.sizeChanged:
      self.Reposition()

  def Reposition(self):
    #print "Field Rects:"
    #print self.GetFieldRect(0)
    #print self.GetFieldRect(1)
    #print self.GetFieldRect(2)

    # layout field 0 with Buttons for quit and clear

    field = self.GetFieldRect(0)
    self.quitButton.SetPosition(wxPoint(field.x, field.y))
    self.quitButton.SetSize(wxSize(field.width / 2, field.height))

    self.clearButton.SetPosition(wxPoint(field.x + field.width / 2, field.y))
    self.clearButton.SetSize(wxSize(field.width / 2, field.height))

    # layout sizing of field 1 w/ check box

    field = self.GetFieldRect(1)
    self.click.SetPosition(wxPoint(field.x + field.width / 2 - 20, field.y))
    self.click.SetSize(wxSize(100, field.height))

    self.sizeChanged = false


############################################################################


class WXCanvas(PiddleWxDc):

  def __init__(self, size=(300, 300), name="piddleWX", status_bar=None, interactive=1,
               show_status=1):
    """Works like all other PIDDLE canvases, except with extra interactive controls.
       interactive is set if the canvas is to use the interactive parts of
       the PIDDLE API, and show_status controls if the default status bar is
       shown"""

    window = wxFrame(NULL, -1, "WXCanvas", wxPyDefaultPosition, wxSize(size[0], size[1]))
    window.Show(true)

    self.window = window

    # Resize the window so the client area is equal to the canvas size
    CSize = window.GetClientSizeTuple()

    # Only leave space for a status bar if it is specified
    if show_status:
      # (should eventually be a call to getStatusBarHeight())
      status_area = 20
    else:
      status_area = 0

    window.SetSize(
      wxSize(size[0] + (size[0] - CSize[0]), size[1] + (size[1] - CSize[1] + status_area)))

    # This bitmap is used to buffer drawing commands.  It is set to the same
    # depth as the screen - explicitly changing it can cause errors when it
    # is blitted to the screen

    bmp = wxEmptyBitmap(size[0], size[1])
    MemDc = wxMemoryDC()
    MemDc.SelectObject(bmp)

    MemDc.Clear()

    self.MemDc = MemDc
    PiddleWxDc.__init__(self, MemDc, size, name)

    self.window = window
    self.size = size

    # Different status bars can be substituted in by overriding self.sb
    if status_bar is None:
      self.sb = _WXCanvasDefaultStatusBar(self, wxPoint(0, size[1]), wxSize(size[0], 20))
    else:
      self.sb = status_bar

    if show_status == 0:
      self.sb.Show(false)

    self.sb.redraw()

    # The default behavior for ignoring events is to pass it only to the
    # status bar.

    # onClick: x,y is Canvas coordinates of mouseclick def
    def ignoreClick(canvas, x, y):
      canvas.sb.OnClick(x, y)

    self.onClick = ignoreClick

    # onOver: x,y is Canvas location of mouse
    def ignoreOver(canvas, x, y):
      canvas.sb.OnOver(x, y)

    self.onOver = ignoreOver

    # onKey: key is printable character or one of the constants above;
    #    modifiers is a tuple containing any of (modshift, modControl)
    def ignoreKey(canvas, key, modifiers):
      pass

    self.onKey = ignoreKey

    # onClickUp: This is an extension:  It registers mouse left-button up
    # events
    def ignoreClickUp(canvas, x, y):
      canvas.sb.OnClickUp(x, y)

    self.onClickUp = ignoreClickUp

    self.interactive = interactive

    EVT_PAINT(window, self._OnPaint)
    EVT_LEFT_DOWN(window, self._OnClick)
    EVT_LEFT_UP(window, self._OnClickUp)
    EVT_MOTION(window, self._OnOver)
    EVT_CHAR(window, self._OnKey)
    # Does not seem to be generated as of wxPython 2.1.13
    EVT_LEAVE_WINDOW(window, self._OnLeaveWindow)

    # onLeaveWindow: This is an extension; it is called when the mouse
    # leaves the canvas window
    def leaveWindow(canvas):
      canvas.sb.OnLeaveWindow()

    self.onLeaveWindow = leaveWindow

  ################################################################
  #  Event Managers for wxPython.  To override event handling, alter the
  #  PIDDLE event handlers, not these

  def _OnClick(self, event):
    if self.interactive == false:
      return None
    if event.GetY() <= self.size[1]:
      self.onClick(self, event.GetX(), event.GetY())

  def _OnClickUp(self, event):
    if self.interactive == false:
      return None
    self.onClickUp(self, event.GetX(), event.GetY())

  def _OnOver(self, event):
    if self.interactive == false:
      return None
    if event.GetY() <= self.size[1]:
      self.onOver(self, event.GetX(), event.GetY())

  def _OnLeaveWindow(self, event):
    if self.interactive == false:
      return None
    self.onLeaveWindow(self)

  def _OnKey(self, event):
    # The following logic is OK for ASCII characters.  Others will need work.
    code = event.KeyCode()
    key = None
    if code >= 0 and code < 256:
      key = chr(event.KeyCode())

    modifier = []

    if event.ControlDown():
      modifier.append('modControl')

    if event.ShiftDown():
      modifier.append('modshift')

    self.onKey(self, key, tuple(modifier))

  def _OnPaint(self, event):
    dc = wxPaintDC(self.window)
    dc.Blit(0, 0, self.size[0], self.size[1], self.MemDc, 0, 0, wxCOPY)
    del dc

  def _OnQuit(self, event):
    """Closes the canvas.  Call to return control your application"""
    self.window.Close()

  def _OnClear(self, event):
    """Clears the canvas by emptying the memory buffer, and redrawing"""
    self.MemDc.Clear()
    dc = wxClientDC(self.window)
    dc.Blit(0, 0, self.size[0], self.size[1], self.MemDc, 0, 0, wxCOPY)

  #############################################################

  #------------ canvas capabilities -------------
  def isInteractive(self):
    """Returns 1 if onClick and onOver events are possible, 0 otherwise."""
    return self.interactive

  def canUpdate(self):
    """Returns 1 if the drawing can be meaningfully updated over time (e.g.,
       screen graphics), 0 otherwise (e.g., drawing to a file)."""
    return 1

  #------------ general management -------------
  def clear(self):
    self.Clear()
    dc = wxClientDC(self.window)
    dc.Blit(0, 0, self.size[0], self.size[1], self.MemDc, 0, 0, wxCOPY)

  def flush(self):
    """Copies the contents of the memory buffer to the screen and enters the
       application main loop"""

    dc = wxClientDC(self.window)
    dc.Blit(0, 0, self.size[0], self.size[1], self.MemDc, 0, 0, wxCOPY)
    del dc

  def setInfoLine(self, s):
    """For interactive Canvases, displays the given string in the 'info
       line' somewhere where the user can probably see it."""
    if self.sb is not None:
      self.sb.SetStatusText(s)
