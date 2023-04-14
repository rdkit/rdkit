from sping import colors
from sping.TK import TKCanvas
from Tkinter import *

# This example program creates a scrolling canvas with sping.TK
# how to tie scrollbars and canvases together. The mechanism

# based on the demo included with python 1.5.2 source under Demo/tkinter/matt/tkcanvas-with-scrollbars.py


class Test(Frame):

  def printit(self):
    print("hi")
    print(self.draw.size)
    print("height = %s" % self.draw.height)

  def saveToPostscript(self):
    "Save the whole canvas to a postscript file"
    print("Saving output to file tkCanvas.eps using sping.PS")
    import sping.PS
    can = sping.PS.PSCanvas(size=self.draw.size, name="tkCanvas.eps")
    self.spingDrawingCommands(can)
    can.save("tkCanvas.eps")

  def createWidgets(self):
    self.question = Label(self, text="Can Find The BLUE Square??????")
    self.question.pack()

    self.QUIT = Button(self, text='QUIT', background='red', height=3, command=self.quit)
    self.QUIT.pack(side=BOTTOM, fill=BOTH)

    Button(self, text='Save To PS', command=self.saveToPostscript).pack()

    spacer = Frame(self, height="0.25i")
    spacer.pack(side=BOTTOM)

    self.draw = TKCanvas(name="SpingTKCanvas", size=(600, 600), master=self,
                         scrollingViewPortSize=(400, 400))

    self.draw.scrollX = Scrollbar(self, orient=HORIZONTAL)
    self.draw.scrollY = Scrollbar(self, orient=VERTICAL)

    # now tie the three together. This is standard boilerplate text
    self.draw['xscrollcommand'] = self.draw.scrollX.set
    self.draw['yscrollcommand'] = self.draw.scrollY.set
    self.draw.scrollX['command'] = self.draw.xview
    self.draw.scrollY['command'] = self.draw.yview

    # draw something. Note that the first square
    # is visible, but you need to scroll to see the second one.

    self.spingDrawingCommands(self.draw)  # hand this a sping Canvas and it draws to it

    # pack 'em up
    self.draw.scrollX.pack(side=BOTTOM, fill=X)
    self.draw.scrollY.pack(side=RIGHT, fill=Y)
    self.draw.pack(side=LEFT)

  def scrollCanvasX(self, *args):
    print("scrolling", args)
    print(self.draw.scrollX.get())

  def spingDrawingCommands(self, anySpingCanvas):
    anySpingCanvas.drawRect(10, 10, 100, 100, edgeColor=colors.blue, fillColor=colors.green)
    anySpingCanvas.drawRect(400, 400, 500, 500, edgeColor=colors.blue, fillColor=colors.lightblue)
    anySpingCanvas.drawRect(30, 400, 130, 500, edgeColor=colors.blue, fillColor=colors.yellow)

  def __init__(self, master=None):
    Frame.__init__(self, master)
    Pack.config(self)
    self.createWidgets()


test = Test()

test.mainloop()
