from sping import colors
from sping.TK import TKCanvas, TKCanvasPIL
from Tkinter import *

# This example program creates a scrolling canvas using TKCanvasPIL as the basis
# The PIL-based canvas may be saved to an image file.
# based on the demo included with python 1.5.2 source under Demo/tkinter/matt/tkcanvas-with-scrollbars.py


class Test(Frame):

  def printit(self):
    print("hi")
    print(self.draw.size)
    print("height = %s" % self.draw.height)

  def saveToJpeg(self):
    print("Saving canvas to file tkCanvasPIL.jpg")
    self.tkpil.save(file="tkCanvasPIL.jpg")

  def createWidgets(self):
    self.question = Label(self, text="Can Find The BLUE Square??????")
    self.question.pack()

    self.QUIT = Button(self, text='QUIT', height=3, command=self.quit)
    self.QUIT.pack(side=BOTTOM)
    Button(self, text='Save to jpeg file', command=self.saveToJpeg).pack(side=BOTTOM)

    spacer = Frame(self, height="0.25i")
    spacer.pack(side=BOTTOM)

    # notice that the scroll region (600 x 600) is larger than
    # displayed size of the widget (400 x 400)

    self.tkpil = TKCanvasPIL(name="SpingTKCanvas", size=(600, 600), master=self,
                             scrollingViewPortSize=(400, 400))

    self.draw = self.tkpil.getTKCanvas()  # retrieve the underlying sping TKCanvas
    # it's a subclass of tk.Canvas

    self.draw.scrollX = Scrollbar(self, orient=HORIZONTAL)
    self.draw.scrollY = Scrollbar(self, orient=VERTICAL)

    # now tie the three together. This is standard boilerplate text
    self.draw['xscrollcommand'] = self.draw.scrollX.set
    self.draw['yscrollcommand'] = self.draw.scrollY.set
    self.draw.scrollX['command'] = self.draw.xview
    self.draw.scrollY['command'] = self.draw.yview

    # draw something into the  canvas Note that the first square
    # is visible, but you need to scroll to see the other ones.

    self.tkpil.drawRect(10, 10, 100, 100, edgeColor=colors.blue, fillColor=colors.green)
    self.tkpil.drawRect(400, 400, 500, 500, edgeColor=colors.blue, fillColor=colors.lightblue)
    self.tkpil.drawRect(30, 400, 130, 500, edgeColor=colors.blue, fillColor=colors.yellow)
    self.tkpil.flush()  # must call this or it won't be visible

    # pack 'em up
    self.draw.scrollX.pack(side=BOTTOM, fill=X)
    self.draw.scrollY.pack(side=RIGHT, fill=Y)
    self.draw.pack(side=LEFT)

  def scrollCanvasX(self, *args):
    print("scrolling", args)
    print(self.draw.scrollX.get())

  def __init__(self, master=None):
    Frame.__init__(self, master)
    Pack.config(self)
    self.createWidgets()


test = Test()

test.mainloop()
