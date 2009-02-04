## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
#  Copyright (C) 2000,2001,2002,2003  greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
""" tools for working with PIL images

"""
#try:
#  from wxPython.wx import *
#except ImportError:
#  hasWx=0
#else:
#  hasWx=1
hasWx=0  
try:
  from qt import *
except ImportError:
  hasQt = 0
else:
  hasQt=1

from PIL import Image
# these are here to help the Installer work:
import PIL.ImageFile
import PIL.GifImagePlugin
import PIL.PngImagePlugin
import PIL.JpegImagePlugin
import PIL.BmpImagePlugin
import PIL.TiffImagePlugin
import PIL.PpmImagePlugin

def ResizeImage(origImg,newSize,filter=Image.BILINEAR,maintainAspect=0,
                priorityX=1):
  """Resizes an image to fit a given space.  The new image is returned.

    **Arguments**

      - origImg: a PIL image

      - newSize: the requested size (either a 2-tuple or a list 2 elements long)

      - filter: the filter to be used in resizing the image

      - maintainAspect: toggles maintaining the aspect ratio of the image

      - prioritiyX: (only meaningful when _maintainAspect_ is nonzero)
           if nonzero, the X size will be given priority in setting the new size,
           otherwise the Y size will take priority

    **Returns**

      a PIL image

    **Notes**

      - if maintainAspect is nonzero, the aspect ratio of the image
       will not be changed.  This implies that the final image may not
       actually be the requested size.

  """
  
  if maintainAspect:
    if priorityX:
      scaleFact = float(origImg.size[0])/newSize[0]
    else:
      scaleFact = float(origImg.size[1])/newSize[1]

    newSize = (int(origImg.size[0]*scaleFact),int(origImg.size[1].scaleFact))

  newImg = origImg.resize(newSize,filter)
  return newImg
    
def FitImage(origImg,newSize,filter=Image.BILINEAR,bgColor=(255,255,255)):
  """Fits an image into a box of a particular size.

    **Arguments**

      - origImg: a PIL image

      - newSize: the requested size (either a 2-tuple or a list 2 elements long)

      - filter: the filter to be used in resizing the image

      - bgColor: the background color to start with
      
    **Returns**

      a PIL image

    **Notes**

      - there may be blank spaces around the original image in the new image,
        these will be filled with _bgColor_
        
  """
  tmpImg = origImg.convert('RGB')
  newImg = Image.new(tmpImg.mode,newSize,bgColor)

  scaleFact = min(float(newSize[0])/origImg.size[0],
                  float(newSize[1])/origImg.size[1])
  
  if scaleFact < 1:
    tImg = origImg.resize((int(origImg.size[0]*scaleFact),int(origImg.size[1]*scaleFact)),filter)
  else:
    tImg = origImg
    

  xDiff = newSize[0] - tImg.size[0]
  if xDiff > 0:
    xLoc = xDiff/2
  else:
    xLoc = 0

  yDiff = newSize[1] - tImg.size[1]
  if yDiff > 0:
    yLoc = yDiff/2
  else:
    yLoc = 0
  
  newImg.paste(tImg,(xLoc,yLoc))
  return newImg.convert(origImg.mode)

def NumericMatrixToImage(data,scaleCols=0,transposeIt=0,
                         minColor=(0,0,0),maxColor=(255,255,255)):
  import numpy
  # copy the data
  data = numpy.array(data,numpy.float)
  if numpy.transpose:
    data = numpy.transpose(data)
  #nRows,nCols = data.shape
  nRows,nCols = data.shape
  if scaleCols:
    minIndices = numpy.argmin(data)
    maxIndices = numpy.argmax(data)
    mins = zeros(nCols)
    maxs = zeros(nCols)
    for i in range(nCols):
      mins[i] = data[i][minIndices[i]]
      maxs[i] = data[i][maxIndices[i]]
    # subtract off the minimum
    data -= mins
    maxs -= mins
    # no zeros here please, we're dividing:
    maxs += numpy.equal(maxs,0.0)

    # and divide:
    data /= maxs
  # okey dokey, get a three D matrix:
  imgMat = numpy.ones((nRows,nCols,3),numpy.integer)

  # start at minColor:
  minColor = numpy.array(minColor)
  maxColor = numpy.array(maxColor)
  imgMat *= minColor
  deltaColor = maxColor-minColor
  # and move to maxColor:
  for i in range(nRows):
    for j in range(nCols):
      imgMat[i,j] += (deltaColor*data[i,j]).astype(numpy.integer)
  d = imgMat.astype('B').tostring()
  img = Image.fromstring('RGB',(nCols,nRows),d)
  return img

if hasWx:
  def PilImgToWxBmp(pilImg):
    """ converts a PIL image into a wxPython bitmap

      **Arguments**

        - pilImg: a PIL image

      **Returns**

        a wxPython bitmap

    """
    wxImg = wxEmptyImage(pilImg.size[0],pilImg.size[1])
    wxImg.SetData(pilImg.tostring())
    bmp = wxImg.ConvertToBitmap()
    return bmp

if hasQt:
  def PilImgToQPixmap(pilImg):
    from StringIO import StringIO
    sio = StringIO()
    pilImg.save(sio,format='png')
    pm = QPixmap()
    pm.loadFromData(sio.getvalue())
    return pm

if __name__ == '__main__':
  if 0:
    boxImg = Image.open('12h.gif').convert('RGB')
    img1 = ResizeImage(boxImg,(100,200))
    img1.save('test1.gif')

    img2 = FitImage(boxImg,(100,200),bgColor=(200,200,200))
    img2.save('test2.gif')

    img3 = FitImage(boxImg,(50,200),bgColor=(200,200,200))
    img3.save('test3.gif')
  else:
    vs = array([[1.,.5,0,0,0,0],[.5,1,.5,1,0,1],[0,.5,1,0,0,1]])
    img = NumericMatrixToImage(vs)
    img = img.resize((200,200))
    img.save('foo.gif')

