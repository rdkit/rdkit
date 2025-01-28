#
#  Copyright (C) 2000  greg Landrum
#
""" functionality for generating an image showing the results of a composite model
voting on a data set

  Uses *Numeric* and *PIL*

"""

import numpy
from PIL import Image, ImageDraw


def CollectVotes(composite, data, badOnly):
  """ collects the votes from _composite_ for the examples in _data_

    **Arguments**

      - composite: a composite model

      - data: a list of examples to run through _composite_

      - badOnly: if set only bad (misclassified) examples will be kept

    **Returns**

      a 4-tuple containing:

        1) the expanded list of vote details (see below)

        2) the list of predicted results

        3) the list of true results

        4) the number of miscounted examples


    **Notes**

pp      - the expanded list of vote details consists of:

        '[ vote1, vote2, ... voteN, 0, res, trueRes]'

        where _res_ is the predicted results and _trueRes_ is the actual result.
        The extra zero is included to allow a line to be drawn between the votes
        and the results.

  """
  res = []
  values = []
  trueValues = []
  misCount = 0
  for pt in data:
    val, _ = composite.ClassifyExample(pt)
    predict = pt[-1]
    if not badOnly or val != predict:
      values.append(val)
      trueValues.append(predict)
      if val != predict:
        misCount = misCount + 1
      res.append(composite.GetVoteDetails() + [0, val, pt[-1]])
  return res, values, trueValues, misCount


def BuildVoteImage(nModels, data, values, trueValues=[], sortTrueVals=0, xScale=10, yScale=2,
                   addLine=1):
  """ constructs the actual image

    **Arguments**

      - nModels: the number of models in the composite

      - data: the results of voting

      - values: predicted values for each example

      - trueValues: true values for each example

      - sortTrueVals: if nonzero the votes will be sorted so
        that the _trueValues_ are in order, otherwise the sort
        is by _values_

      - xScale: number of pixels per vote in the x direction

      - yScale: number of pixels per example in the y direction

      - addLine: if nonzero, a purple line is drawn separating
         the votes from the examples

    **Returns**

      a PIL image

  """
  nData = len(data)
  data = numpy.array(data, numpy.integer)
  if sortTrueVals and trueValues != []:
    order = numpy.argsort(trueValues)
  else:
    order = numpy.argsort(values)
  data = [data[x] for x in order]
  maxVal = max(numpy.ravel(data))
  data = data * 255 / maxVal
  datab = data.astype('B')
  img = Image.frombytes('L', (nModels, nData), datab.tobytes())

  if addLine:
    img = img.convert('RGB')
    canvas = ImageDraw.Draw(img)
    if trueValues != []:
      canvas.line([(nModels - 3, 0), (nModels - 3, nData)], fill=(128, 0, 128))
    else:
      canvas.line([(nModels - 2, 0), (nModels - 2, nData)], fill=(128, 0, 128))
  img = img.resize((nModels * xScale, nData * yScale))
  return img


def VoteAndBuildImage(composite, data, badOnly=0, sortTrueVals=0, xScale=10, yScale=2, addLine=1):
  """ collects votes on the examples and constructs an image

    **Arguments**

      - composte: a composite model

      - data: the examples to be voted upon

      - badOnly: if nonzero only the incorrect votes will be shown

      - sortTrueVals: if nonzero the votes will be sorted so
        that the _trueValues_ are in order, otherwise the sort
        is by _values_

      - xScale: number of pixels per vote in the x direction

      - yScale: number of pixels per example in the y direction

      - addLine: if nonzero, a purple line is drawn separating
         the votes from the examples

    **Returns**

      a PIL image

  """
  nModels = len(composite) + 3
  print('nModels:', nModels - 3)

  res, values, trueValues, misCount = CollectVotes(composite, data, badOnly)
  print('%d examples were misclassified' % misCount)
  img = BuildVoteImage(nModels, res, values, trueValues, sortTrueVals, xScale, yScale, addLine)
  return img


def Usage():
  """ provides a list of arguments for when this is used from the command line

  """
  import sys

  print('Usage: VoteImg.py [optional arguments] <modelfile.pkl> <datafile.qdat>')
  print('Optional Arguments:')
  print('\t-o outfilename: the name of the output image file.')
  print('\t                The extension determines the type of image saved.')
  print('\t-b: only include bad (misclassified) examples')
  print('\t-t: sort the results by the true (input) classification')
  print('\t-x scale: scale the image along the x axis (default: 10)')
  print('\t-y scale: scale the image along the y axis (default: 2)')
  print('\t-d databasename: instead of using a qdat file, pull the data from')
  print('\t                 a database.  In this case the filename argument')
  print('\t                 is used to indicate the name of the table in the database.')

  sys.exit(-1)


if __name__ == '__main__':
  import getopt
  import pickle
  import sys

  from rdkit.ML.Data import DataUtils

  args, extra = getopt.getopt(sys.argv[1:], 'o:bthx:y:d:')
  if len(extra) < 2:
    Usage()
  badOnly = 0
  sortTrueVals = 0
  xScale = 10
  yScale = 2
  dbName = ''
  outFileName = 'foo.png'
  for arg, val in args:
    if arg == '-b':
      badOnly = 1
    elif arg == '-t':
      sortTrueVals = 1
    elif arg == '-o':
      outFileName = val
    elif arg == '-x':
      xScale = int(val)
    elif arg == '-y':
      yScale = int(val)
    elif arg == '-d':
      dbName = val
    elif arg == '-h':
      Usage()
    else:
      Usage()
  modelFile = open(extra[0], 'rb')
  model = pickle.load(modelFile)

  fName = extra[1]
  if dbName == '':
    data = DataUtils.BuildQuantDataSet(fName)
  else:
    data = DataUtils.DBToQuantData(dbName, fName)  # Function no longer defined

  dataSet = data.GetNamedData()

  img = VoteAndBuildImage(model, dataSet, badOnly=badOnly, sortTrueVals=sortTrueVals, xScale=xScale,
                          yScale=yScale)
  img.save(outFileName)
