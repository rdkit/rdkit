#
#  Copyright (C) 2001-2004  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#
""" The "parser" for compound descriptors.

I almost hesitate to document this, because it's not the prettiest
thing the world has ever seen... but it does work (for at least some
definitions of the word).

Rather than getting into the whole mess of writing a parser for the
compound descriptor expressions, I'm just using string substitutions
and python's wonderful ability to *eval* code.

It would probably be a good idea at some point to replace this with a
real parser, if only for the flexibility and intelligent error
messages that would become possible.

The general idea is that we're going to deal with expressions where
atomic descriptors have some kind of method applied to them which
reduces them to a single number for the entire composition.  Compound
descriptors (those applicable to the compound as a whole) are not
operated on by anything in particular (except for standard math stuff).

Here's the general flow of things:

  1) Composition descriptor references ($a, $b, etc.) are replaced with the
     corresponding descriptor names using string substitution.
     (*_SubForCompoundDescriptors*)

  2) Atomic descriptor references ($1, $2, etc) are replaced with lookups
     into the atomic dict with "DEADBEEF" in place of the atom name.
     (*_SubForAtomicVars*)

  3) Calls to Calculator Functions are augmented with a reference to
     the composition and atomic dictionary
     (*_SubMethodArgs*)

**NOTE:**

  anytime we don't know the answer for a descriptor, rather than
  throwing a (completely incomprehensible) exception, we just return
  -666.  So bad descriptor values should stand out like sore thumbs.

"""

# The wildcard import is required to make functions available for the eval statement
from math import *

from rdkit import RDConfig

__DEBUG = False

# we do this to allow the use of stuff in the math module

# ----------------------
# atomic descriptor section
# ----------------------
# these are the methods which can be applied to ATOMIC descriptors.
knownMethods = ['SUM', 'MIN', 'MAX', 'MEAN', 'AVG', 'DEV', 'HAS']


def HAS(strArg, composList, atomDict):
  """ *Calculator Method*

    does a string search

    **Arguments**

      - strArg: the arguments in string form

      - composList: the composition vector

      - atomDict: the atomic dictionary

    **Returns**

      1 or 0

  """
  splitArgs = strArg.split(',')
  if len(splitArgs) > 1:
    for atom, _ in composList:
      tStr = splitArgs[0].replace('DEADBEEF', atom)
      where = eval(tStr)
      what = eval(splitArgs[1])
      if what in where:
        return 1
    return 0
  else:
    return -666


def SUM(strArg, composList, atomDict):
  """ *Calculator Method*

    calculates the sum of a descriptor across a composition

    **Arguments**

      - strArg: the arguments in string form

      - compos: the composition vector

      - atomDict: the atomic dictionary

    **Returns**

      a float

  """
  accum = 0.0
  for atom, num in composList:
    tStr = strArg.replace('DEADBEEF', atom)
    accum = accum + eval(tStr) * num
  return accum


def MEAN(strArg, composList, atomDict):
  """ *Calculator Method*

    calculates the average of a descriptor across a composition

    **Arguments**

      - strArg: the arguments in string form

      - compos: the composition vector

      - atomDict: the atomic dictionary

    **Returns**

      a float

  """
  accum = 0.0
  nSoFar = 0
  for atom, num in composList:
    tStr = strArg.replace('DEADBEEF', atom)
    accum = accum + eval(tStr) * num
    nSoFar = nSoFar + num
  return accum / nSoFar


AVG = MEAN


def DEV(strArg, composList, atomDict):
  """ *Calculator Method*

    calculates the average deviation of a descriptor across a composition

    **Arguments**

      - strArg: the arguments in string form

      - compos: the composition vector

      - atomDict: the atomic dictionary

    **Returns**

      a float

  """
  avg = MEAN(strArg, composList, atomDict)
  accum = 0.0
  nSoFar = 0.0
  for atom, num in composList:
    tStr = strArg.replace('DEADBEEF', atom)
    accum = accum + abs(eval(tStr) - avg) * num
    nSoFar = nSoFar + num
  return accum / nSoFar


def MIN(strArg, composList, atomDict):
  """ *Calculator Method*

    calculates the minimum value of a descriptor across a composition

    **Arguments**

      - strArg: the arguments in string form

      - compos: the composition vector

      - atomDict: the atomic dictionary

    **Returns**

      a float

  """
  accum = []
  for atom, _ in composList:
    tStr = strArg.replace('DEADBEEF', atom)
    accum.append(eval(tStr))
  return min(accum)


def MAX(strArg, composList, atomDict):
  """ *Calculator Method*

    calculates the maximum value of a descriptor across a composition

    **Arguments**

      - strArg: the arguments in string form

      - compos: the composition vector

      - atomDict: the atomic dictionary

    **Returns**

      a float

  """
  accum = []
  for atom, _ in composList:
    tStr = strArg.replace('DEADBEEF', atom)
    accum.append(eval(tStr))
  return max(accum)


# ------------------
#  string replacement routines
#    these are not intended to be called by clients
# ------------------


def _SubForAtomicVars(cExpr, varList, dictName):
  """ replace atomic variables with the appropriate dictionary lookup

   *Not intended for client use*

  """
  for i in range(len(varList)):
    cExpr = cExpr.replace('$%d' % (i + 1), '%s["DEADBEEF"]["%s"]' % (dictName, varList[i]))
  return cExpr


def _SubForCompoundDescriptors(cExpr, varList, dictName):
  """ replace compound variables with the appropriate list index

   *Not intended for client use*

  """
  for i in range(len(varList)):
    cExpr = cExpr.replace('$%s' % chr(ord('a') + i), '%s["%s"]' % (dictName, varList[i]))
  return cExpr


def _SubMethodArgs(cExpr, knownMethods):
  """ alters the arguments of calls to calculator methods

  *Not intended for client use*

  This is kind of putrid (and the code ain't so pretty either)
  The general idea is that the various special methods for atomic
  descriptors need two extra arguments (the composition and the atomic
  dict).  Rather than make the user type those in, we just find
  invocations of these methods and fill out the function calls using
  string replacements.
  """
  res = cExpr
  for method in knownMethods:
    p = 0
    while p != -1 and p < len(res):
      p = res.find(method, p)
      if p != -1:
        p = p + len(method) + 1
        start = p
        parenCount = 1
        while parenCount and p < len(res):
          if res[p] == ')':
            parenCount = parenCount - 1
          elif res[p] == '(':
            parenCount = parenCount + 1
          p = p + 1
        if p <= len(res):
          res = res[0:start] + "'%s',compos,atomDict" % (res[start:p - 1]) + res[p - 1:]
  return res


def CalcSingleCompoundDescriptor(compos, argVect, atomDict, propDict):
  """ calculates the value of the descriptor for a single compound

    **ARGUMENTS:**

      - compos: a vector/tuple containing the composition
         information... in the form:
         '[("Fe",1.),("Pt",2.),("Rh",0.02)]'

      - argVect: a vector/tuple with three elements:

           1) AtomicDescriptorNames:  a list/tuple of the names of the
             atomic descriptors being used. These determine the
             meaning of $1, $2, etc. in the expression

           2) CompoundDescriptorNames:  a list/tuple of the names of the
             compound descriptors being used. These determine the
             meaning of $a, $b, etc. in the expression

           3) Expr: a string containing the expression to be used to
             evaluate the final result.

      - atomDict:
           a dictionary of atomic descriptors.  Each atomic entry is
           another dictionary containing the individual descriptors
           and their values

      - propVect:
           a list of descriptors for the composition.

    **RETURNS:**

      the value of the descriptor, -666 if a problem was encountered

    **NOTE:**

      - because it takes rather a lot of work to get everything set
          up to calculate a descriptor, if you are calculating the
          same descriptor for multiple compounds, you probably want to
          be calling _CalcMultipleCompoundsDescriptor()_.

  """
  try:
    atomVarNames = argVect[0]
    compositionVarNames = argVect[1]
    formula = argVect[2]
    formula = _SubForCompoundDescriptors(formula, compositionVarNames, 'propDict')
    formula = _SubForAtomicVars(formula, atomVarNames, 'atomDict')
    evalTarget = _SubMethodArgs(formula, knownMethods)
  except Exception:
    if __DEBUG:
      import traceback
      print('Sub Failure!')
      traceback.print_exc()
      print(evalTarget)
      print(propDict)
      raise RuntimeError('Failure 1')
    else:
      return -666

  try:
    v = eval(evalTarget)
  except Exception:
    if __DEBUG:
      import traceback
      outF = open(RDConfig.RDCodeDir + '/ml/descriptors/log.txt', 'a+')
      outF.write('#------------------------------\n')
      outF.write('formula: %s\n' % repr(formula))
      outF.write('target: %s\n' % repr(evalTarget))
      outF.write('propDict: %s\n' % (repr(propDict)))

      outF.write('keys: %s\n' % (repr(sorted(atomDict))))
      outF.close()
      print('ick!')
      print('formula:', formula)
      print('target:', evalTarget)
      print('propDict:', propDict)
      print('keys:', atomDict.keys())
      traceback.print_exc()
      raise RuntimeError('Failure 2')
    else:
      v = -666
  return v


def CalcMultipleCompoundsDescriptor(composVect, argVect, atomDict, propDictList):
  """ calculates the value of the descriptor for a list of compounds

    **ARGUMENTS:**

      - composVect: a vector of vector/tuple containing the composition
         information.
         See _CalcSingleCompoundDescriptor()_ for an explanation of the elements.

      - argVect: a vector/tuple with three elements:

           1) AtomicDescriptorNames:  a list/tuple of the names of the
             atomic descriptors being used. These determine the
             meaning of $1, $2, etc. in the expression

           2) CompoundDsscriptorNames:  a list/tuple of the names of the
             compound descriptors being used. These determine the
             meaning of $a, $b, etc. in the expression

           3) Expr: a string containing the expression to be used to
             evaluate the final result.

      - atomDict:
           a dictionary of atomic descriptors.  Each atomic entry is
           another dictionary containing the individual descriptors
           and their values

      - propVectList:
         a vector of vectors of descriptors for the composition.

    **RETURNS:**

      a vector containing the values of the descriptor for each
      compound.  Any given entry will be -666 if problems were
      encountered

  """
  res = [-666] * len(composVect)
  try:
    atomVarNames = argVect[0]
    compositionVarNames = argVect[1]
    formula = argVect[2]
    formula = _SubForCompoundDescriptors(formula, compositionVarNames, 'propDict')
    formula = _SubForAtomicVars(formula, atomVarNames, 'atomDict')
    evalTarget = _SubMethodArgs(formula, knownMethods)
  except Exception:
    return res
  for i in range(len(composVect)):
    propDict = propDictList[i]
    compos = composVect[i]
    try:
      v = eval(evalTarget)
    except Exception:
      v = -666
    res[i] = v
  return res


# ------------
#  Demo/testing code
# ------------
def _exampleCode():  # pragma: nocover
  piece1 = [['d1', 'd2', 's1'], ['d1', 'd2', 's1']]
  aDict = {'Fe': {'d1': 1., 'd2': 2., 's1': 'abc'}, 'Pt': {'d1': 10., 'd2': 20., 's1': 'def'}}
  pDict = {'d1': 100., 'd2': 200.}
  compos = [('Fe', 1), ('Pt', 1)]

  cExprs = [
    "SUM($1)", "SUM($1)+SUM($2)", "SUM($1)+SUM($1)", "MEAN($1)", "DEV($2)", "MAX($1)",
    "MIN($1)/MAX($1)", "MIN($2)", "SUM($1)/$a", "sqrt($a+$b)", "SUM((3.*$1)/($2))", 'HAS($3,"def")',
    'HAS($3,"xyz")', "foo"
  ]

  for cExpr in cExprs:
    argVect = piece1 + [cExpr]
    print(cExpr)
    print(CalcSingleCompoundDescriptor(compos, argVect, aDict, pDict))
    print(CalcMultipleCompoundsDescriptor([compos, compos], argVect, aDict, [pDict, pDict]))


if __name__ == '__main__':  # pragma: nocover
  _exampleCode()
