#
#  Copyright (C) 2001,2002  greg Landrum and Rational Discovery LLC
#
""" descriptor calculator for compounds defined by a composition alone
  (only the composition is required)

"""

from rdkit import RDConfig
from rdkit.ML.Descriptors import Descriptors, Parser
from rdkit.utils import chemutils

# the list of possible ways to count valence electrons that we know
countOptions = [('NVAL', 'total number of valence electrons'),
                ('NVAL_NO_FULL_F', 'number of valence electrons neglecting filled f shells'),
                ('NVAL_NO_FULL_D', 'number of valence electrons neglecting filled d shells'),
                ('NVAL_NO_FULL', 'number of valence electrons neglecting filled f and d shells')]


def GetAllDescriptorNames(db, tbl1, tbl2, user='sysdba', password='masterkey'):
  """ gets possible descriptor names from a database

    **Arguments**

      - db: the name of the database to use

      - tbl1: the name of the table to be used for reading descriptor values

      - tbl2: the name of the table to be used for reading notes about the
        descriptors (*descriptions of the descriptors if you like*)

      - user: the user name for DB access

      - password: the password for DB access

    **Returns**

      a 2-tuple containing:

        1) a list of column names

        2) a list of column descriptors

    **Notes**

      - this uses _Dbase.DbInfo_  and Dfunctionality for querying the database

      - it is assumed that tbl2 includes 'property' and 'notes' columns

  """
  from rdkit.Dbase.DbConnection import DbConnect
  conn = DbConnect(db, user=user, password=password)

  colNames = conn.GetColumnNames(table=tbl1)
  colDesc = map(lambda x: (x[0].upper(), x[1]), conn.GetColumns('property,notes', table=tbl2))
  for name, desc in countOptions:
    colNames.append(name)
    colDesc.append((name, desc))
  return colNames, colDesc


class CompoundDescriptorCalculator(Descriptors.DescriptorCalculator):
  """ used for calculating descriptors

   This is the central point for descriptor calculation

   **Notes**

   - There are two kinds of descriptors this cares about:

      1) *Simple Descriptors* can be calculated solely using atomic descriptor
         values and the composition of the compound.  The full list of possible
         simple descriptors is determined by the types of *Calculator Methods*
         (see below) and the contents of an atomic database.

         Simple Descriptors can be marked as *nonZeroDescriptors*.  These are used
         to winnow out atom types where particular atomic descriptors are zero
         (usually indicating that the value is unknown)

         Simple Descriptors are maintained locally in the _simpleList_

      2) *Compound Descriptors* may rely upon more complicated computation schemes
         and descriptors for the compound as a whole (e.g. structural variables, etc.).
         The full list of compound descriptors is limitless.  They are calculated using
         the _ML.Descriptors.Parser_ module.

         Compound Descriptors are maintained locally in the _compoundList_

   - This class has a some special methods which are labelled as *Calculator Method*
     These are used internally to take atomic descriptors and reduce them to a single
     simple descriptor value for a composition.  They are primarily intended for internal use.

   - a *composition vector* is a list of 2-tuples: '[(atom1name,atom1Num),...]'
     where atom1Num is the contribution of the atom to the stoichiometry of the
     compound. No assumption is made about the stoichiometries (i.e. they don't
     have to be either integral or all sum to one).

  """

  # ------------
  #  methods used to calculate descriptors
  # ------------

  def SUM(self, desc, compos):
    """ *Calculator Method*

      sums the descriptor values across the composition

      **Arguments**

        - desc: the name of the descriptor

        - compos: the composition vector

      **Returns**

        a float

    """
    res = 0.0
    for atom, num in compos:
      res = res + self.atomDict[atom][desc] * num
    return res

  def MEAN(self, desc, compos):
    """ *Calculator Method*

      averages the descriptor values across the composition

      **Arguments**

        - desc: the name of the descriptor

        - compos: the composition vector

      **Returns**

        a float

    """
    res = 0.0
    nSoFar = 0.0
    for atom, num in compos:
      res = res + self.atomDict[atom][desc] * num
      nSoFar = nSoFar + num
    return res / nSoFar

  def DEV(self, desc, compos):
    """ *Calculator Method*

      average deviation of the descriptor values across the composition

      **Arguments**

        - desc: the name of the descriptor

        - compos: the composition vector

      **Returns**

        a float

    """
    mean = self.MEAN(desc, compos)
    res = 0.0
    nSoFar = 0.0
    for atom, num in compos:
      res = res + abs(self.atomDict[atom][desc] - mean) * num
      nSoFar = nSoFar + num
    return res / nSoFar

  def MIN(self, desc, compos):
    """ *Calculator Method*

      minimum of the descriptor values across the composition

      **Arguments**

        - desc: the name of the descriptor

        - compos: the composition vector

      **Returns**

        a float

    """
    return min(map(lambda x, y=desc, z=self: z.atomDict[x[0]][y], compos))

  def MAX(self, desc, compos):
    """ *Calculator Method*

      maximum of the descriptor values across the composition

      **Arguments**

        - desc: the name of the descriptor

        - compos: the composition vector

      **Returns**

        a float

    """
    return max(map(lambda x, y=desc, z=self: z.atomDict[x[0]][y], compos))

  # ------------
  #  Other methods
  # ------------

  def ProcessSimpleList(self):
    """ Handles the list of simple descriptors

      This constructs the list of _nonZeroDescriptors_ and _requiredDescriptors_.

      There's some other magic going on that I can't decipher at the moment.

    """
    global countOptions

    self.nonZeroDescriptors = []
    lCopy = self.simpleList[:]
    tList = map(lambda x: x[0], countOptions)
    for entry in lCopy:
      if 'NONZERO' in entry[1]:
        if entry[0] not in tList:
          self.nonZeroDescriptors.append('%s != 0' % entry[0])
        if len(entry[1]) == 1:
          self.simpleList.remove(entry)
        else:
          self.simpleList[self.simpleList.index(entry)][1].remove('NONZERO')
    self.requiredDescriptors = map(lambda x: x[0], self.simpleList)
    for entry in tList:
      if entry in self.requiredDescriptors:
        self.requiredDescriptors.remove(entry)

  def ProcessCompoundList(self):
    """ Adds entries from the _compoundList_ to the list of _requiredDescriptors_

      Each compound descriptor is surveyed.  Any atomic descriptors it requires
      are added to the list of _requiredDescriptors_ to be pulled from the database.

    """
    # add in the atomic descriptors we will need
    for entry in self.compoundList:
      for atomicDesc in entry[1]:
        if atomicDesc != '' and atomicDesc not in self.requiredDescriptors:
          self.requiredDescriptors.append(atomicDesc)

  def BuildAtomDict(self):
    """ builds the local atomic dict

     We don't want to keep around all descriptor values for all atoms, so this
     method takes care of only pulling out the descriptors in which we are
     interested.

     **Notes**

       - this uses _chemutils.GetAtomicData_ to actually pull the data

    """
    self.ProcessSimpleList()
    self.ProcessCompoundList()

    self.atomDict = {}
    whereString = ' and '.join(self.nonZeroDescriptors)
    if whereString != '':
      whereString = 'where ' + whereString
    chemutils.GetAtomicData(self.atomDict, self.requiredDescriptors, self.dbName, self.dbTable,
                            whereString, self.dbUser, self.dbPassword, includeElCounts=1)

  def CalcSimpleDescriptorsForComposition(self, compos='', composList=None):
    """ calculates all simple descriptors for a given composition

      **Arguments**

        - compos: a string representation of the composition

        - composList: a *composVect*

        The client must provide either _compos_ or _composList_.  If both are
        provided, _composList_ takes priority.

      **Returns**
        the list of descriptor values

      **Notes**

        - when _compos_ is provided, this uses _chemutils.SplitComposition_
          to split the composition into its individual pieces

        - if problems are encountered because of either an unknown descriptor or
          atom type, a _KeyError_ will be raised.

    """
    if composList is None:
      composList = chemutils.SplitComposition(compos)
    try:
      res = []
      for descName, targets in self.simpleList:
        for target in targets:
          try:
            method = getattr(self, target)
          except AttributeError:
            print('Method %s does not exist' % (target))
          else:
            res.append(method(descName, composList))
    except KeyError as msg:
      print('composition %s caused problems' % composList)
      raise KeyError(msg)
    return res

  def CalcCompoundDescriptorsForComposition(self, compos='', composList=None, propDict={}):
    """ calculates all simple descriptors for a given composition

      **Arguments**

        - compos: a string representation of the composition

        - composList: a *composVect*

        - propDict: a dictionary containing the properties of the composition
          as a whole (e.g. structural variables, etc.)

        The client must provide either _compos_ or _composList_.  If both are
        provided, _composList_ takes priority.

      **Returns**
        the list of descriptor values

      **Notes**

        - when _compos_ is provided, this uses _chemutils.SplitComposition_
          to split the composition into its individual pieces

    """
    if composList is None:
      composList = chemutils.SplitComposition(compos)
    res = []
    for cl in self.compoundList:
      val = Parser.CalcSingleCompoundDescriptor(composList, cl[1:], self.atomDict, propDict)
      res.append(val)
    return res

  def CalcDescriptorsForComposition(self, composVect, propDict):
    """ calculates all descriptors for a given composition

      **Arguments**

        - compos: a string representation of the composition

        - propDict: a dictionary containing the properties of the composition
          as a whole (e.g. structural variables, etc.). These are used to
          generate Compound Descriptors

      **Returns**
        the list of all descriptor values

      **Notes**

        - this uses _chemutils.SplitComposition_
          to split the composition into its individual pieces

    """
    composList = chemutils.SplitComposition(composVect[0])
    try:
      r1 = self.CalcSimpleDescriptorsForComposition(composList=composList)
    except KeyError:
      res = []
    else:
      r2 = self.CalcCompoundDescriptorsForComposition(composList=composList, propDict=propDict)
      res = r1 + r2

    return tuple(res)

  CalcDescriptors = CalcDescriptorsForComposition

  def GetDescriptorNames(self):
    """ returns a list of the names of the descriptors this calculator generates

    """
    if self.descriptorNames is not None:
      return self.descriptorNames
    else:
      res = []
      for descName, targets in self.simpleList:
        for target in targets:
          if hasattr(self, target):
            res.append('%s_%s' % (target, descName))
          else:
            print('Method %s does not exist' % (target))
      for entry in self.compoundList:
        res.append(entry[0])
      self.descriptorNames = res[:]
      return tuple(res)

  def __init__(self, simpleList, compoundList=None, dbName=None, dbTable='atomic_data',
               dbUser='sysdba', dbPassword='masterkey'):
    """ Constructor

      **Arguments**

        - simpleList: list of simple descriptors to be calculated
              (see below for format)

        - compoundList: list of compound descriptors to be calculated
              (see below for format)

        - dbName: name of the atomic database to be used

        - dbTable: name the table in _dbName_ which has atomic data

        - dbUser: user name for DB access

        - dbPassword: password for DB access

      **Note**

        - format of simpleList:
           a list of 2-tuples containing:

              1) name of the atomic descriptor

              2) a list of operations on that descriptor (e.g. NonZero, Max, etc.)
                 These must correspond to the *Calculator Method* names above.

        - format of compoundList:
           a list of 2-tuples containing:

              1) name of the descriptor to be calculated

              2) list of selected atomic descriptor names (define $1, $2, etc.)

              3) list of selected compound descriptor names (define $a, $b, etc.)

              4) text formula defining the calculation (see _Parser_)

    """

    if dbName is None:
      dbName = RDConfig.RDDataDatabase

    Descriptors.DescriptorCalculator.__init__(self)
    self.simpleList = [(x[0].upper(), [y.upper() for y in x[1]]) for x in simpleList]
    self.descriptorNames = None
    self.compoundList = compoundList
    if self.compoundList is None:
      self.compoundList = []
    self.dbName = dbName
    self.dbTable = dbTable
    self.dbUser = dbUser
    self.dbPassword = dbPassword


def _exampleCode():
  d = [('DED', ['NonZero', 'Mean', 'Dev']), ('M_B_electroneg', ['NonZero']),
       ('Cov_rad', ['Max', 'Min'])]
  o = CompoundDescriptorCalculator(d)
  o.BuildAtomDict()
  print('len:', len(o.atomDict.keys()))
  for key in list(o.atomDict)[-4:-1]:
    print(key, o.atomDict[key])

  print('descriptors:', o.GetDescriptorNames())
  composList = ['Nb', 'Nb3', 'NbPt', 'Nb2Pt']
  for compos in composList:
    descs = o.CalcSimpleDescriptorsForComposition(compos)
    print(compos, descs)


if __name__ == '__main__':  # pragma: nocover
  _exampleCode()
