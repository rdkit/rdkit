# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" contains a class to store parameters for and results from
Composite building

"""
from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.Dbase import DbModule


def SetDefaults(runDetails):
  """  initializes a details object with default values

      **Arguments**

        - details:  (optional) a _CompositeRun.CompositeRun_ object.
          If this is not provided, the global _runDetails will be used.

      **Returns**

        the initialized _CompositeRun_ object.


  """
  runDetails.nRuns = 1
  runDetails.nModels = 10
  runDetails.outName = ''
  runDetails.badName = ''
  runDetails.splitRun = 0
  runDetails.splitFrac = 0.7
  runDetails.lockRandom = 0
  runDetails.randomActivities = 0
  runDetails.shuffleActivities = 0
  runDetails.replacementSelection = 0

  #
  # Tree Parameters
  #
  runDetails.useTrees = 1
  runDetails.pruneIt = 0
  runDetails.lessGreedy = 0
  runDetails.limitDepth = -1
  runDetails.recycleVars = 0
  runDetails.randomDescriptors = 0  # toggles growing of random forests

  #
  # KNN Parameters
  #
  runDetails.useKNN = 0
  runDetails.knnDistFunc = ''
  runDetails.knnNeighs = 0

  #
  # SigTree Parameters
  #
  runDetails.useSigTrees = 0
  runDetails.useCMIM = 0
  runDetails.allowCollections = False

  #
  # Naive Bayes Classifier Parameters
  #
  runDetails.useNaiveBayes = 0
  runDetails.mEstimateVal = -1.0
  runDetails.useSigBayes = 0

  #   #
  #   # SVM Parameters
  #   #
  #   runDetails.useSVM = 0
  #   runDetails.svmKernel = SVM.radialKernel
  #   runDetails.svmType = SVM.cSVCType
  #   runDetails.svmGamma = None
  #   runDetails.svmCost = None
  #   runDetails.svmWeights = None
  #   runDetails.svmDataType = 'float'
  #   runDetails.svmDegree = 3
  #   runDetails.svmCoeff = 0.0
  #   runDetails.svmEps = 0.001
  #   runDetails.svmNu = 0.5
  #   runDetails.svmCache = 40
  #   runDetails.svmShrink = 1
  #   runDetails.svmDataType='float'

  runDetails.bayesModel = 0
  runDetails.dbName = ''
  runDetails.dbUser = RDConfig.defaultDBUser
  runDetails.dbPassword = RDConfig.defaultDBPassword
  runDetails.dbWhat = '*'
  runDetails.dbWhere = ''
  runDetails.dbJoin = ''
  runDetails.qTableName = ''
  runDetails.qBounds = []
  runDetails.qBoundCount = ''
  runDetails.activityBounds = []
  runDetails.activityBoundsVals = ''
  runDetails.detailedRes = 0
  runDetails.noScreen = 0
  runDetails.threshold = 0.0
  runDetails.filterFrac = 0.0
  runDetails.filterVal = 0.0
  runDetails.modelFilterVal = 0.0
  runDetails.modelFilterFrac = 0.0
  runDetails.internalHoldoutFrac = 0.3
  runDetails.pickleDataFileName = ''
  runDetails.startAt = None
  runDetails.persistTblName = ''
  runDetails.randomSeed = (23, 42)
  runDetails.note = ''

  return runDetails


class CompositeRun:
  """ class to store parameters for and results from Composite building

   This class has a default set of fields which are added to the database.

   By default these fields are stored in a tuple, so they are immutable.  This
     is probably what you want.


  """
  fields = (("rundate", "varchar(32)"),
            ("dbName", "varchar(200)"),
            ("dbWhat", "varchar(200)"),
            ("dbWhere", "varchar(200)"),
            ("dbJoin", "varchar(200)"),
            ("tableName", "varchar(80)"),
            ("note", "varchar(120)"),
            ("shuffled", "smallint"),
            ("randomized", "smallint"),
            ("overall_error", "float"),
            ("holdout_error", "float"),
            ("overall_fraction_dropped", "float"),
            ("holdout_fraction_dropped", "float"),
            ("overall_correct_conf", "float"),
            ("overall_incorrect_conf", "float"),
            ("holdout_correct_conf", "float"),
            ("holdout_incorrect_conf", "float"),
            ("overall_result_matrix", "varchar(256)"),
            ("holdout_result_matrix", "varchar(256)"),
            ("threshold", "float"),
            ("splitFrac", "float"),
            ("filterFrac", "float"),
            ("filterVal", "float"),
            ("modelFilterVal", "float"),
            ("modelFilterFrac", "float"),
            ("nModels", "int"),
            ("limitDepth", "int"),
            ("bayesModels", "int"),
            ("qBoundCount", "varchar(3000)"),
            ("activityBoundsVals", "varchar(200)"),
            ("cmd", "varchar(500)"),
            ("model", DbModule.binaryTypeName), )

  def _CreateTable(self, cn, tblName):
    """ *Internal Use only*

    """
    names = map(lambda x: x.strip().upper(), cn.GetTableNames())
    if tblName.upper() not in names:
      curs = cn.GetCursor()
      fmt = []
      for name, value in self.fields:
        fmt.append('%s %s' % (name, value))
      fmtStr = ','.join(fmt)
      curs.execute('create table %s (%s)' % (tblName, fmtStr))
      cn.Commit()
    else:
      heads = [x.upper() for x in cn.GetColumnNames()]
      curs = cn.GetCursor()
      for name, value in self.fields:
        if name.upper() not in heads:
          curs.execute('alter table %s add %s %s' % (tblName, name, value))
      cn.Commit()

  def Store(self, db='models.gdb', table='results', user='sysdba', password='masterkey'):
    """ adds the result to a database

      **Arguments**

        - db: name of the database to use

        - table: name of the table to use

        - user&password: connection information

    """
    cn = DbConnect(db, table, user, password)
    curs = cn.GetCursor()
    self._CreateTable(cn, table)

    cols = []
    vals = []
    for name, _ in self.fields:
      try:
        v = getattr(self, name)
      except AttributeError:
        pass
      else:
        cols.append('%s' % name)
        vals.append(v)

    nToDo = len(vals)
    qs = ','.join([DbModule.placeHolder] * nToDo)
    vals = tuple(vals)

    cmd = 'insert into %s (%s) values (%s)' % (table, ','.join(cols), qs)
    curs.execute(cmd, vals)
    cn.Commit()

  def GetDataSet(self, **kwargs):
    """ Returns a MLDataSet pulled from a database using our stored
    values.

    """
    from rdkit.ML.Data import DataUtils
    data = DataUtils.DBToData(self.dbName, self.tableName, user=self.dbUser,
                              password=self.dbPassword, what=self.dbWhat, where=self.dbWhere,
                              join=self.dbJoin, **kwargs)

    return data

  def GetDataSetInfo(self, **kwargs):
    """ Returns a MLDataSet pulled from a database using our stored
    values.

    """
    conn = DbConnect(self.dbName, self.tableName)
    res = conn.GetColumnNamesAndTypes(join=self.dbJoin, what=self.dbWhat, where=self.dbWhere)
    return res
