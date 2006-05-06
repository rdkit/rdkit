# $Id: MatOps.py 5030 2006-03-02 18:46:51Z glandrum $
#
#  Copyright (C) 2000-2006  greg Landrum
#
#   @@ All Rights Reserved  @@
#
""" Matrix operations which may or may not come in handy some day


  **NOTE**: the two functions defined here have been moved to ML.Data.Stats

"""
from ML.Data import Stats

FormCovarianceMatrix = Stats.FormCovarianceMatrix
PrincipalComponents = Stats.PrincipalComponents

if __name__ == '__main__':
  import sys
  import files

  fileN = sys.argv[1]
  iV,dV = files.ReadDataFile(fileN)
  eVals,eVects=PrincipalComponents(iV)
  print 'eVals: ', eVals
  print 'eVects:', eVects
  

  

