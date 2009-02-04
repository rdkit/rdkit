# $Id$
#
# Copyright (C) 2001-2008  greg Landrum
#
#   @@ All Rights Reserved  @@
#
""" contains code for standardization of data matrices for clustering


"""
from ML.Data import Stats

def StdDev(mat):
  """ the standard deviation classifier

   This uses _ML.Data.Stats.StandardizeMatrix()_ to do the work
   
  """
  return Stats.StandardizeMatrix(mat)
    
methods = [
  ("None",lambda x:x,"No Standardization"),
  ("Standard Deviation",StdDev,"Use the standard deviation"),
  ]

