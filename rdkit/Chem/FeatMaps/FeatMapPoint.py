# $Id$
#
# Copyright (C) 2006 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit.Chem import ChemicalFeatures


class FeatMapPoint(ChemicalFeatures.FreeChemicalFeature):
  weight = 0.0
  featDirs = None

  def __init__(self, *args, **kwargs):
    ChemicalFeatures.FreeChemicalFeature.__init__(self, *args, **kwargs)
    self.featDirs = []

  def initFromFeat(self, feat):
    """
    >>> from rdkit import Geometry
    >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
    >>> fmp = FeatMapPoint()
    >>> fmp.initFromFeat(sfeat)
    >>> fmp.GetFamily()==sfeat.GetFamily()
    True
    >>> fmp.GetType()==sfeat.GetType()
    True
    >>> list(fmp.GetPos())
    [0.0, 0.0, 0.0]
    >>> fmp.featDirs == []
    True

    >>> sfeat.featDirs = [Geometry.Point3D(1.0,0,0)]
    >>> fmp.initFromFeat(sfeat)
    >>> len(fmp.featDirs)
    1

    """
    self.SetFamily(feat.GetFamily())
    self.SetType(feat.GetType())
    self.SetPos(feat.GetPos())
    if hasattr(feat, 'featDirs'):
      self.featDirs = feat.featDirs[:]

  def GetDist2(self, other):
    """
    >>> from rdkit import Geometry
    >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
    >>> fmp = FeatMapPoint()
    >>> fmp.initFromFeat(sfeat)
    >>> fmp.GetDist2(sfeat)
    0.0
    >>> sfeat.SetPos(Geometry.Point3D(2,0,0))
    >>> fmp.GetDist2(sfeat)
    4.0
    """
    return (self.GetPos() - other.GetPos()).LengthSq()

  def GetDirMatch(self, other, useBest=True):
    """
    >>> from rdkit import Geometry
    >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
    >>> fmp = FeatMapPoint()
    >>> fmp.initFromFeat(sfeat)
    >>> fmp.GetDirMatch(sfeat)
    1.0

    >>> sfeat.featDirs=[Geometry.Point3D(0,0,1),Geometry.Point3D(0,0,-1)]
    >>> fmp.featDirs=[Geometry.Point3D(0,0,1),Geometry.Point3D(1,0,0)]
    >>> fmp.GetDirMatch(sfeat)
    1.0
    >>> fmp.GetDirMatch(sfeat,useBest=True)
    1.0
    >>> fmp.GetDirMatch(sfeat,useBest=False)
    0.0

    >>> sfeat.featDirs=[Geometry.Point3D(0,0,1)]
    >>> fmp.GetDirMatch(sfeat,useBest=False)
    0.5

    >>> sfeat.featDirs=[Geometry.Point3D(0,0,1)]
    >>> fmp.featDirs=[Geometry.Point3D(0,0,-1)]
    >>> fmp.GetDirMatch(sfeat)
    -1.0
    >>> fmp.GetDirMatch(sfeat,useBest=False)
    -1.0


    """
    if not self.featDirs or not other.featDirs:
      return 1.0

    if not useBest:
      accum = 0.0
    else:
      accum = -100000.0
    for sDir in self.featDirs:
      for oDir in other.featDirs:
        d = sDir.DotProduct(oDir)
        if useBest:
          if d > accum:
            accum = d
        else:
          accum += d

    if not useBest:
      accum /= len(self.featDirs) * len(other.featDirs)

    return accum


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
