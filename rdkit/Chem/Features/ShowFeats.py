# $Id: ShowFeats.py 537 2007-08-20 14:54:35Z landrgr1 $
#
# Created by Greg Landrum Aug 2006
#
#

_version = "0.3.2"

_usage = """
   ShowFeats [optional args] <filenames>

  if "-" is provided as a filename, data will be read from stdin (the console)
"""

_welcomeMessage = "This is ShowFeats version %s" % (_version)

import math

#set up the logger:
from rdkit import RDLogger as logging

logger = logging.logger()
logger.setLevel(logging.INFO)

from rdkit import Geometry
from rdkit.Chem.Features import FeatDirUtilsRD as FeatDirUtils

_featColors = {
  'Donor': (0, 1, 1),
  'Acceptor': (1, 0, 1),
  'NegIonizable': (1, 0, 0),
  'PosIonizable': (0, 0, 1),
  'ZnBinder': (1, .5, .5),
  'Aromatic': (1, .8, .2),
  'LumpedHydrophobe': (.5, .25, 0),
  'Hydrophobe': (.5, .25, 0),
}


def _getVectNormal(v, tol=1e-4):
  if math.fabs(v.x) > tol:
    res = Geometry.Point3D(v.y, -v.x, 0)
  elif math.fabs(v.y) > tol:
    res = Geometry.Point3D(-v.y, v.x, 0)
  elif math.fabs(v.z) > tol:
    res = Geometry.Point3D(1, 0, 0)
  else:
    raise ValueError('cannot find normal to the null vector')
  res.Normalize()
  return res


_canonArrowhead = None


def _buildCanonArrowhead(headFrac, nSteps, aspect):
  global _canonArrowhead
  startP = RDGeometry.Point3D(0, 0, headFrac)
  _canonArrowhead = [startP]

  scale = headFrac * aspect
  baseV = RDGeometry.Point3D(scale, 0, 0)
  _canonArrowhead.append(baseV)

  twopi = 2 * math.pi
  for i in range(1, nSteps):
    v = RDGeometry.Point3D(scale * math.cos(i * twopi), scale * math.sin(i * twopi), 0)
    _canonArrowhead.append(v)


_globalArrowCGO = []
_globalSphereCGO = []
# taken from pymol's cgo.py
BEGIN = 2
END = 3
TRIANGLE_FAN = 6
COLOR = 6
VERTEX = 4
NORMAL = 5
SPHERE = 7
CYLINDER = 9
ALPHA = 25


def _cgoArrowhead(viewer, tail, head, radius, color, label, headFrac=0.3, nSteps=10, aspect=.5):
  global _globalArrowCGO
  delta = head - tail
  normal = _getVectNormal(delta)
  delta.Normalize()

  dv = head - tail
  dv.Normalize()
  dv *= headFrac
  startP = head

  normal *= headFrac * aspect

  cgo = [
    BEGIN, TRIANGLE_FAN, COLOR, color[0], color[1], color[2], NORMAL, dv.x, dv.y, dv.z, VERTEX,
    head.x + dv.x, head.y + dv.y, head.z + dv.z
  ]
  base = [
    BEGIN, TRIANGLE_FAN, COLOR, color[0], color[1], color[2], NORMAL, -dv.x, -dv.y, -dv.z, VERTEX,
    head.x, head.y, head.z
  ]
  v = startP + normal
  cgo.extend([NORMAL, normal.x, normal.y, normal.z])
  cgo.extend([VERTEX, v.x, v.y, v.z])
  base.extend([VERTEX, v.x, v.y, v.z])
  for i in range(1, nSteps):
    v = FeatDirUtils.ArbAxisRotation(360. / nSteps * i, delta, normal)
    cgo.extend([NORMAL, v.x, v.y, v.z])
    v += startP
    cgo.extend([VERTEX, v.x, v.y, v.z])
    base.extend([VERTEX, v.x, v.y, v.z])

  cgo.extend([NORMAL, normal.x, normal.y, normal.z])
  cgo.extend([VERTEX, startP.x + normal.x, startP.y + normal.y, startP.z + normal.z])
  base.extend([VERTEX, startP.x + normal.x, startP.y + normal.y, startP.z + normal.z])
  cgo.append(END)
  base.append(END)
  cgo.extend(base)

  #viewer.server.renderCGO(cgo,label)
  _globalArrowCGO.extend(cgo)


def ShowArrow(viewer, tail, head, radius, color, label, transparency=0, includeArrowhead=True):
  global _globalArrowCGO
  if transparency:
    _globalArrowCGO.extend([ALPHA, 1 - transparency])
  else:
    _globalArrowCGO.extend([ALPHA, 1])
  _globalArrowCGO.extend([
    CYLINDER,
    tail.x,
    tail.y,
    tail.z,
    head.x,
    head.y,
    head.z,
    radius * .10,
    color[0],
    color[1],
    color[2],
    color[0],
    color[1],
    color[2],
  ])

  if includeArrowhead:
    _cgoArrowhead(viewer, tail, head, radius, color, label)


def ShowMolFeats(mol, factory, viewer, radius=0.5, confId=-1, showOnly=True, name='',
                 transparency=0.0, colors=None, excludeTypes=[], useFeatDirs=True, featLabel=None,
                 dirLabel=None, includeArrowheads=True, writeFeats=False, showMol=True,
                 featMapFile=False):
  global _globalSphereCGO
  if not name:
    if mol.HasProp('_Name'):
      name = mol.GetProp('_Name')
    else:
      name = 'molecule'
  if not colors:
    colors = _featColors

  if showMol:
    viewer.ShowMol(mol, name=name, showOnly=showOnly, confId=confId)

  molFeats = factory.GetFeaturesForMol(mol)
  if not featLabel:
    featLabel = f'{name}-feats'
    viewer.server.resetCGO(featLabel)
  if not dirLabel:
    dirLabel = featLabel + "-dirs"
    viewer.server.resetCGO(dirLabel)

  for feat in molFeats:
    family = feat.GetFamily()
    if family in excludeTypes:
      continue
    pos = feat.GetPos(confId)
    color = colors.get(family, (.5, .5, .5))

    if transparency:
      _globalSphereCGO.extend([ALPHA, 1 - transparency])
    else:
      _globalSphereCGO.extend([ALPHA, 1])
    _globalSphereCGO.extend(
      [COLOR, color[0], color[1], color[2], SPHERE, pos.x, pos.y, pos.z, radius])
    if writeFeats:
      aidText = ' '.join([str(x + 1) for x in feat.GetAtomIds()])
      print(f'{family}\t{pos.x:.3f}\t{pos.y:.3f}\t{pos.z:.3f}\t1.0\t# {aidText}')

    if featMapFile:
      print(f"  family={family} pos=({pos.x:.3f}, {pos.y:.3f}, {pos.z:.3f}) weight=1.0", end='',
            file=featMapFile)

    if useFeatDirs:
      ps = []
      if family == 'Aromatic':
        ps, _ = FeatDirUtils.GetAromaticFeatVects(mol.GetConformer(confId), feat.GetAtomIds(), pos,
                                                  scale=1.0)

      elif family == 'Donor':
        aids = feat.GetAtomIds()
        if len(aids) == 1:
          FeatVectsDictMethod = {
            1: FeatDirUtils.GetDonor1FeatVects,
            2: FeatDirUtils.GetDonor2FeatVects,
            3: FeatDirUtils.GetDonor3FeatVects,
          }
          featAtom = mol.GetAtomWithIdx(aids[0])
          numHvyNbrs = len([1 for x in featAtom.GetNeighbors() if x.GetAtomicNum() > 1])
          ps, _ = FeatVectsDictMethod[numHvyNbrs](mol.GetConformer(confId), aids, scale=1.0)

      elif family == 'Acceptor':
        aids = feat.GetAtomIds()
        if len(aids) == 1:
          FeatVectsDictMethod = {
            1: FeatDirUtils.GetDonor1FeatVects,
            2: FeatDirUtils.GetDonor2FeatVects,
            3: FeatDirUtils.GetDonor3FeatVects,
          }
          featAtom = mol.GetAtomWithIdx(aids[0])
          numHvyNbrs = len([x for x in featAtom.GetNeighbors() if x.GetAtomicNum() > 1])
          ps, _ = FeatVectsDictMethod[numHvyNbrs](mol.GetConformer(confId), aids, scale=1.0)

      for tail, head in ps:
        ShowArrow(viewer, tail, head, radius, color, dirLabel, transparency=transparency,
                  includeArrowhead=includeArrowheads)
        if featMapFile:
          vect = head - tail
          print(f'dir=({vect.x:.3f}, {vect.y:.3f}, {vect.z:.3f})', end='', file=featMapFile)

    if featMapFile:
      aidText = ' '.join([str(x + 1) for x in feat.GetAtomIds()])
      print(f'# {aidText}', file=featMapFile)


import os
# --- ----  --- ----  --- ----  --- ----  --- ----  --- ----
import sys
from optparse import OptionParser

from rdkit import RDConfig

parser = OptionParser(_usage, version='%prog ' + _version)

parser.add_option('-x', '--exclude', default='',
                  help='provide a list of feature names that should be excluded')
parser.add_option('-f', '--fdef', default=os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'),
                  help='provide the name of the feature definition (fdef) file.')
parser.add_option('--noDirs', '--nodirs', dest='useDirs', default=True, action='store_false',
                  help='do not draw feature direction indicators')
parser.add_option('--noHeads', dest='includeArrowheads', default=True, action='store_false',
                  help='do not draw arrowheads on the feature direction indicators')
parser.add_option('--noClear', '--noClear', dest='clearAll', default=False, action='store_true',
                  help='do not clear PyMol on startup')
parser.add_option('--noMols', '--nomols', default=False, action='store_true',
                  help='do not draw the molecules')
parser.add_option('--writeFeats', '--write', default=False, action='store_true',
                  help='print the feature information to the console')
parser.add_option('--featMapFile', '--mapFile', default='',
                  help='save a feature map definition to the specified file')
parser.add_option('--verbose', default=False, action='store_true', help='be verbose')

if __name__ == '__main__':
  from rdkit import Chem
  from rdkit.Chem import AllChem
  from rdkit.Chem.PyMol import MolViewer

  options, args = parser.parse_args()
  if len(args) < 1:
    parser.error('please provide either at least one sd or mol file')

  try:
    v = MolViewer()
  except Exception:
    logger.error(
      'Unable to connect to PyMol server.\nPlease run ~landrgr1/extern/PyMol/launch.sh to start it.'
    )
    sys.exit(1)
  if options.clearAll:
    v.DeleteAll()

  try:
    fdef = open(options.fdef, 'r').read()
  except IOError:
    logger.error('ERROR: Could not open fdef file %s' % options.fdef)
    sys.exit(1)

  factory = AllChem.BuildFeatureFactoryFromString(fdef)

  if options.writeFeats:
    print('# Family   \tX    \tY   \tZ   \tRadius\t # Atom_ids')

  if options.featMapFile:
    if options.featMapFile == '-':
      options.featMapFile = sys.stdout
    else:
      options.featMapFile = file(options.featMapFile, 'w+')
    print('# Feature map generated by ShowFeats v%s' % _version, file=options.featMapFile)
    print("ScoreMode=All", file=options.featMapFile)
    print("DirScoreMode=Ignore", file=options.featMapFile)
    print("BeginParams", file=options.featMapFile)
    for family in factory.GetFeatureFamilies():
      print("   family=%s width=1.0 radius=3.0" % family, file=options.featMapFile)
    print("EndParams", file=options.featMapFile)
    print("BeginPoints", file=options.featMapFile)

  i = 1
  for midx, molN in enumerate(args):
    if molN != '-':
      featLabel = f'{molN}_Feats' % molN
    else:
      featLabel = f'Mol{midx + 1}_Feats'

    v.server.resetCGO(featLabel)
    # this is a big of kludgery to work around what seems to be a pymol cgo bug:
    v.server.sphere((0, 0, 0), .01, (1, 0, 1), featLabel)
    dirLabel = featLabel + "-dirs"

    v.server.resetCGO(dirLabel)
    # this is a big of kludgery to work around what seems to be a pymol cgo bug:
    v.server.cylinder((0, 0, 0), (.01, .01, .01), .01, (1, 0, 1), dirLabel)

    if molN != '-':
      try:
        ms = Chem.SDMolSupplier(molN)
      except Exception:
        logger.error('Problems reading input file: %s' % molN)
        ms = []
    else:
      ms = Chem.SDMolSupplier()
      ms.SetData(sys.stdin.read())

    for m in ms:
      nm = f'Mol_{i}'
      if m.HasProp('_Name'):
        nm += '_' + m.GetProp('_Name')
      if options.verbose:
        if m.HasProp('_Name'):
          print("#Molecule: ", m.GetProp('_Name'))
        else:
          print("#Molecule: ", nm)
      ShowMolFeats(m, factory, v, transparency=0.25, excludeTypes=options.exclude, name=nm,
                   showOnly=False, useFeatDirs=options.useDirs, featLabel=featLabel,
                   dirLabel=dirLabel, includeArrowheads=options.includeArrowheads,
                   writeFeats=options.writeFeats, showMol=not options.noMols,
                   featMapFile=options.featMapFile)
      i += 1
      if not i % 100:
        logger.info(f"Done {i} poses")

    if ms:
      v.server.renderCGO(_globalSphereCGO, featLabel, 1)
      if options.useDirs:
        v.server.renderCGO(_globalArrowCGO, dirLabel, 1)

  if options.featMapFile:
    print("EndPoints", file=options.featMapFile)

  sys.exit(0)
