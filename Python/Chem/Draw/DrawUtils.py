# $Id$
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
import os
from MolDrawing import MolDrawing
from pyRDKit.sping.SVG.pidSVG import SVGCanvas as Canvas
from pyRDKit.sping import pid
import tempfile

_batikLoc=os.environ.get('BATIK_HOME','/usr/local/src/batik-1.5.1/')
_batikJar='batik-rasterizer.jar'
_javaExe='/usr/java/j2sdk1.4.2_05/bin/java'

def Init(batikLoc=None,batikJar=None,javaExe=None):
  global _batikLoc,_batikJar,_javaExe
  if batikLoc: _batikLoc = batikLoc
  if batikJar: _batikJar = batikJar
  if javaExe:  _javaExe  = javaExe

  if os.path.exists(os.path.join(_batikLoc,_batikJar)):
    return 1
  else:
    return 0
    
def MolsToJpeg(mols,filenames,size=(200,200),dblSize=0,frame=0,
              highlightAtoms=None,
              verbose=0):
  if isinstance(filenames,basestring):
    mols = [mols]
    filenames = [filenames]
    if highlightAtoms: highlightAtoms = [highlightAtoms]

  cmd = _javaExe
  batik = os.path.join(_batikLoc,_batikJar)
  javaArgs='-Djava.awt.headless=true -jar %s'%batik
  canvas = Canvas(size=size)
  drawing = MolDrawing(canvas)

  if dblSize:
    drawing.atomLabelFontSize=18
    drawing.bondLineWidth *= 2
    drawing.dblBondOffset *= 2
    
  if verbose: sys.stderr.write('generating SVG\n')
  svgNames = []
  for i,mol in enumerate(mols):
    filename = filenames[i]
    if highlightAtoms:
      highlight=highlightAtoms[i]
    else:
      highlight=None
    canvas.clear()
    drawing.AddMol(mol,canvas=canvas,highlightAtoms=highlight)
    svgName = '%s.svg'%(filename.split('.jpg')[0])
    if frame:
      canvas.drawRect(0,0,size[0]-1,size[1]-1,edgeColor=pid.black)
    canvas.save(svgName)
    svgNames.append(svgName)

  if verbose: sys.stderr.write('converting to jpg:\n')
  fNames=' %s'%(' '.join(svgNames))
  fmtArgs='-m image/jpeg -dpi 100 -w %d -h %d -q .80'%(size[0],size[1])

  try:
    #res = os.spawnlp(os.P_WAIT,cmd,javaArgs,'-m image/jpeg','-dpi 100',
    #          '-q .80','-w',size[0],'-h',size[1],fNames)
    arg = '%s %s %s %s &> /dev/null'%(cmd,javaArgs,fmtArgs,fNames)
    os.system(arg)
  except:
    import traceback
    if verbose: traceback.print_exc()
    ok = 0
  else:
    ok = 1

    for name in svgNames:
      try:
        #os.unlink(name)
        pass
      except:
        import traceback
        if verbose: traceback.print_exc()
  return ok
def SmilesToJpeg(smiles,filenames,**kwargs):
  from pyRDKit import Chem
  if isinstance(filenames,basestring):
    smiles = [smiles]
    filenames = [filenames]
  mols = []
  for smi in smiles:
    mol = Chem.MolFromSmiles(smi)
    Chem.Compute2DCoords(mol)
    mols.append(mol)
  return MolsToJpeg(mols,filenames,**kwargs)
    
if __name__=='__main__':
  from pyRDKit import Chem
  import tempfile
  fN1 = tempfile.mktemp('.jpg')
  print fN1
  SmilesToJpeg(['c1ccccc1C(=O)O','c1ccccc1C(=O)N'],
               [fN1,'foo2.jpg'],size=(300,300),verbose=1,
               highlightAtoms=[range(3),range(4)],dblSize=1,
               frame=1)
  
