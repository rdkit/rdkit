# $Id$
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" 

Dumping CDXML from python

"""
from elementtree.SimpleXMLWriter import XMLWriter
from pyRDKit import Chem

header = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.camsoft.com/xml/cdxml.dtd" >
"""
_fontsStartId=0

def DumpAtom(atom,w,highlightAtoms=[],
             highlightParms={'color':'4','size':'18'}):
  global _fontsStartId
  attrs = {'id':str(atom.GetIdx())}
  if atom.GetAtomicNum()!=6:
    attrs['Element'] = str(atom.GetAtomicNum())
  if atom.GetFormalCharge()!=0:
    attrs['Charge'] = str(atom.GetFormalCharge())
  #if atom.ExplicitHydrogenCount()!=0:
  #  attrs['NumHydrogens'] = str(atom.ExplicitHydrogenCount())
  #  attrs['ImplicitHydrogens'] = "0"
  w.start('n',attrs)
  if highlightAtoms and atom.GetIdx() in highlightAtoms:
    attrs = {}
    attrs.update(highlightParms)
    attrs['font']=attrs.get('font',str(_fontsStartId))
    w.start('t')
    w.start('s',attrs)
    w.data(atom.GetSymbol())
    w.end()
    w.end()
  w.end()

def DumpBond(bond,w,baseIdx,highlightBonds=[],
             highlightParms={'color':'4','LineWidth':'2'}):
  attrs = {'id':str(bond.GetIdx()+baseIdx)}
  if bond.GetBondType()!=Chem.BondType.AROMATIC:
    attrs['Order'] = str(int(bond.GetBondType()));
  else:
    attrs['Order'] = str(1.5)
    attrs['Display'] = "Dash"
  attrs['B'] = str(bond.GetBeginAtomIdx())
  attrs['E'] = str(bond.GetEndAtomIdx())
  bondDir = bond.GetBondDir()
  if bondDir==Chem.BondDir.BEGINWEDGE:
    attrs['Display'] = "WedgeBegin"
  elif bondDir==Chem.BondDir.ENDWEDGE:
    attrs['Display'] = "WedgeEnd"
  if bond.GetIdx() in highlightBonds:
    attrs.update(highlightParms)
  w.start('b',attrs)
  w.end()
  
def DumpFrag(frag,w,extrasId,highlightBonds=[],highlightAtoms=[]):
  w.start('fragment',id=str(extrasId))
  extrasId+=1
  for at in frag.GetAtoms():
    DumpAtom(at,w,highlightAtoms=highlightAtoms)
  for bond in frag.GetBonds():
    DumpBond(bond,w,frag.GetNumAtoms()+1,highlightBonds=highlightBonds)
  w.end()


def MolToCDXML(mol,outF,highlightBonds=[],highlightAtoms=[]):
  global header,_fontsStartId
  outF.write(header)
  w = XMLWriter(outF)
  extrasId = mol.GetNumAtoms()+mol.GetNumBonds()+2
  cdxml = w.start('CDXML')
  w.start('colortable')
  w.element('color',r="1",g="1",b="1")
  w.element('color',r="0",g="0",b="0")
  w.element('color',r="1",g="0",b="0")
  w.element('color',r="0",g="0",b="1")
  w.element('color',r="1",g="1",b="0")
  w.element('color',r="1",g="0",b="1")
  w.end()
  w.start('fonttable')
  _fontsStartId=extrasId
  w.element('font',id=str(extrasId),charset="iso-8859-1",name="Arial")
  extrasId+=1
  w.element('font',id=str(extrasId),charset="iso-8859-1",name="Times New Roman")
  extrasId+=1
  w.end()
  w.start('page',id=str(extrasId))
  extrasId+=1
  DumpFrag(mol,w,extrasId,highlightBonds=highlightBonds,
           highlightAtoms=highlightAtoms)
  w.end()
  w.close(cdxml)

if __name__ == '__main__':
  import sys,cStringIO
  try:
    from pyRDKit.utils import chemdraw
  except:
    chemdraw = None
  from pyRDKit import Chem
  for smi in ['c1cnccc1','C1=CN=CC=C1']:
    print '---------------------------------'
    m = Chem.MolFromSmiles(smi)
    io = cStringIO.StringIO()
    MolToCDXML(m,io,highlightBonds=[1,2,3],highlightAtoms=[3])
    txt = io.getvalue()
    txt = txt.replace('><','>\n<')
    cdxml = txt
    if chemdraw:
      chemdraw.CDXDisplay(cdxml)
    print cdxml
