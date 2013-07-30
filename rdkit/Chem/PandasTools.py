'''
Importing pandasTools enables several features that allow for using RDKit molecules as columns of a Pandas dataframe.
If the dataframe is containing a molecule format in a column (e.g. smiles), like in this example:
>>> from rdkit.Chem import PandasTools
>>> import pandas as pd
>>> import os
>>> from rdkit import RDConfig
>>> antibiotics = pd.DataFrame(columns=['Name','Smiles'])
>>> antibiotics = antibiotics.append({'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C','Name':'Penicilline G'}, ignore_index=True)#Penicilline G
>>> antibiotics = antibiotics.append({'Smiles':'CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O','Name':'Tetracycline'}, ignore_index=True)#Tetracycline
>>> antibiotics = antibiotics.append({'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C','Name':'Ampicilline'}, ignore_index=True)#Ampicilline
>>> print antibiotics.columns
Index([Name, Smiles], dtype=object)
>>> print antibiotics
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
1   Tetracycline  CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4...
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...

a new column can be created holding the respective RDKit molecule objects. The fingerprint can be included to accelerate substructure searches on the dataframe.

>>> PandasTools.AddMoleculeColumnToFrame(antibiotics,'Smiles','Molecule',includeFingerprints=True)
>>> print antibiotics.columns
Index([Name, Smiles, Molecule], dtype=object)

A substructure filter can be applied on the dataframe using the RDKit molecule column, because the ">=" operator has been modified to work as a substructure check.
Such the antibiotics containing the beta-lactam ring "C1C(=O)NC1" can be obtained by

>>> beta_lactam = Chem.MolFromSmiles('C1C(=O)NC1')
>>> beta_lactam_antibiotics = antibiotics[antibiotics['Molecule'] >= beta_lactam]
>>> print beta_lactam_antibiotics[['Name','Smiles']]
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...


It is also possible to load an SDF file can be load into a dataframe.

>>> sdfFile = os.path.join(RDConfig.RDDataDir,'NCI/first_200.props.sdf')
>>> frame = PandasTools.LoadSDF(sdfFile,smilesName='SMILES',molColName='Molecule',includeFingerprints=True)
>>> frame.info # doctest: +SKIP
<bound method DataFrame.info of <class 'pandas.core.frame.DataFrame'>
Int64Index: 200 entries, 0 to 199
Data columns:
AMW                       200  non-null values
CLOGP                     200  non-null values
CP                        200  non-null values
CR                        200  non-null values
DAYLIGHT.FPG              200  non-null values
DAYLIGHT_CLOGP            200  non-null values
FP                        200  non-null values
ID                        200  non-null values
ISM                       200  non-null values
LIPINSKI_VIOLATIONS       200  non-null values
NUM_HACCEPTORS            200  non-null values
NUM_HDONORS               200  non-null values
NUM_HETEROATOMS           200  non-null values
NUM_LIPINSKIHACCEPTORS    200  non-null values
NUM_LIPINSKIHDONORS       200  non-null values
NUM_RINGS                 200  non-null values
NUM_ROTATABLEBONDS        200  non-null values
P1                        30  non-null values
SMILES                    200  non-null values
Molecule                  200  non-null values
dtypes: object(20)>

In order to support rendering the molecules as images in the HTML export of the dataframe, the __str__ method is monkey-patched to return a base64 encoded PNG:
>>> molX = Chem.MolFromSmiles('Fc1cNc2ccccc12')
>>> print molX # doctest: +SKIP
<img src="data:image/png;base64,..." alt="Mol"/>
This can be reverted using the ChangeMoleculeRendering method
>>> ChangeMoleculeRendering(renderer='String')
>>> print molX # doctest: +SKIP
<rdkit.Chem.rdchem.Mol object at 0x10d179440>
>>> ChangeMoleculeRendering(renderer='PNG')
>>> print molX # doctest: +SKIP
<img src="data:image/png;base64,..." alt="Mol"/>
'''


from rdkit import Chem
from rdkit.Chem import Draw
from base64 import b64encode
from StringIO import StringIO
import types
try:
  import pandas as pd
  if 'display.width' in  pd.core.config._registered_options:
    pd.set_option('display.width',100000000000)
  if 'display.height' in  pd.core.config._registered_options:
    pd.set_option('display.height',100000000000)
  if 'display.max_colwidth' in  pd.core.config._registered_options:
    pd.set_option('display.max_colwidth',100000000000)
  #saves the default pandas rendering to allow restauration
  defPandasRendering = pd.core.frame.DataFrame.to_html
except ImportError:
  pd = None



def patchPandasHTMLrepr(self):
  '''
  Patched default escaping of HTML control characters to allow molecule image rendering dataframes
  '''
  formatter = pd.core.format.DataFrameFormatter(self,buf=None,columns=None,col_space=None,colSpace=None,header=True,index=True,
                                               na_rep='NaN',formatters=None,float_format=None,sparsify=None,index_names=True,
                                               justify = None, force_unicode=None,bold_rows=True,classes=None,escape=False)
  formatter.to_html()
  html = formatter.buf.getvalue()
  return html

def patchPandasHeadMethod(self,n=5):
  '''Ensure inheritance of patched to_html in "head" subframe
  '''
  df = self[:n]
  df.to_html = types.MethodType(patchPandasHTMLrepr,df)
  df.head = types.MethodType(patchPandasHeadMethod,df)
  return df

def _get_image(x):
  """displayhook function for PIL Images, rendered as PNG"""
  import pandas as pd
  sio = StringIO()    
  x.save(sio,format='PNG')
  s = b64encode(sio.getvalue())
  pd.set_option('display.max_columns',len(s)+1000)
  pd.set_option('display.max_rows',len(s)+1000)
  if len(s)+100 > pd.get_option("display.max_colwidth"):
    pd.set_option("display.max_colwidth",len(s)+1000)
  return s

from rdkit import DataStructs

try:
  from rdkit.Avalon import pyAvalonTools as pyAvalonTools
  _fingerprinter=lambda x,y:pyAvalonTools.GetAvalonFP(x,isQuery=y,bitFlags=pyAvalonTools.avalonSSSBits)
except ImportError:
  _fingerprinter=lambda x,y:Chem.PatternFingerprint(x,fpSize=2048)

def _molge(x,y):
  """Allows for substructure check using the >= operator (X has substructure Y -> X >= Y) by
  monkey-patching the __ge__ function 
  This has the effect that the pandas/numpy rowfilter can be used for substructure filtering (filtered = dframe[dframe['RDKitColumn'] >= SubstructureMolecule])
  """
  if x is None or y is None: return False
  if hasattr(x,'_substructfp'):
    if not hasattr(y,'_substructfp'):
      y._substructfp=_fingerprinter(y,True)
    if not DataStructs.AllProbeBitsMatch(y._substructfp,x._substructfp):
      return False
  return x.HasSubstructMatch(y)
            
Chem.Mol.__ge__ = _molge # lambda x,y: x.HasSubstructMatch(y)

def PrintAsBase64PNGString(x,renderer = None):
  '''returns the molecules as base64 encoded PNG image
  '''
  return '<img src="data:image/png;base64,%s" alt="Mol"/>'%_get_image(Draw.MolToImage(x))

def PrintDefaultMolRep(x):
  return str(x.__repr__())
	
#Chem.Mol.__str__ = lambda x: '<img src="data:image/png;base64,%s" alt="Mol"/>'%get_image(Draw.MolToImage(x))
Chem.Mol.__str__ = PrintAsBase64PNGString

def _MolPlusFingerprintFromSmiles(smi):
  '''Precomputes fingerprints and stores results in molecule objects to accelerate substructure matching
  '''
  m = Chem.MolFromSmiles(smi)
  if m is not None:
    m._substructfp=_fingerprinter(m,False)
  return m
  
def RenderImagesInAllDataFrames(images=True):
  '''Changes the default dataframe rendering to not escape HTML characters, thus allowing rendered images in all dataframes.
  IMPORTANT: THIS IS A GLOBAL CHANGE THAT WILL AFFECT TO COMPLETE PYTHON SESSION. If you want to change the rendering only 
  for a single dataframe use the "ChangeMoleculeRendering" method instead.
  '''
  if images:
    pd.core.frame.DataFrame.to_html = patchPandasHTMLrepr
  else:
    pd.core.frame.DataFrame.to_html = defPandasRendering
	    

def AddMoleculeColumnToFrame(frame, smilesCol='Smiles', molCol = 'ROMol',includeFingerprints=False):
  '''Converts the molecules contains in "smilesCol" to RDKit molecules and appends them to the dataframe "frame" using the specified column name.
  If desired, a fingerprint can be computed and stored with the molecule objects to accelerate substructure matching
  '''
  if not includeFingerprints:
    frame[molCol]=frame.apply(lambda x: Chem.MolFromSmiles(x[smilesCol]), axis=1)
  else:
    frame[molCol]=frame.apply(lambda x: _MolPlusFingerprintFromSmiles(x[smilesCol]), axis=1) 
  RenderImagesInAllDataFrames(images=True)
  #frame.to_html = types.MethodType(patchPandasHTMLrepr,frame)
  #frame.head = types.MethodType(patchPandasHeadMethod,frame)
  
  
def ChangeMoleculeRendering(frame=None, renderer='PNG'):
  '''Allows to change the rendering of the molecules between base64 PNG images and string representations. 
  This serves two purposes: First it allows to avoid the generation of images if this is not desired and, secondly, it allows to enable image rendering for
  newly created dataframe that already contains molecules, without having to rerun the time-consuming AddMoleculeColumnToFrame. Note: this behaviour is, because some pandas methods, e.g. head()
  returns a new dataframe instance that uses the default pandas rendering (thus not drawing images for molecules) instead of the monkey-patched one.
  '''
  if renderer == 'String':
	Chem.Mol.__str__ = PrintDefaultMolRep
  else:
	Chem.Mol.__str__ = PrintAsBase64PNGString
  if frame is not None:
	frame.to_html = types.MethodType(patchPandasHTMLrepr,frame)
    
def LoadSDF(filename, smilesName='SMILES', idName='ID',molColName = 'ROMol',includeFingerprints=False):
  """ Read file in SDF format and return as Panda data frame """
  df = None
  with open(filename, 'rU') as f:
    for i, mol in enumerate(Chem.ForwardSDMolSupplier(f)):
      if mol is None: continue
      row = dict((k, mol.GetProp(k)) for k in mol.GetPropNames())
      if mol.HasProp('_Name'): row[idName] = mol.GetProp('_Name')
      row[smilesName] = Chem.MolToSmiles(mol)
      row = pd.DataFrame(row, index=[i])
      if df is None:
        df = row
      else:
        df = df.append(row)
  AddMoleculeColumnToFrame(df, smilesCol=smilesName, molCol = molColName,includeFingerprints=includeFingerprints)
  return df  



if __name__ == "__main__":
  import sys
  if pd is None:
    print >>sys.stderr,"pandas installation not found, skipping tests"
  else:
    v = pd.version.version.split('.')
    if v[0]=='0' and int(v[1])<10:
      print >>sys.stderr,"pandas installation >=0.10 not found, skipping tests"
    else:
      import doctest
      failed,tried=doctest.testmod(optionflags=doctest.ELLIPSIS+doctest.NORMALIZE_WHITESPACE)
      if failed:
        sys.exit(failed)

# $Id$
#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
