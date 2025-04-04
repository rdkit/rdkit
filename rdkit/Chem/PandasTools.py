#
#  Copyright (c) 2013-2022, Novartis Institutes for BioMedical Research Inc.
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
'''
Importing pandasTools enables several features that allow for using RDKit molecules as columns of a
Pandas dataframe.
If the dataframe is containing a molecule format in a column (e.g. smiles), like in this example:

>>> import pandas as pd
>>> from rdkit.Chem import PandasTools
>>> PandasTools.InstallPandasTools() # <- only necessary during testing, you don't need to do this
>>> import os
>>> from rdkit import RDConfig
>>> antibiotics = pd.DataFrame(columns=['Name','Smiles'])
>>> antibiotics = pd.concat([antibiotics, pd.DataFrame.from_records([{'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',
...   'Name':'Penicilline G'}])], ignore_index=True) #Penicilline G
>>> antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
...   'Smiles':'CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O',
...   'Name':'Tetracycline'}])], ignore_index=True) #Tetracycline
>>> antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
...   'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C',
...   'Name':'Ampicilline'}])], ignore_index=True) #Ampicilline
>>> print([str(x) for x in  antibiotics.columns])
['Name', 'Smiles']
>>> print(antibiotics)
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
1   Tetracycline  CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4...
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...

a new column can be created holding the respective RDKit molecule objects. The fingerprint can be
included to accelerate substructure searches on the dataframe.

>>> PandasTools.AddMoleculeColumnToFrame(antibiotics,'Smiles','Molecule',includeFingerprints=True)
>>> print([str(x) for x in  antibiotics.columns])
['Name', 'Smiles', 'Molecule']

A substructure filter can be applied on the dataframe using the RDKit molecule column,
because the ">=" operator has been modified to work as a substructure check.
Such the antibiotics containing the beta-lactam ring "C1C(=O)NC1" can be obtained by

>>> beta_lactam = Chem.MolFromSmiles('C1C(=O)NC1')
>>> beta_lactam_antibiotics = antibiotics[antibiotics['Molecule'] >= beta_lactam] 
>>> print(beta_lactam_antibiotics[['Name','Smiles']])
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...


It is also possible to load an SDF file can be load into a dataframe.

>>> sdfFile = os.path.join(RDConfig.RDDataDir,'NCI/first_200.props.sdf')
>>> frame = PandasTools.LoadSDF(sdfFile,smilesName='SMILES',molColName='Molecule',
...            includeFingerprints=True)
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

The standard ForwardSDMolSupplier keywords are also available:

>>> sdfFile = os.path.join(RDConfig.RDDataDir,'NCI/first_200.props.sdf')
>>> frame = PandasTools.LoadSDF(sdfFile, smilesName='SMILES', molColName='Molecule',
...            includeFingerprints=True, removeHs=False, strictParsing=True)

Conversion to html is quite easy:

>>> PandasTools.molRepresentation = 'svg' #< the default is 'png', for png representation
>>> htm = frame.to_html()
...
>>> str(htm[:36])
'<table border="1" class="dataframe">'
>>> PandasTools.molRepresentation = 'png' #< switch back to the default


In order to support rendering the molecules as images in the HTML export of the
dataframe, we use a custom formatter for columns containing RDKit molecules,
and also disable escaping of HTML where needed.
'''

import logging
import sys
from base64 import b64encode

import numpy as np

import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, SDWriter, rdchem
from rdkit.Chem.Scaffolds import MurckoScaffold

InteractiveRenderer = None
drawOptions = None
if hasattr(rdkit, 'IPythonConsole'):
  try:
    from rdkit.Chem.Draw.IPythonConsole import InteractiveRenderer, drawOptions
  except ImportError:
    pass

from io import BytesIO
from xml.dom import minidom

log = logging.getLogger(__name__)

_originalSettings = {}
highlightSubstructures = True
molRepresentation = "png"  # supports also SVG
molSize = (200, 200)
molJustify = "center"  # supports also left, right


def _molge(x, y):
  """Allows for substructure check using the >= operator (X has substructure Y -> X >= Y) by
    monkey-patching the __ge__ function
    This has the effect that the pandas/numpy rowfilter can be used for substructure filtering
    (filtered = dframe[dframe['RDKitColumn'] >= SubstructureMolecule])
    """
  if x is None or y is None:
    return False
  if hasattr(x, '_substructfp'):
    if not hasattr(y, '_substructfp'):
      y._substructfp = _fingerprinter(y, True)
    if not DataStructs.AllProbeBitsMatch(y._substructfp, x._substructfp):
      return False
  match = x.GetSubstructMatch(y)
  x.__sssAtoms = []
  if match:
    if highlightSubstructures:
      x.__sssAtoms = list(match)
    return True
  else:
    return False


def PrintAsImageString(x):
  """Returns the molecules as base64 encoded PNG image or as SVG"""
  if highlightSubstructures and hasattr(x, '__sssAtoms'):
    highlightAtoms = x.__sssAtoms
  else:
    highlightAtoms = []
  useSVG = (molRepresentation.lower() == 'svg')
  if InteractiveRenderer and InteractiveRenderer.isEnabled(x):
    size = [max(30, s) for s in molSize]
    return InteractiveRenderer.generateHTMLBody(x, size, useSVG=useSVG)
  else:
    if useSVG:
      svg = Draw._moltoSVG(x, molSize, highlightAtoms, "", kekulize=True, drawOptions=drawOptions)
      svg = minidom.parseString(svg)
      svg = svg.getElementsByTagName('svg')[0]
      svg.attributes['viewbox'] = f'0 0 {molSize[0]} {molSize[1]}'
      svg.attributes['style'] = f'max-width: {molSize[0]}px; height: {molSize[1]}px;'
      svg.attributes['data-content'] = 'rdkit/molecule'
      return svg.toxml()
    else:
      data = Draw._moltoimg(x, molSize, highlightAtoms, "", returnPNG=True, kekulize=True,
                            drawOptions=drawOptions)
      return (
        f'<div style="width: {molSize[0]}px; height: {molSize[1]}px" data-content="rdkit/molecule">'
        '<img src="data:image/png;base64,%s" alt="Mol"/>'
        '</div>' % _get_image(data))


try:
  from rdkit.Chem import PandasPatcher
  PandasPatcher.PrintAsImageString = PrintAsImageString
  PandasPatcher.molJustify = molJustify
  PandasPatcher.InteractiveRenderer = InteractiveRenderer
except Exception:
  log.warning("Failed to patch pandas - PandasTools will have limited functionality")


def RenderImagesInAllDataFrames(images=True):
  '''Changes the default dataframe rendering to not escape HTML characters, thus allowing
    rendered images in all dataframes.
    IMPORTANT: THIS IS A GLOBAL CHANGE THAT WILL AFFECT TO COMPLETE PYTHON SESSION. If you want
    to change the rendering only for a single dataframe use the "ChangeMoleculeRendering" method
    instead.
    '''
  try:
    PandasPatcher.renderImagesInAllDataFrames(images)
  except NameError:
    log.warning("Failed to patch pandas - unable to change molecule rendering")


try:
  import pandas as pd
except ImportError:
  pass
else:

  def LoadSDF(filename, idName='ID', molColName='ROMol', includeFingerprints=False,
              isomericSmiles=True, smilesName=None, embedProps=False, removeHs=True,
              strictParsing=True, sanitize=True):
    '''Read file in SDF format and return as Pandas data frame.
      If embedProps=True all properties also get embedded in Mol objects in the molecule column.
      If molColName=None molecules would not be present in resulting DataFrame (only properties
      would be read).
      
      Sanitize boolean is passed on to Chem.ForwardSDMolSupplier sanitize. 
      If neither molColName nor smilesName are set, sanitize=false.
      '''
    if isinstance(filename, str):
      if filename.lower()[-3:] == ".gz":
        import gzip
        f = gzip.open(filename, "rb")
      else:
        f = open(filename, 'rb')
      close = f.close
    else:
      f = filename
      close = None  # don't close an open file that was passed in
    records = []
    indices = []
    if molColName is None and smilesName is None:
      sanitize = False
    for i, mol in enumerate(
        Chem.ForwardSDMolSupplier(f, sanitize=sanitize, removeHs=removeHs,
                                  strictParsing=strictParsing)):
      if mol is None:
        continue
      row = dict((k, mol.GetProp(k)) for k in mol.GetPropNames())
      if molColName is not None and not embedProps:
        for prop in mol.GetPropNames():
          mol.ClearProp(prop)
      if mol.HasProp('_Name'):
        row[idName] = mol.GetProp('_Name')
      if smilesName is not None:
        try:
          row[smilesName] = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
        except Exception:
          log.warning('No valid smiles could be generated for molecule %s', i)
          row[smilesName] = None
      if molColName is not None and not includeFingerprints:
        row[molColName] = mol
      elif molColName is not None:
        row[molColName] = _MolPlusFingerprint(mol)
      records.append(row)
      indices.append(i)

    if close is not None:
      close()
    df = pd.DataFrame(records, index=indices)
    ChangeMoleculeRendering(df)
    return df

  def RGroupDecompositionToFrame(groups, mols, include_core=False, redraw_sidechains=False):
    """ returns a dataframe with the results of R-Group Decomposition

    >>> from rdkit import Chem
    >>> from rdkit.Chem import rdRGroupDecomposition
    >>> from rdkit.Chem import PandasTools
    >>> import pandas as pd
    >>> scaffold = Chem.MolFromSmiles('c1ccccn1')
    >>> mols = [Chem.MolFromSmiles(smi) for smi in 'c1c(F)cccn1 c1c(Cl)c(C)ccn1 c1c(O)cccn1 c1c(F)c(C)ccn1 c1cc(Cl)c(F)cn1'.split()]
    >>> groups,_ = rdRGroupDecomposition.RGroupDecompose([scaffold],mols,asSmiles=False,asRows=False)
    >>> df = PandasTools.RGroupDecompositionToFrame(groups,mols,include_core=False)
    >>> list(df.columns)
    ['Mol', 'R1', 'R2']
    >>> df = PandasTools.RGroupDecompositionToFrame(groups,mols,include_core=True)
    >>> list(df.columns)
    ['Mol', 'Core', 'R1', 'R2']
    >>> len(df)
    5
    >>> df.columns() # doctest: +SKIP
    <class 'pandas*...*DataFrame'>
    RangeIndex: 5 entries, 0 to 4
    Data columns (total 4 columns):
    Mol     5 non-null object
    Core    5 non-null object
    R1      5 non-null object
    R2      5 non-null object
    dtypes: object(4)
    memory usage: *...*

    """
    groups = groups.copy()
    cols = ['Mol'] + list(groups.keys())
    if redraw_sidechains:
      from rdkit.Chem import rdDepictor
      for k, vl in groups.items():
        if k == 'Core':
          continue
        for i, v in enumerate(vl):
          vl[i] = Chem.RemoveHs(v)
          rdDepictor.Compute2DCoords(vl[i])

    if not include_core:
      cols.remove('Core')
      del groups['Core']
    groups['Mol'] = mols
    frame = pd.DataFrame(groups, columns=cols)
    ChangeMoleculeRendering(frame)
    return frame


def _get_image(x):
  """displayhook function for PNG data"""
  return b64encode(x).decode('ascii')


try:
  from rdkit.Avalon import pyAvalonTools as pyAvalonTools

  # Calculate the Avalon fingerprint


  def _fingerprinter(x, y):
    return pyAvalonTools.GetAvalonFP(x, isQuery=y, bitFlags=pyAvalonTools.avalonSSSBits)
except ImportError:
  # Calculate fingerprint using SMARTS patterns
  def _fingerprinter(x, y):
    return Chem.PatternFingerprint(x, fpSize=2048)


def _MolPlusFingerprint(m):
  '''Precomputes fingerprints and stores results in molecule objects to accelerate
       substructure matching
    '''
  if m is not None:
    m._substructfp = _fingerprinter(m, False)
  return m


def AddMoleculeColumnToFrame(frame, smilesCol='Smiles', molCol='ROMol', includeFingerprints=False):
  '''Converts the molecules contains in "smilesCol" to RDKit molecules and appends them to the
    dataframe "frame" using the specified column name.
    If desired, a fingerprint can be computed and stored with the molecule objects to accelerate
    substructure matching
    '''
  if not includeFingerprints:
    frame[molCol] = frame[smilesCol].map(Chem.MolFromSmiles)
  else:
    frame[molCol] = frame[smilesCol].map(
      lambda smiles: _MolPlusFingerprint(Chem.MolFromSmiles(smiles)))
  ChangeMoleculeRendering(frame)


def ChangeMoleculeRendering(frame=None, renderer='image'):
  '''Allows to change the rendering of the molecules between image and string
    representations.
    This serves two purposes: First it allows to avoid the generation of images if this is
    not desired and, secondly, it allows to enable image rendering for newly created dataframe
    that already contains molecules, without having to rerun the time-consuming
    AddMoleculeColumnToFrame. Note: this behaviour is, because some pandas methods, e.g. head()
    returns a new dataframe instance that uses the default pandas rendering (thus not drawing
    images for molecules) instead of the monkey-patched one.
    '''
  # if frame is None this is is a no-op; we keep just for backwards compatibility
  if frame is None:
    return
  try:
    PandasPatcher.changeMoleculeRendering(frame, renderer)
  except NameError:
    log.warning("Failed to patch pandas - unable to change molecule rendering")


def WriteSDF(df, out, molColName='ROMol', idName=None, properties=None, allNumeric=False,
             forceV3000=False):
  '''Write an SD file for the molecules in the dataframe. Dataframe columns can be exported as
    SDF tags if specified in the "properties" list. "properties=list(df.columns)" would export
    all columns.
    The "allNumeric" flag allows to automatically include all numeric columns in the output.
    User has to make sure that correct data type is assigned to column.
    "idName" can be used to select a column to serve as molecule title. It can be set to
    "RowID" to use the dataframe row key as title.
    '''
  close = None
  if isinstance(out, str):
    if out.lower()[-3:] == ".gz":
      import gzip
      out = gzip.open(out, "wt")
      close = out.close

  writer = SDWriter(out)
  if forceV3000:
    writer.SetForceV3000(True)

  if properties is None:
    properties = []
  else:
    properties = list(properties)
  if allNumeric:
    properties.extend([
      dt for dt in df.dtypes.keys()
      if (np.issubdtype(df.dtypes[dt], np.floating) or np.issubdtype(df.dtypes[dt], np.integer))
    ])

  if molColName in properties:
    properties.remove(molColName)
  if idName in properties:
    properties.remove(idName)
  writer.SetProps(properties)
  for row in df.iterrows():
    # make a local copy I can modify
    mol = Chem.Mol(row[1][molColName])

    if idName is not None:
      if idName == 'RowID':
        mol.SetProp('_Name', str(row[0]))
      else:
        mol.SetProp('_Name', str(row[1][idName]))
    for p in properties:
      cell_value = row[1][p]
      # Make sure float does not get formatted in E notation
      if np.issubdtype(type(cell_value), np.floating):
        s = '{:f}'.format(cell_value).rstrip("0")  # "f" will show 7.0 as 7.00000
        if s[-1] == ".":
          s += "0"  # put the "0" back on if it's something like "7."
        mol.SetProp(p, s)
      else:
        mol.SetProp(p, str(cell_value))
    writer.write(mol)
  writer.close()
  if close is not None:
    close()


_saltRemover = None


def RemoveSaltsFromFrame(frame, molCol='ROMol'):
  '''
    Removes salts from mols in pandas DataFrame's ROMol column
    '''
  global _saltRemover
  if _saltRemover is None:
    from rdkit.Chem import SaltRemover
    _saltRemover = SaltRemover.SaltRemover()
  frame[molCol] = frame.apply(lambda x: _saltRemover.StripMol(x[molCol]), axis=1)


def SaveSMILESFromFrame(frame, outFile, molCol='ROMol', NamesCol='', isomericSmiles=False):
  '''
    Saves smi file. SMILES are generated from column with RDKit molecules. Column
    with names is optional.
    '''
  w = Chem.SmilesWriter(outFile, isomericSmiles=isomericSmiles)
  if NamesCol != '':
    for m, n in zip(frame[molCol], (str(c) for c in frame[NamesCol])):
      m.SetProp('_Name', n)
      w.write(m)
    w.close()
  else:
    for m in frame[molCol]:
      w.write(m)
    w.close()


def SaveXlsxFromFrame(frame, outFile, molCol='ROMol', size=(300, 300), formats=None):
  """
      Saves pandas DataFrame as a xlsx file with embedded images.
      molCol can be either a single column label or a list of column labels.
      It maps numpy data types to excel cell types:
      int, float -> number
      datetime -> datetime
      object -> string (limited to 32k character - xlsx limitations)

      The formats parameter can be optionally set to a dict of XlsxWriter
      formats (https://xlsxwriter.readthedocs.io/format.html#format), e.g.:
      {
        'write_string':  {'text_wrap': True}
      }
      Currently supported keys for the formats dict are:
      'write_string', 'write_number', 'write_datetime'.

      Cells with compound images are a bit larger than images due to excel.
      Column width weirdness explained (from xlsxwriter docs):
      The width corresponds to the column width value that is specified in Excel.
      It is approximately equal to the length of a string in the default font of Calibri 11.
      Unfortunately, there is no way to specify "AutoFit" for a column in the Excel file format.
      This feature is only available at runtime from within Excel.
      """

  import xlsxwriter  # don't want to make this a RDKit dependency

  cols = list(frame.columns)
  if isinstance(molCol, str):
    molCol = [molCol]
  molCol = list(set(molCol))
  dataTypes = dict(frame.dtypes)
  molCol_indices = [cols.index(mc) for mc in molCol]

  workbook = xlsxwriter.Workbook(outFile)  # New workbook
  cell_formats = {}
  formats = formats or {}
  for key in ['write_string', 'write_number', 'write_datetime']:
    format = formats.get(key, None)
    if format is not None:
      format = workbook.add_format(format)
    cell_formats[key] = format
  worksheet = workbook.add_worksheet()  # New work sheet

  # Write first row with column names
  for col_idx, col in enumerate(cols):
    worksheet.write_string(0, col_idx, col)

  for row_idx, (_, row) in enumerate(frame.iterrows()):
    row_idx_actual = row_idx + 1

    worksheet.set_row(row_idx_actual, height=size[1])  # looks like height is not in px?

    for col_idx, col in enumerate(cols):
      if col_idx in molCol_indices:
        image_data = BytesIO()
        m = row[col]
        img = Draw.MolToImage(m if isinstance(m, Chem.Mol) else Chem.Mol(), size=size,
                              options=drawOptions)
        img.save(image_data, format='PNG')
        worksheet.insert_image(row_idx_actual, col_idx, "f", {'image_data': image_data})
        worksheet.set_column(col_idx, col_idx,
                             width=size[0] / 6.)  # looks like height is not in px?
      elif str(dataTypes[col]) == "object":
        # string length is limited in xlsx
        worksheet.write_string(row_idx_actual, col_idx,
                               str(row[col])[:32000], cell_formats['write_string'])
      elif ('float' in str(dataTypes[col])) or ('int' in str(dataTypes[col])):
        if (row[col] != np.nan) or (row[col] != np.inf):
          worksheet.write_number(row_idx_actual, col_idx, row[col], cell_formats['write_number'])
      elif 'datetime' in str(dataTypes[col]):
        worksheet.write_datetime(row_idx_actual, col_idx, row[col], cell_formats['write_datetime'])

  workbook.close()
  image_data.close()


def FrameToGridImage(frame, column='ROMol', legendsCol=None, **kwargs):
  '''
    Draw grid image of mols in pandas DataFrame.
    '''
  if legendsCol:
    if legendsCol == frame.index.name:
      kwargs['legends'] = [str(c) for c in frame.index]
    else:
      kwargs['legends'] = [str(c) for c in frame[legendsCol]]
  return Draw.MolsToGridImage(list(frame[column]), **kwargs)


def AddMurckoToFrame(frame, molCol='ROMol', MurckoCol='Murcko_SMILES', Generic=False):
  '''
    Adds column with SMILES of Murcko scaffolds to pandas DataFrame.

    Generic set to true results in SMILES of generic framework.
    '''
  if Generic:

    def func(x):
      return Chem.MolToSmiles(
        MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(x[molCol])))
  else:

    def func(x):
      return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(x[molCol]))

  frame[MurckoCol] = frame.apply(func, axis=1)


def AlignMol(mol, scaffold):
  """
    Aligns mol (RDKit mol object) to scaffold (SMILES string)
    """
  scaffold = Chem.MolFromSmiles(scaffold)
  AllChem.Compute2DCoords(scaffold)
  AllChem.GenerateDepictionMatching2DStructure(mol, scaffold)
  return mol


def AlignToScaffold(frame, molCol='ROMol', scaffoldCol='Murcko_SMILES'):
  '''
    Aligns molecules in molCol to scaffolds in scaffoldCol
    '''
  frame[molCol] = frame.apply(lambda x: AlignMol(x[molCol], x[scaffoldCol]), axis=1)


# ==========================================================================================
# Monkey patching RDkit functionality
def InstallPandasTools():
  """ Monkey patch an RDKit method of Chem.Mol and pandas """
  try:
    PandasPatcher.patchPandas()
  except NameError:
    pass
  if 'Chem.Mol.__ge__' not in _originalSettings:
    _originalSettings['Chem.Mol.__ge__'] = rdchem.Mol.__ge__
    rdchem.Mol.__ge__ = _molge


def UninstallPandasTools():
  """ Unpatch an RDKit method of Chem.Mol and pandas """
  try:
    PandasPatcher.unpatchPandas()
  except NameError:
    pass
  if 'Chem.Mol.__ge__' in _originalSettings:
    Chem.Mol.__ge__ = _originalSettings['Chem.Mol.__ge__']
  _originalSettings.clear()


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE,
                              verbose=verbose)
  if (failed):
    sys.exit(failed)


InstallPandasTools()

try:
  import xlsxwriter
except ImportError:
  xlsxwriter = None
import unittest


class TestCase(unittest.TestCase):

  @unittest.skipIf(xlsxwriter is None or pd is None, 'pandas/xlsxwriter not installed')
  def testGithub1507(self):
    import os
    from rdkit import RDConfig
    sdfFile = os.path.join(RDConfig.RDDataDir, 'NCI/first_200.props.sdf')
    frame = LoadSDF(sdfFile)
    SaveXlsxFromFrame(frame, 'foo.xlsx', formats={'write_string': {'text_wrap': True}})

  @unittest.skipIf(pd is None, 'pandas not installed')
  def testGithub3701(self):
    ' problem with update to pandas v1.2.0 '
    df = pd.DataFrame({"name": ["ethanol", "furan"], "smiles": ["CCO", "c1ccoc1"]})
    AddMoleculeColumnToFrame(df, 'smiles', 'molecule')
    self.assertEqual(len(df.molecule), 2)


if __name__ == '__main__':  # pragma: nocover
  try:
    import pandas as pd

    def _getPandasVersion():
      """ Get the pandas version as a tuple """
      import re
      try:
        v = pd.__version__
      except AttributeError:
        v = pd.version.version
      v = re.split(r'[^0-9,.]', v)[0].split('.')
      return tuple(int(vi) for vi in v)
  except ImportError:
    pd = None

  class TestCase(unittest.TestCase):

    @unittest.skipIf(xlsxwriter is None or pd is None, 'pandas/xlsxwriter not installed')
    def testGithub1507(self):
      import os

      from rdkit import RDConfig
      sdfFile = os.path.join(RDConfig.RDDataDir, 'NCI/first_200.props.sdf')
      frame = LoadSDF(sdfFile)
      SaveXlsxFromFrame(frame, 'foo.xlsx', formats={'write_string': {'text_wrap': True}})

    @unittest.skipIf(pd is None, 'pandas not installed')
    def testGithub3701(self):
      ' problem with update to pandas v1.2.0 '
      df = pd.DataFrame({"name": ["ethanol", "furan"], "smiles": ["CCO", "c1ccoc1"]})
      AddMoleculeColumnToFrame(df, 'smiles', 'molecule')
      self.assertEqual(len(df.molecule), 2)

  runTests = True
  try:
    pd
    runTests = not _getPandasVersion() < (0, 19)
  except NameError:
    runTests = False

  if runTests:
    _runDoctests()
    unittest.main()
  else:
    log.warning("pandas installation >=0.19 not found, skipping tests")
