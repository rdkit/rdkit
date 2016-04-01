# $Id: run.py 3578 2010-05-21 15:16:34Z landrgr1 $
#
# Created by Greg Landrum, May 2010
""" contains the controlling code for the R-group decomposition web service

"""
import os,sys,base64,time,subprocess
from mod_python import apache
from cStringIO import StringIO
import Extract
from rdkit import Chem


def Decompose(req):
    options,foo = Extract.parser.parse_args([])
    if req.form.has_key('core'):
        options.core=str(req.form['core'])
    if req.form.has_key('symmetrize'):
        options.symmetrize=True
    if req.form.has_key('labelledCores'):
        options.labelledCores=True
    if req.form.has_key('requireLabels'):
        options.requireLabels=True
    if req.form.has_key('coreFormat'):
        options.coreFormat=str(req.form['coreFormat'])
    if req.form.has_key('dataFormat'):
        options.dataFormat=str(req.form['dataFormat'])
    if req.form.has_key('outputFormat'):
        options.outputFormat=str(req.form['outputFormat'])

    if req.form.has_key('filename') and req.form['filename'].filename:
        molData = req.form['filename'].value
    else:
        molData = str(req.form.get('molData',''))

    if not options.core:
        raise ValueError,'no core definition'
    if options.coreFormat=='smarts':
      options.cores=[Chem.MolFromSmarts(x) for x in options.core.split('.')]
    elif options.coreFormat=='smiles':
      options.cores=[Chem.MolFromSmiles(x) for x in options.core.split('.')]
    else:
      options.cores=[]
    while options.cores.count(None): options.cores.remove(None)
    if not options.cores: 
        raise ValueError,'no valid cores found'
    
    from rdkit.Novartis.Plumbing.Input.MolLoader import LoadMols
    from rdkit.Novartis.Plumbing.Types.Mol import getRDKMol
    mols = LoadMols(molData)
    mols = [getRDKMol(x) for x in mols]
    sio = StringIO()
    options.outFile=sio
    mols,res = Extract.RunDecomposition(mols,options)
    Extract.CreateOutput(mols,res,options)
    return sio.getvalue()



    smis=req.form.get('smiles','')
    smis = smis.split(',')
    mols = [Chem.MolFromSmiles(x) for x in smis]
    while mols.count(None): mols.remove(None)


    minTani=float(req.form.get('MinTanimotoSim',0.6))
    minClusterSize=int(req.form.get('MinClusterSize',4))
    minMCSSize=int(req.form.get('MinMCSSize',6))
    ringMatches=not int(req.form.get('AllowRingChainMatches','0'))
    completeRings=not int(req.form.get('AllowIncompleteRings','0'))
    atomTyping=req.form.get('AtomTyping','basic')        
    timeoutLimit=int(req.form.get('Timeout',30))

    page="""%(res)s"""%(locals())
    return page

def Help(req):
    res = """
Available methods:
- Decompose : this is callable as:
http://www.dv.cadd.eu.novartis.intra/CaddServices/SideChains/run.py/Decompose


Until real documentation arrives, here's a sample HTML form that uses
this. Explanations of some of the parameters are below:

<form action="run.py/Decompose" method="post" enctype="multipart/form-data" id="mainForm">
<input type="submit" id="submit" />
<br />core: <input type="text" name="core" id="core" />
<br />symmetrize?: <input type="checkbox" name="symmetrize" id="symmetrize" checked=1/>
<br />labelled cores?: <input type="checkbox" name="labelledCores" id="labelledCores" checked=1 />
      require labels?: <input type="checkbox" name="requireLabels" id="requireLabels" checked=1 />
<p/>provide either a:
<br />Filename: <input type="file" name="filename" id="filename" />
<br />or some data: <br /><textarea name="molData" id="molData" cols="80" rows="20"></textarea>


<input type="hidden" name="coreFormat" id="coreFormat" value="smarts" />
<input type="hidden" name="dataFormat" id="dataFormat" value="other" />
<input type="hidden" name="outputFormat" id="outputFormat" value="csv" />
</form>

Explanations of some of the parameters:
- Molecule data can be provided as either SDF or delimited text. The
tool tries to be smart about figuring the format out, but if you try
you can easily break this.
- if provided, the "filename" element takes precedence over the
"molData" element.
- In either SMARTS or SMILES, R groups should be marked with dummy
atoms (symbol "*")
- If they are provided as SMARTS (the default), R labels in the cores
can be provided using atom mapping, i.e.[*:1]
- If you want to provide multiple cores as SMILES or SMARTS, separate
the individual cores with dots. i.e. c1cnc([*:1])cc1.c1cccn1[*:1]
- supported output formats are "sdf" and "csv"

"""
    return res
