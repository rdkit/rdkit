import RDLogger as logging
logger = logging.logger()
logger.setLevel(logging.INFO)
import Chem
from Chem import Crippen
from Chem import AllChem
from Chem.ChemUtils.AlignDepict import AlignDepict

_version="0.2.2"
_greet="This is TemplateExpand version %s"%_version

_usage="""
Usage: TemplateExpand [options] template <sidechains>

 Unless otherwise indicated, the template and sidechains are assumed to be
   Smiles

 Each sidechain entry should be:
   [-r] SMARTS filename
     The SMARTS pattern is used to recognize the attachment point,
     if the -r argument is not provided, then atoms matching the pattern
     will be removed from the sidechains.
   or
   -n filename
     where the attachment atom is the first atom in each molecule
   The filename provides the list of potential sidechains.

 options:
   -o filename.sdf: provides the name of the output file, otherwise
                    stdout is used
   --sdf : sidechains should be in SD files
   --molTemplate: the template is a Mol file, a new depiction
                  will not be generated

"""
def Usage():
  import sys
  print _usage
  sys.exit(-1)

nDumped=0 
def _exploder(mol,depth,sidechains,outF,core,chainIndices):
  global nDumped
  ourChains = sidechains[depth]
  patt = '[%d*]'%(depth+1)
  patt = Chem.MolFromSmarts(patt)
  for i,chain in enumerate(ourChains):
    tchain = chainIndices[:]
    tchain.append(i+1)
    rs = Chem.ReplaceSubstructs(mol,patt,chain,replaceAll=True)
    if rs:
      r = rs[0]
      if depth<len(sidechains)-1:
        _exploder(r,depth+1,sidechains,outF,core,tchain)
      else:
        Chem.SanitizeMol(r)
        if r.HasSubstructMatch(core):
          AlignDepict(r,core)
        else:
          print >>sys.stderr,'>>> no match'
          AllChem.Compute2DCoords(r)
        Chem.Kekulize(r)
        r.SetProp("_Name","TemplateEnum: Mol_%d"%(nDumped+1))
        mb = Chem.MolToMolBlock(r)
        outF.write(mb)
        print >>outF,">  <seq_num>\n%d\n"%(nDumped+1)
        print >>outF,">  <reagent_indices>\n%s\n"%('_'.join([str(x) for x in tchain]))
        print >>outF,'$$$$'
        nDumped += 1
        if not nDumped%100:
          logger.info('Done %d molecules'%nDumped)
          
def Explode(template,sidechains,outF):
  chainIndices=[]
  core = Chem.DeleteSubstructs(template,Chem.MolFromSmiles('[*]'))
  _exploder(template,0,sidechains,outF,core,chainIndices)
                              

def ConstructSidechains(suppl,sma=None,replace=True):
  if sma:
    try:
      patt = Chem.MolFromSmarts(sma)
    except:
      logger.error('could not construct pattern from smarts: %s'%sma,
                   exc_info=True)
      return None
  else:
    patt = None

  if replace:
    replacement = Chem.MolFromSmiles('[*]')

  res = []
  for mol in suppl:
    if not mol:
      continue
    if patt:
      if replace:
        tmp = list(Chem.ReplaceSubstructs(mol,patt,replacement))
        for i,entry in enumerate(tmp):
          match = entry.GetSubstructMatch(replacement)
          if match:
            idx = match[0]
            smi = Chem.MolToSmiles(entry,rootedAtAtom=idx)
            entry = Chem.MolFromSmiles(smi)
            # entry now has [*] as atom 0 and the neighbor
            # as atom 1. Cleave the [*]:
            entry = Chem.DeleteSubstructs(entry,replacement)
            # now we have a molecule with the atom to be joined
            # in position zero; Keep that:
            tmp[i] = entry
      else:
        # no replacement, use the pattern to reorder
        # atoms only:
        match = mol.GetSubstructMatch(patt)
        if match:
          smi = Chem.MolToSmiles(mol,rootedAtAtom=match[0])
          entry = Chem.MolFromSmiles(smi)
          # now we have a molecule with the atom to be joined
          # in position zero; Keep that:
          tmp = [entry]
        else:
          tmp = None
    else:
      tmp = [mol]
    if tmp:
      res.extend(tmp)
  return res

if __name__=='__main__':
  import getopt,sys
  print >>sys.stderr,_greet
  
  try:
    args,extras = getopt.getopt(sys.argv[1:],'o:h',[
      'sdf',
      'molTemplate',
      'force'
      ])
  except:
    import traceback
    traceback.print_exc()
    Usage()

  if len(extras)<3:
    Usage()
  
  tooLong=1000
  sdLigands=False
  molTemplate=False
  redrawTemplate=True
  outF=None
  forceIt=False
  for arg,val in args:
    if arg=='-o':
      outF=val
    elif arg=='--sdf':
      sdLigands=True
    elif arg=='--molTemplate':
      molTemplate=True
    elif arg=='--force':
      forceIt=True
    elif arg=='-h':
      Usage()
      sys.exit(0)

  if molTemplate:
    try:
      template = Chem.MolFromMolFile(extras[0])
    except:
      logger.error('Could not construct template from mol file: %s'%extras[0],
                   exc_info=True)
      sys.exit(1)
    if redrawTemplate:
      AllChem.Compute2DCoords(template)
  else:
    try:
      template = Chem.MolFromSmiles(extras[0])
      AllChem.Compute2DCoords(template)
    except:
      logger.error('Could not construct template from smiles: %s'%extras[0],
                   exc_info=True)
      sys.exit(1)

  sidechains = []
  pos = 1
  while pos<len(extras):
    if extras[pos]=='-r':
      replaceIt=False
      pos += 1
    else:
      replaceIt=True
    if extras[pos]=='-n':
      sma = None
    else:
      sma = extras[pos]
    pos += 1
    try:
      dat = extras[pos]
    except IndexError:
      logger.error('missing a sidechain filename')
      sys.exit(-1)
    pos += 1
    if sdLigands:
      try:
        suppl = Chem.SDMolSupplier(dat)
      except:
        logger.error('could not construct supplier from SD file: %s'%dat,
                     exc_info=True)
        suppl = []
    else:
      try:
        suppl = Chem.SmilesMolSupplier(dat)
      except:
        logger.error('could not construct supplier from smiles file: %s'%dat,
                     exc_info=True)
        suppl = []
      suppl = [x for x in suppl]
    chains = ConstructSidechains(suppl,sma=sma,replace=replaceIt)
    if chains:
      sidechains.append(chains)

  count = 1
  for chain in sidechains:
    count *= len(chain)
  if not sidechains or not count:
    print >>sys.stderr,"No molecules to be generated."
    sys.exit(0)

                                 

  if not forceIt and count>tooLong:
    print >>sys.stderr,"This will generate %d molecules."%count
    ans = raw_input("Continue anyway? [no]")
    if ans not in ('y','yes','Y','YES'):
      sys.exit(0)

  if outF:
    try:
      outF = file(outF,'w+')
    except IOError:
      logger.error('could not open file %s for writing'%(outF),
                   exc_info=True)
  else:
    outF = sys.stdout


  Explode(template,sidechains,outF)
  



  
