# $Id: TemplateExpand.py 1053 2008-07-30 12:03:29Z landrgr1 $
#
#  Created by Greg Landrum August, 2006
#
#

from rdkit import RDLogger as logging

logger = logging.logger()
logger.setLevel(logging.INFO)
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen
from rdkit.Chem.ChemUtils.AlignDepict import AlignDepict

_version = "0.8.0"
_greet = "This is TemplateExpand version %s" % _version

_usage = """
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
   -o filename.sdf:      provides the name of the output file, otherwise
                         stdout is used

   --sdf :               expect the sidechains to be in SD files

   --moltemplate:        the template(s) are in a mol/SD file, new depiction(s)
                         will not be generated unless the --redraw argument is also
                         provided

   --smilesFileTemplate: extract the template(s) from a SMILES file instead of 
                         expecting SMILES on the command line.

   --redraw:             generate a new depiction for the molecular template(s)

   --useall:
     or
   --useallmatches:      generate a product for each possible match of the attachment
                         pattern to each sidechain. If this is not provided, the first
                         match (not canonically defined) will be used.

   --force:              by default, the program prompts the user if the library is 
                         going to contain more than 1000 compounds. This argument 
                         disables the prompt.
   
   --templateSmarts="smarts":  provides a space-delimited list containing the SMARTS 
                               patterns to be used to recognize attachment points in
                               the template
             
   --autoNames:          when set this toggle causes the resulting compounds to be named
                         based on there sequence id in the file, e.g. 
                         "TemplateEnum: Mol_1", "TemplateEnum: Mol_2", etc.
                         otherwise the names of the template and building blocks (from
                         the input files) will be combined to form a name for each
                         product molecule.

   --3D :                Generate 3d coordinates for the product molecules instead of 2d coordinates,
                         requires the --moltemplate option

   --tether :            refine the 3d conformations using a tethered minimization


"""


def Usage():
  print(_usage, file=sys.stderr)
  sys.exit(-1)


nDumped = 0


def _exploder(mol, depth, sidechains, core, chainIndices, autoNames=True, templateName='',
              resetCounter=True, do3D=False, useTethers=False):
  global nDumped
  if resetCounter:
    nDumped = 0
  ourChains = sidechains[depth]
  patt = Chem.MolFromSmiles(f'[{depth + 1}*]')
  for i, (chainIdx, chain) in enumerate(ourChains):
    tchain = chainIndices[:]
    tchain.append((i, chainIdx))
    rs = Chem.ReplaceSubstructs(mol, patt, chain, replaceAll=True)
    if rs:
      r = rs[0]
      if depth < len(sidechains) - 1:
        for entry in _exploder(r, depth + 1, sidechains, core, tchain, autoNames=autoNames,
                               templateName=templateName, resetCounter=0, do3D=do3D,
                               useTethers=useTethers):
          yield entry
      else:
        try:
          Chem.SanitizeMol(r)
        except ValueError:
          import traceback
          traceback.print_exc()
          continue
        if not do3D:
          if r.HasSubstructMatch(core):
            try:
              AlignDepict(r, core)
            except Exception:
              import traceback
              traceback.print_exc()
              print(Chem.MolToSmiles(r), file=sys.stderr)
          else:
            print('>>> no match', file=sys.stderr)
            AllChem.Compute2DCoords(r)
        else:
          r = Chem.AddHs(r)
          AllChem.ConstrainedEmbed(r, core, useTethers)
        Chem.Kekulize(r)
        if autoNames:
          tName = "TemplateEnum: Mol_%d" % (nDumped + 1)
        else:
          tName = templateName
          for bbI, bb in enumerate(tchain):
            bbMol = sidechains[bbI][bb[0]][1]
            if bbMol.HasProp('_Name'):
              bbNm = bbMol.GetProp('_Name')
            else:
              bbNm = str(bb[1])
            tName += '_' + bbNm

        r.SetProp("_Name", tName)
        r.SetProp('seq_num', str(nDumped + 1))
        r.SetProp('reagent_indices', '_'.join([str(x[1]) for x in tchain]))
        for bbI, bb in enumerate(tchain):
          bbMol = sidechains[bbI][bb[0]][1]
          if bbMol.HasProp('_Name'):
            bbNm = bbMol.GetProp('_Name')
          else:
            bbNm = str(bb[1])
          r.SetProp('building_block_%d' % (bbI + 1), bbNm)
          r.SetIntProp('_idx_building_block_%d' % (bbI + 1), bb[1])
          for propN in bbMol.GetPropNames():
            r.SetProp('building_block_%d_%s' % (bbI + 1, propN), bbMol.GetProp(propN))
        nDumped += 1
        if not nDumped % 100:
          logger.info('Done %d molecules' % nDumped)
        yield r


def Explode(template, sidechains, outF, autoNames=True, do3D=False, useTethers=False):
  chainIndices = []
  core = Chem.DeleteSubstructs(template, Chem.MolFromSmiles('[*]'))
  try:
    templateName = template.GetProp('_Name')
  except KeyError:
    templateName = "template"

  for mol in _exploder(template, 0, sidechains, core, chainIndices, autoNames=autoNames,
                       templateName=templateName, do3D=do3D, useTethers=useTethers):
    outF.write(Chem.MolToMolBlock(mol))
    for pN in mol.GetPropNames():
      print('>  <%s>\n%s\n' % (pN, mol.GetProp(pN)), file=outF)
    print('$$$$', file=outF)


def MoveDummyNeighborsToBeginning(mol, useAll=False):
  dummyPatt = Chem.MolFromSmiles('[*]')
  matches = mol.GetSubstructMatches(dummyPatt)
  res = []
  for match in matches:
    smi = Chem.MolToSmiles(mol, rootedAtAtom=match[0])
    entry = Chem.MolFromSmiles(smi)
    # entry now has [*] as atom 0 and the neighbor
    # as atom 1. Cleave the [*]:
    entry = Chem.DeleteSubstructs(entry, dummyPatt)
    for propN in mol.GetPropNames():
      entry.SetProp(propN, mol.GetProp(propN))

    # now we have a molecule with the atom to be joined
    # in position zero; Keep that:
    res.append(entry)
    if not useAll:
      break
  return res


def ConstructSidechains(suppl, sma=None, replace=True, useAll=False):
  if sma:
    patt = Chem.MolFromSmarts(sma)
    if patt is None:
      logger.error('could not construct pattern from smarts: %s' % sma, exc_info=True)
      return None
  else:
    patt = None

  if replace:
    replacement = Chem.MolFromSmiles('[*]')

  res = []
  for idx, mol in enumerate(suppl):
    if not mol:
      continue
    if patt:
      if not mol.HasSubstructMatch(patt):
        logger.warning(
          'The substructure pattern did not match sidechain %d. This may result in errors.' %
          (idx + 1))
      if replace:
        tmp = list(Chem.ReplaceSubstructs(mol, patt, replacement))
        if not useAll:
          tmp = [tmp[0]]
        for i, entry in enumerate(tmp):
          entry = MoveDummyNeighborsToBeginning(entry)
          if not entry:
            continue
          entry = entry[0]

          for propN in mol.GetPropNames():
            entry.SetProp(propN, mol.GetProp(propN))
          # now we have a molecule with the atom to be joined
          # in position zero; Keep that:
          tmp[i] = (idx + 1, entry)
      else:
        # no replacement, use the pattern to reorder
        # atoms only:
        matches = mol.GetSubstructMatches(patt)
        if matches:
          tmp = [0] * len(matches)
          for i, match in enumerate(matches):
            smi = Chem.MolToSmiles(mol, rootedAtAtom=match[0])
            entry = Chem.MolFromSmiles(smi)
            for propN in mol.GetPropNames():
              entry.SetProp(propN, mol.GetProp(propN))

            # now we have a molecule with the atom to be joined
            # in position zero; Keep that:
            tmp[i] = (idx + 1, entry)

        else:
          tmp = None
    else:
      tmp = [(idx + 1, mol)]

    if tmp:
      res.extend(tmp)

  return res


if __name__ == '__main__':
  import getopt
  print(_greet, file=sys.stderr)

  try:
    args, extras = getopt.getopt(sys.argv[1:], 'o:h', [
      'sdf',
      'moltemplate',
      'molTemplate',
      'smilesFileTemplate',
      'templateSmarts=',
      'redraw',
      'force',
      'useall',
      'useallmatches',
      'autoNames',
      '3D',
      '3d',
      'tethers',
      'tether',
    ])
  except Exception:
    import traceback
    traceback.print_exc()
    Usage()

  if len(extras) < 3:
    Usage()

  tooLong = 1000
  sdLigands = False
  molTemplate = False
  redrawTemplate = False
  outF = None
  forceIt = False
  useAll = False
  templateSmarts = []
  smilesFileTemplate = False
  autoNames = False
  do3D = False
  useTethers = False
  for arg, val in args:
    if arg == '-o':
      outF = val
    elif arg == '--sdf':
      sdLigands = True
    elif arg in ('--moltemplate', '--molTemplate'):
      molTemplate = True
    elif arg == '--smilesFileTemplate':
      smilesFileTemplate = True
    elif arg == '--templateSmarts':
      templateSmarts = val
    elif arg == '--redraw':
      redrawTemplate = True
    elif arg == '--force':
      forceIt = True
    elif arg == '--autoNames':
      autoNames = True
    elif arg in ('--useall', '--useallmatches'):
      useAll = True
    elif arg in ('--3D', '--3d'):
      do3D = True
    elif arg in ('--tethers', '--tether'):
      useTethers = True
    elif arg == '-h':
      Usage()
      sys.exit(0)

  if do3D:
    if not molTemplate:
      raise ValueError('the --3D option is only useable in combination with --moltemplate')
    if redrawTemplate:
      logger.warning(
        '--redrawTemplate does not make sense in combination with --molTemplate. removing it')
      redrawTemplate = False

  if templateSmarts:
    splitL = templateSmarts.split(' ')
    templateSmarts = []
    for i, sma in enumerate(splitL):
      patt = Chem.MolFromSmarts(sma)
      if not patt:
        raise ValueError('could not convert smarts "%s" to a query' % sma)
      if i >= 4:
        i += 1
      replace = Chem.MolFromSmiles('[%d*]' % (i + 1))
      templateSmarts.append((patt, replace))

  if molTemplate:
    removeHs = not do3D
    try:
      s = Chem.SDMolSupplier(extras[0], removeHs=removeHs)
      templates = [x for x in s]
    except Exception:
      logger.error('Could not construct templates from input file: %s' % extras[0], exc_info=True)
      sys.exit(1)
    if redrawTemplate:
      for template in templates:
        AllChem.Compute2DCoords(template)
  else:
    if not smilesFileTemplate:
      try:
        templates = [Chem.MolFromSmiles(extras[0])]
      except Exception:
        logger.error('Could not construct template from smiles: %s' % extras[0], exc_info=True)
        sys.exit(1)
    else:
      try:
        s = Chem.SmilesMolSupplier(extras[0], titleLine=False)
        templates = [x for x in s]
      except Exception:
        logger.error('Could not construct templates from input file: %s' % extras[0], exc_info=True)
        sys.exit(1)
    for template in templates:
      AllChem.Compute2DCoords(template)
  if templateSmarts:
    finalTs = []
    for i, template in enumerate(templates):
      for j, (patt, replace) in enumerate(templateSmarts):
        if not template.HasSubstructMatch(patt):
          logger.error('template %d did not match sidechain pattern %d, skipping it' %
                       (i + 1, j + 1))
          template = None
          break
        template = Chem.ReplaceSubstructs(template, patt, replace)[0]
      if template:
        Chem.SanitizeMol(template)
        finalTs.append(template)
    templates = finalTs

  sidechains = []
  pos = 1
  while pos < len(extras):
    if extras[pos] == '-r':
      replaceIt = False
      pos += 1
    else:
      replaceIt = True
    if extras[pos] == '-n':
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
      except Exception:
        logger.error('could not construct supplier from SD file: %s' % dat, exc_info=True)
        suppl = []
    else:
      tmpF = file(dat, 'r')
      inL = tmpF.readline()
      if len(inL.split(' ')) < 2:
        nmCol = -1
      else:
        nmCol = 1
      try:
        suppl = Chem.SmilesMolSupplier(dat, nameColumn=nmCol)
      except Exception:
        logger.error('could not construct supplier from smiles file: %s' % dat, exc_info=True)
        suppl = []
      suppl = [x for x in suppl]
    chains = ConstructSidechains(suppl, sma=sma, replace=replaceIt, useAll=useAll)
    if chains:
      sidechains.append(chains)
  count = 1
  for chain in sidechains:
    count *= len(chain)
  count *= len(templates)
  if not sidechains or not count:
    print("No molecules to be generated.", file=sys.stderr)
    sys.exit(0)

  if not forceIt and count > tooLong:
    print("This will generate %d molecules." % count, file=sys.stderr)
    print("Continue anyway? [no] ", file=sys.stderr, end='')
    sys.stderr.flush()
    ans = sys.stdin.readline().strip()
    if ans not in ('y', 'yes', 'Y', 'YES'):
      sys.exit(0)

  if outF and outF != "-":
    try:
      outF = file(outF, 'w+')
    except IOError:
      logger.error('could not open file %s for writing' % (outF), exc_info=True)
  else:
    outF = sys.stdout

  for template in templates:
    Explode(template, sidechains, outF, autoNames=autoNames, do3D=do3D, useTethers=useTethers)
