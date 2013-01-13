
RDKit Cookbook
%%%%%%%%%%%%%%



What is this?
*************

This document provides examples of how to carry out particular tasks 
using the RDKit functionality from Python. The contents have been
contributed by the RDKit community.

If you find mistakes, or have suggestions for improvements, please
either fix them yourselves in the source document (the .rst file) or
send them to the mailing list: rdkit-discuss@lists.sourceforge.net 
(you will need to subscribe first)


Miscellaneous Topics
********************

Using a different aromaticity model
-----------------------------------

By default, the RDKit applies its own model of aromaticity (explained
in the RDKit Theory Book) when it reads in molecules. It is, however,
fairly easy to override this and use your own aromaticity model.

The easiest way to do this is it provide the molecules as SMILES with
the aromaticity set as you would prefer to have it. For example,
consider indole: 

.. image:: images/indole1.png
 
By default the RDKit considers both rings to be aromatic:

>>> from rdkit import Chem
>>> m = Chem.MolFromSmiles('N1C=Cc2ccccc12')
>>> m.GetSubstructMatches(Chem.MolFromSmarts('c'))
((1,), (2,), (3,), (4,), (5,), (6,), (7,), (8,))

If you'd prefer to treat the five-membered ring as aliphatic, which is
how the input SMILES is written, you just need to do a partial
sanitization that skips the kekulization and aromaticity perception
steps: 

>>> m2 = Chem.MolFromSmiles('N1C=Cc2ccccc12',sanitize=False)
>>> Chem.SanitizeMol(m2,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
>>> m2.GetSubstructMatches(Chem.MolFromSmarts('c'))
((3,), (4,), (5,), (6,), (7,), (8,))

It is, of course, also possible to write your own aromaticity
perception function, but that is beyond the scope of this document.


Manipulating Molecules
**********************

Cleaning up heterocycles
------------------------

Mailing list discussions:

*  http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01185.html
*  http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01162.html
*  http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01900.html   
*  http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01901.html   

The code:

.. testcode::

  """ sanifix4.py
   
    Contribution from James Davidson
  """
  from rdkit import Chem
  from rdkit.Chem import AllChem

  def _FragIndicesToMol(oMol,indices):
      em = Chem.EditableMol(Chem.Mol())

      newIndices={}
      for i,idx in enumerate(indices):
          em.AddAtom(oMol.GetAtomWithIdx(idx))
          newIndices[idx]=i

      for i,idx in enumerate(indices):
          at = oMol.GetAtomWithIdx(idx)
          for bond in at.GetBonds():
              if bond.GetBeginAtomIdx()==idx:
                  oidx = bond.GetEndAtomIdx()
              else:
                  oidx = bond.GetBeginAtomIdx()
              # make sure every bond only gets added once:
              if oidx<idx:
                  continue
              em.AddBond(newIndices[idx],newIndices[oidx],bond.GetBondType())
      res = em.GetMol()
      res.ClearComputedProps()
      Chem.GetSymmSSSR(res)
      res.UpdatePropertyCache(False)
      res._idxMap=newIndices
      return res

  def _recursivelyModifyNs(mol,matches,indices=None):
      if indices is None:
          indices=[]
      res=None
      while len(matches) and res is None:
          tIndices=indices[:]
          nextIdx = matches.pop(0)
          tIndices.append(nextIdx)
          nm = Chem.Mol(mol.ToBinary())
          nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
          nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
          cp = Chem.Mol(nm.ToBinary())
          try:
              Chem.SanitizeMol(cp)
          except ValueError:
              res,indices = _recursivelyModifyNs(nm,matches,indices=tIndices)
          else:
              indices=tIndices
              res=cp
      return res,indices

  def AdjustAromaticNs(m,nitrogenPattern='[n&D2&H0;r5,r6]'):
      """
         default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
         to fix: O=c1ccncc1
      """
      Chem.GetSymmSSSR(m)
      m.UpdatePropertyCache(False)

      # break non-ring bonds linking rings:
      em = Chem.EditableMol(m)
      linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
      plsFix=set()
      for a,b in linkers:
          em.RemoveBond(a,b)
          plsFix.add(a)
          plsFix.add(b)
      nm = em.GetMol()
      for at in plsFix:
          at=nm.GetAtomWithIdx(at)
          if at.GetIsAromatic() and at.GetAtomicNum()==7:
              at.SetNumExplicitHs(1)
              at.SetNoImplicit(True)

      # build molecules from the fragments:
      fragLists = Chem.GetMolFrags(nm)
      frags = [_FragIndicesToMol(nm,x) for x in fragLists]

      # loop through the fragments in turn and try to aromatize them:
      ok=True
      for i,frag in enumerate(frags):
          cp = Chem.Mol(frag.ToBinary())
          try:
              Chem.SanitizeMol(cp)
          except ValueError:
              matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
              lres,indices=_recursivelyModifyNs(frag,matches)
              if not lres:
                  #print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                  ok=False
                  break
              else:
                  revMap={}
                  for k,v in frag._idxMap.iteritems():
                      revMap[v]=k
                  for idx in indices:
                      oatom = m.GetAtomWithIdx(revMap[idx])
                      oatom.SetNoImplicit(True)
                      oatom.SetNumExplicitHs(1)
      if not ok:
          return None
      return m


Examples of using it:

.. testcode:: 

  smis= ('O=c1ccc2ccccc2n1',
         'Cc1nnnn1C',
         'CCc1ccc2nc(=O)c(cc2c1)Cc1nnnn1C1CCCCC1',
         'c1cnc2cc3ccnc3cc12',
         'c1cc2cc3ccnc3cc2n1',
         'O=c1ccnc(c1)-c1cnc2cc3ccnc3cc12',
         'O=c1ccnc(c1)-c1cc1',
         )
  for smi in smis:
        m = Chem.MolFromSmiles(smi,False)  
        try:
            m.UpdatePropertyCache(False)
            cp = Chem.Mol(m.ToBinary())
            Chem.SanitizeMol(cp)
            m = cp
            print 'fine:',Chem.MolToSmiles(m)
        except ValueError:
            nm=AdjustAromaticNs(m)
            if nm is not None:
                Chem.SanitizeMol(nm)
                print 'fixed:',Chem.MolToSmiles(nm)
            else:
                print 'still broken:',smi

This produces:

.. testoutput::

    fixed: O=c1ccc2ccccc2[nH]1
    fine: Cc1nnnn1C
    fixed: CCc1ccc2[nH]c(=O)c(Cc3nnnn3C3CCCCC3)cc2c1
    fine: C1=Cc2cc3c(cc2=N1)C=CN=3
    fine: C1=Cc2cc3c(cc2=N1)N=CC=3
    fixed: O=c1cc[nH]c(C2=CN=c3cc4c(cc32)=NC=C4)c1
    still broken: O=c1ccnc(c1)-c1cc1

Parallel conformation generation
--------------------------------

Mailing list discussion:
http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02648.html

The code::

  """ contribution from Andrew Dalke """
  import sys
  from rdkit import Chem
  from rdkit.Chem import AllChem

  # Download this from http://pypi.python.org/pypi/futures
  from concurrent import futures

  # Download this from http://pypi.python.org/pypi/progressbar
  import progressbar

  ## On my machine, it takes 39 seconds with 1 worker and 10 seconds with 4.
  ## 29.055u 0.102s 0:28.68 101.6%   0+0k 0+3io 0pf+0w
  #max_workers=1

  ## With 4 threads it takes 11 seconds.
  ## 34.933u 0.188s 0:10.89 322.4%   0+0k 125+1io 0pf+0w
  max_workers=4

  # (The "u"ser time includes time spend in the children processes.
  #  The wall-clock time is 28.68 and 10.89 seconds, respectively.)

  # This function is called in the subprocess.
  # The parameters (molecule and number of conformers) are passed via a Python 
  def generateconformations(m, n):
      m = Chem.AddHs(m)
      ids=AllChem.EmbedMultipleConfs(m, numConfs=n)
      for id in ids:
          AllChem.UFFOptimizeMolecule(m, confId=id)
      # EmbedMultipleConfs returns a Boost-wrapped type which
      # cannot be pickled. Convert it to a Python list, which can.
      return m, list(ids)

  smi_input_file, sdf_output_file = sys.argv[1:3]

  n = int(sys.argv[3])

  writer = Chem.SDWriter(sdf_output_file)

  suppl = Chem.SmilesMolSupplier(smi_input_file, titleLine=False)

  with futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
      # Submit a set of asynchronous jobs
      jobs = []
      for mol in suppl:
          if mol:
              job = executor.submit(generateconformations, mol, n)
              jobs.append(job)

      widgets = ["Generating conformations; ", progressbar.Percentage(), " ", 
                 progressbar.ETA(), " ", progressbar.Bar()]
      pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(jobs))
      for job in pbar(futures.as_completed(jobs)):
          mol,ids=job.result()
          for id in ids:
              writer.write(mol, confId=id)
  writer.close()

Neutralizing Charged Molecules
------------------------------

Mailing list discussion:
http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02648.html

Wiki page: 
http://code.google.com/p/rdkit/wiki/NeutralisingCompounds

The code:

.. testcode::

  """ contribution from Hans de Winter """
  from rdkit import Chem
  from rdkit.Chem import AllChem

  def _InitialiseNeutralisationReactions():
      patts= (
          # Imidazoles
          ('[n+;H]','n'),
          # Amines
          ('[N+;!H0]','N'),
          # Carboxylic acids and alcohols
          ('[$([O-]);!$([O-][#7])]','O'),
          # Thiols
          ('[S-;X1]','S'),
          # Sulfonamides
          ('[$([N-;X2]S(=O)=O)]','N'),
          # Enamines
          ('[$([N-;X2][C,N]=C)]','N'),
          # Tetrazoles
          ('[n-]','[nH]'),
          # Sulfoxides
          ('[$([S-]=O)]','S'),
          # Amides
          ('[$([N-]C=O)]','N'),
          )
      return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

  _reactions=None
  def NeutraliseCharges(smiles, reactions=None):
      global _reactions
      if reactions is None:
          if _reactions is None:
              _reactions=_InitialiseNeutralisationReactions()
          reactions=_reactions
      mol = Chem.MolFromSmiles(smiles)
      replaced = False
      for i,(reactant, product) in enumerate(reactions):
          while mol.HasSubstructMatch(reactant):
              replaced = True
              rms = AllChem.ReplaceSubstructs(mol, reactant, product)
              mol = rms[0]
      if replaced:
          return (Chem.MolToSmiles(mol,True), True)
      else:
          return (smiles, False)

Examples of using it:

.. testcode:: 

  smis=("c1cccc[nH+]1",
        "C[N+](C)(C)C","c1ccccc1[NH3+]",
        "CC(=O)[O-]","c1ccccc1[O-]",
        "CCS",
        "C[N-]S(=O)(=O)C",
        "C[N-]C=C","C[N-]N=C",
        "c1ccc[n-]1",
        "CC[N-]C(=O)CC")
  for smi in smis:
      (molSmiles, neutralised) = NeutraliseCharges(smi)
      print smi,"->",molSmiles

This produces:

.. testoutput::

    c1cccc[nH+]1 -> c1ccncc1
    C[N+](C)(C)C -> C[N+](C)(C)C
    c1ccccc1[NH3+] -> Nc1ccccc1
    CC(=O)[O-] -> CC(=O)O
    c1ccccc1[O-] -> Oc1ccccc1
    CCS -> CCS
    C[N-]S(=O)(=O)C -> CNS(C)(=O)=O
    C[N-]C=C -> C=CNC
    C[N-]N=C -> C=NNC
    c1ccc[n-]1 -> c1cc[nH]c1
    CC[N-]C(=O)CC -> CCNC(=O)CC




License
*******

This document is copyright (C) 2012 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 3.0 License.
To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.


The intent of this license is similar to that of the RDKit itself.
In simple words: “Do whatever you want with it, but please give us some credit.”
