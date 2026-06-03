#
#  Copyright (C) 2003-2021  Greg Landrum and other RDKit contributors
#         All Rights Reserved
#
""" This is a rough coverage test of the python wrapper

it's intended to be shallow, but broad

"""

import doctest
import gc
import gzip
import logging
import os
import sys
import tempfile
import unittest
from io import StringIO

from rdkit import Chem


class TestCase(unittest.TestCase):

  def test_cdxml(self):
    cdxml = """<?xml version="1.0" encoding="UTF-8" ?>
        <!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
        <CDXML
         CreationProgram="ChemDraw JS 2.0.0.9"
         Name="ACS Document 1996"
         BoundingBox="94.75 178.16 154.89 211.17"
         WindowPosition="0 0"
         WindowSize="0 0"
         FractionalWidths="yes"
         InterpretChemically="yes"
         ShowAtomQuery="yes"
         ShowAtomStereo="no"
         ShowAtomEnhancedStereo="yes"
         ShowAtomNumber="no"
         ShowResidueID="no"
         ShowBondQuery="yes"
         ShowBondRxn="yes"
         ShowBondStereo="no"
         ShowTerminalCarbonLabels="no"
         ShowNonTerminalCarbonLabels="no"
         HideImplicitHydrogens="no"
         Magnification="666"
         LabelFont="24"
         LabelSize="10"
         LabelFace="96"
         CaptionFont="24"
         CaptionSize="10"
         HashSpacing="2.50"
         MarginWidth="1.60"
         LineWidth="0.60"
         BoldWidth="2"
         BondLength="14.40"
         BondSpacing="18"
         ChainAngle="120"
         LabelJustification="Auto"
         CaptionJustification="Left"
         AminoAcidTermini="HOH"
         ShowSequenceTermini="yes"
         ShowSequenceBonds="yes"
         ShowSequenceUnlinkedBranches="no"
         ResidueWrapCount="40"
         ResidueBlockCount="10"
         ResidueZigZag="yes"
         NumberResidueBlocks="no"
         PrintMargins="36 36 36 36"
         MacPrintInfo="0003000001200120000000000B6608A0FF84FF880BE309180367052703FC0002000001200120000000000B6608A0000100000064000000010001010100000001270F000100010000000000000000000000000002001901900000000000400000000000000000000100000000000000000000000000000000"
         ChemPropName=""
         ChemPropFormula="Chemical Formula: "
         ChemPropExactMass="Exact Mass: "
         ChemPropMolWt="Molecular Weight: "
         ChemPropMOverZ="m/z: "
         ChemPropAnalysis="Elemental Analysis: "
         ChemPropBoilingPt="Boiling Point: "
         ChemPropMeltingPt="Melting Point: "
         ChemPropCritTemp="Critical Temp: "
         ChemPropCritPres="Critical Pres: "
         ChemPropCritVol="Critical Vol: "
         ChemPropGibbs="Gibbs Energy: "
         ChemPropLogP="Log P: "
         ChemPropMR="MR: "
         ChemPropHenry="Henry&apos;s Law: "
         ChemPropEForm="Heat of Form: "
         ChemProptPSA="tPSA: "
         ChemPropID=""
         ChemPropFragmentLabel=""
         color="0"
         bgcolor="1"
         RxnAutonumberStart="1"
         RxnAutonumberConditions="no"
         RxnAutonumberStyle="Roman"
        ><colortable>
        <color r="1" g="1" b="1"/>
        <color r="0" g="0" b="0"/>
        <color r="1" g="0" b="0"/>
        <color r="1" g="1" b="0"/>
        <color r="0" g="1" b="0"/>
        <color r="0" g="1" b="1"/>
        <color r="0" g="0" b="1"/>
        <color r="1" g="0" b="1"/>
        </colortable><fonttable>
        <font id="24" charset="utf-8" name="Arial"/>
        </fonttable><page
         id="32"
         BoundingBox="0 0 542 354"
         Width="542"
         Height="354"
         HeaderPosition="36"
         FooterPosition="36"
         PageOverlap="0"
         PrintTrimMarks="yes"
         HeightPages="1"
         WidthPages="1"
         DrawingSpace="poster"
        ><fragment
         id="10"
         BoundingBox="94.75 178.16 154.89 211.17"
         Z="4"
        ><n
         id="7"
         p="95.05 187.47"
         Z="1"
         AS="N"
        /><n
         id="9"
         p="95.05 201.87"
         Z="3"
         AS="N"
        /><n
         id="11"
         p="106.31 210.84"
         Z="5"
         AS="N"
        /><n
         id="13"
         p="120.35 207.64"
         Z="7"
         AS="N"
        /><n
         id="15"
         p="126.59 194.67"
         Z="9"
         AS="N"
        /><n
         id="17"
         p="120.35 181.69"
         Z="11"
         AS="N"
        /><n
         id="19"
         p="106.31 178.49"
         Z="13"
         AS="N"
        /><n
         id="28"
         p="140.99 194.67"
         Z="22"
         NodeType="Nickname"
         NeedsClean="yes"
         AS="N"
        ><fragment
         id="33"
        ><n
         id="34"
         p="148.17 207.09"
         Element="8"
         NumHydrogens="0"
        /><n
         id="35"
         p="162.52 207.09"
        /><n
         id="36"
         p="176.87 207.09"
        /><n
         id="37"
         p="169.69 194.67"
        /><n
         id="38"
         p="169.69 219.52"
        /><n
         id="39"
         p="140.99 194.67"
        /><n
         id="40"
         p="148.17 182.24"
         Element="8"
         NumHydrogens="0"
        /><n
         id="41"
         p="126.64 194.67"
         NodeType="ExternalConnectionPoint"
        /><b
         id="42"
         B="39"
         E="40"
         Order="2"
        /><b
         id="43"
         B="35"
         E="38"
        /><b
         id="44"
         B="35"
         E="36"
        /><b
         id="45"
         B="35"
         E="37"
        /><b
         id="46"
         B="34"
         E="35"
        /><b
         id="47"
         B="34"
         E="39"
        /><b
         id="48"
         B="41"
         E="39"
        /></fragment><t
         p="137.66 198.28"
         BoundingBox="137.66 189.64 154.89 198.28"
         LabelJustification="Left"
         LabelAlignment="Left"
        ><s font="24" size="9.95" color="0" face="96">Boc</s></t></n><b
         id="21"
         Z="15"
         B="7"
         E="9"
         BS="N"
        /><b
         id="22"
         Z="16"
         B="9"
         E="11"
         BS="N"
        /><b
         id="23"
         Z="17"
         B="11"
         E="13"
         BS="N"
        /><b
         id="24"
         Z="18"
         B="13"
         E="15"
         BS="N"
        /><b
         id="25"
         Z="19"
         B="15"
         E="17"
         BS="N"
        /><b
         id="26"
         Z="20"
         B="17"
         E="19"
         BS="N"
        /><b
         id="27"
         Z="21"
         B="19"
         E="7"
         BS="N"
        /><b
         id="29"
         Z="23"
         B="15"
         E="28"
         BS="N"
        /></fragment></page></CDXML>"""
    mols = Chem.MolsFromCDXML(cdxml)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmiles(mols[0]), "CC(C)(C)OC(=O)C1CCCCCC1")

    mols = Chem.MolsFromCDXML(cdxml, True, False)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmiles(mols[0]), "CC(C)(C)OC(=O)C1CCCCCC1")

    mols = Chem.MolsFromCDXML(cdxml, False, False)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmiles(mols[0]), "CC(C)(C)OC(=O)C1CCCCCC1")

    params = Chem.CDXMLParserParams()
    mols = Chem.MolsFromCDXML(cdxml, params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmiles(mols[0]), "CC(C)(C)OC(=O)C1CCCCCC1")

    params.sanitize = True
    params.removeHs = False
    mols = Chem.MolsFromCDXML(cdxml, params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmiles(mols[0]), "CC(C)(C)OC(=O)C1CCCCCC1")

    params.sanitize = False
    params.removeHs = False

    mols = Chem.MolsFromCDXML(cdxml, params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmiles(mols[0]), "CC(C)(C)OC(=O)C1CCCCCC1")
    
    
    
  def test_cdxml(self):
    try: from rdkit.Chem import rdChemDraw
    except:
      return

    rdbase = os.environ['RDBASE']
    cdxfilename = os.path.join(rdbase,
                            'Code/GraphMol/test_data/CDXML/ring-stereo1.cdx')
    mols = Chem.MolsFromCDXMLFile(cdxfilename)
    filename = os.path.join(rdbase,
                            'Code/GraphMol/test_data/CDXML/ring-stereo1.cdxml')
    mols2 = Chem.MolsFromCDXMLFile(filename)
    smi1 = [Chem.MolToSmiles(m) for m in mols]
    smi2 = [Chem.MolToSmiles(m) for m in mols2]
    self.assertEqual(smi1, smi2)

    self.assertEqual(smi1, ['C1CC[C@H]2CCCC[C@H]2C1'])
    with open(cdxfilename, 'rb') as f:
      data = f.read()
    params = Chem.CDXMLParserParams(True, True, Chem.CDXMLFormat.CDX)
    mols3 = Chem.MolsFromCDXML(data, params)
    smi3 = [Chem.MolToSmiles(m) for m in mols3]
    self.assertEqual(smi1, smi3)
    if Chem.HasChemDrawCDXSupport():
      # ensure we can round trip through CDXML, CDX
      for smi, mol in zip(smi3, mols3):
        cdxml = Chem.MolToCDXMLBlock(mol)
        cdxml2 = Chem.MolToCDXMLBlock(mol, Chem.CDXMLFormat.CDXML)
        # check default is cdxml
        self.assertEqual(cdxml, cdxml2)
        cdx = Chem.MolToCDXMLBlock(mol, Chem.CDXMLFormat.CDX)
        self.assertEqual(type(cdx), bytes)

        self.assertEqual(Chem.MolToSmiles(Chem.MolsFromCDXML(cdxml)[0]),
                         smi)
        self.assertEqual(Chem.MolToSmiles(Chem.MolsFromCDXML(cdx, params)[0]),
                         smi)
        
    

  def test_formats(self):
    try:
      from rdkit.Chem import rdChemDraw
      self.assertEqual(Chem.HasChemDrawCDXSupport(),True)
    except:
      self.assertEqual(Chem.HasChemDrawCDXSupport(),False)
      return

    rdbase = os.environ['RDBASE']
    cdxfilename = os.path.join(rdbase,
                            'Code/GraphMol/test_data/CDXML/ring-stereo1.cdx')
    mols = Chem.MolsFromCDXMLFile(cdxfilename)
    cdxmlfilename = os.path.join(rdbase,
                            'Code/GraphMol/test_data/CDXML/ring-stereo1.cdxml')

    tests = [
      # we can deduce extensions from filenames, but not from streams (yet!)
      # filename, Stream, IsCDX, CDX res, CDXML res, Auto Res
      (cdxfilename, True, True, True, False, False),
      (cdxfilename, False, True, True, False, True),
      (cdxmlfilename, True, False, False, True, True),
      (cdxmlfilename, False, False, False, True, True),
      ]

    for filename, stream, iscdx, cdxres, cdxmlres, autores in tests:
        for format, res in zip([Chem.CDXMLFormat.CDX, Chem.CDXMLFormat.CDXML, Chem.CDXMLFormat.Auto],
                               [cdxres, cdxmlres, autores]):
          if stream:
            with open(filename, 'rb') as f:
              data = f.read()
              try:
                mols = Chem.MolsFromCDXML(data, Chem.CDXMLParserParams(True, True, format))
                if res: assert mols
              except RuntimeError:
                assert res == False
          else:
            mols = Chem.MolsFromCDXMLFile(filename, Chem.CDXMLParserParams(True, True, format))
            if res: assert mols

  def test_cdxml_query_api(self):
    rdbase = os.environ['RDBASE']
    filename = os.path.join(rdbase, 'rdkit/Chem/test_data/benzene.cdxml')

    params = Chem.CDXMLParserParams()
    self.assertFalse(params.parseQueries)
    self.assertFalse(params.strictQueryParsing)

    params.parseQueries = True
    params.strictQueryParsing = True

    file_mols = Chem.MolsFromCDXMLFile(filename, params)
    self.assertEqual(len(file_mols), 1)
    self.assertEqual(Chem.MolToSmiles(file_mols[0]), 'c1ccccc1')

    with open(filename) as inf:
      block = inf.read()

    block_mols = Chem.MolsFromCDXML(block, params)
    self.assertEqual(len(block_mols), 1)
    self.assertEqual(Chem.MolToSmiles(block_mols[0]), 'c1ccccc1')

  def test_cdxml_bond_queries(self):
    rdbase = os.environ['RDBASE']
    query_base = os.path.join(rdbase, 'Code/GraphMol/test_data/CDXML/queries')
    cases = [
      ('anybond', os.path.join(rdbase, 'Code/GraphMol/test_data/CDXML/anybond.cdxml'),
       'file', '[#6]1~[#6]-[#6]-[#6]-[#6]-[#6]-1'),
      ('single-or-double', os.path.join(query_base, 'furan_sd.cdxml'),
       'file', '[#8]1:[#6]-,=[#6]:[#6]-,=[#6]:1'),
      ('single-or-aromatic', os.path.join(query_base, 'furan_sa.cdxml'),
       'file', '[#8]1:[#6][#6]:[#6][#6]:1'),
      ('double-or-aromatic', os.path.join(query_base, 'furan_da.cdxml'),
       'file', '[#8]1:[#6]=,:[#6]:[#6]=,:[#6]:1'),
      ('ring-topology', os.path.join(query_base, 'CCOC_Rng.cdxml'),
       'file', '[#8](-[#6]-&@[#6])-[#6]'),
      ('chain-topology', os.path.join(query_base, 'CCOC_Chn.cdxml'),
       'file', '[#8](-[#6]-&!@[#6])-[#6]'),
    ]
    params = Chem.CDXMLParserParams()
    params.parseQueries = True

    for _, filename, mode, expected_smarts in cases:
      if mode == 'file':
        mols = Chem.MolsFromCDXMLFile(filename, params)
      else:
        with open(filename) as inf:
          mols = Chem.MolsFromCDXML(inf.read(), params)

      self.assertEqual(len(mols), 1)
      self.assertEqual(Chem.MolToSmarts(mols[0]), expected_smarts)

  def test_cdxml_hydrogen_bond_query(self):
    rdbase = os.environ['RDBASE']
    query_base = os.path.join(rdbase, 'Code/GraphMol/test_data/CDXML/queries')
    params = Chem.CDXMLParserParams()
    params.parseQueries = True
    mols = Chem.MolsFromCDXMLFile(
      os.path.join(query_base, 'qbond_hydrogen.cdxml'),
      params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmarts(mols[0]), '[#8][#8]')
    bond = mols[0].GetBondWithIdx(0)
    self.assertEqual(bond.GetBondType(), Chem.BondType.HYDROGEN)
    self.assertIn('H:0.0', Chem.MolToCXSmarts(mols[0]))

  def test_cdxml_multiattachment_query(self):
    rdbase = os.environ['RDBASE']
    params = Chem.CDXMLParserParams()
    params.sanitize = False
    params.parseQueries = True
    mols = Chem.MolsFromCDXMLFile(
      os.path.join(rdbase, 'rdkit/Chem/test_data/ferrocene.cdxml'),
      params)
    self.assertEqual(len(mols), 1)

    dummy_atoms = [atom for atom in mols[0].GetAtoms() if atom.GetAtomicNum() == 0]
    self.assertEqual(len(dummy_atoms), 2)

    endpoint_sets = []
    for bond in mols[0].GetBonds():
      if not bond.HasProp('_MolFileBondEndPts'):
        continue
      self.assertTrue(bond.HasProp('_MolFileBondAttach'))
      self.assertEqual(bond.GetProp('_MolFileBondAttach'), 'ANY')
      endpoint_sets.append(bond.GetProp('_MolFileBondEndPts'))

    self.assertEqual(sorted(endpoint_sets), ['(5 1 2 3 4 5)', '(5 6 7 8 9 10)'])

  def test_cdxml_external_connection_fragment_query(self):
    rdbase = os.environ['RDBASE']
    params = Chem.CDXMLParserParams()
    params.parseQueries = True
    mols = Chem.MolsFromCDXMLFile(
      os.path.join(rdbase, 'External/ChemDraw/test_data/atom-to-fragment.cdxml'),
      params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmarts(mols[0]), '[#6]-[#6]=[#6]=[#6](-[#6])-[#6]')

  def test_cdxml_spiro_ring_bond_count_query(self):
    rdbase = os.environ['RDBASE']
    query_base = os.path.join(rdbase, 'Code/GraphMol/test_data/CDXML/queries')
    params = Chem.CDXMLParserParams()
    params.parseQueries = True
    mols = Chem.MolsFromCDXMLFile(
      os.path.join(query_base, 'qrestrict_ringbond_spiro.cdxml'),
      params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmarts(mols[0]), '[#6]1-[#6]-[#7]-[#6&x{4-}]-[#6]-[#6]-1')

    self.assertFalse(Chem.MolFromSmiles('N1CCCCC1').HasSubstructMatch(mols[0]))
    self.assertFalse(Chem.MolFromSmiles('N1CCC2(CC1)CCCC2').HasSubstructMatch(mols[0]))
    self.assertTrue(Chem.MolFromSmiles('N1C2(CCCCC2)CCCC1').HasSubstructMatch(mols[0]))
    self.assertFalse(Chem.MolFromSmiles('N1CCC2CCCCC2C1').HasSubstructMatch(mols[0]))

  def test_cdxml_explicit_hydrogen_minimum_queries(self):
    rdbase = os.environ['RDBASE']
    params = Chem.CDXMLParserParams()
    params.parseQueries = True
    params.strictQueryParsing = True

    mols = Chem.MolsFromCDXMLFile(
      os.path.join(rdbase, 'Code/GraphMol/test_data/CDXML/query-atoms.cdxml'), params)
    self.assertEqual(len(mols), 3)
    self.assertIn('!H0', Chem.MolToSmarts(mols[0]))
    self.assertTrue(Chem.MolFromSmiles('Cc1ccccc1').HasSubstructMatch(mols[0]))
    self.assertFalse(Chem.MolFromSmiles('Cc1ccc(C)cc1').HasSubstructMatch(mols[0]))

    mols = Chem.MolsFromCDXMLFile(
      os.path.join(rdbase, 'Code/GraphMol/test_data/CDXML/chirality1.cdxml'), params)
    self.assertEqual(len(mols), 1)
    smarts = Chem.MolToSmarts(mols[0])
    self.assertNotIn('!H0', smarts)
    self.assertNotIn('!H1', smarts)

  def test_cdxml_atom_restriction_queries(self):
    rdbase = os.environ['RDBASE']
    query_base = os.path.join(rdbase, 'Code/GraphMol/test_data/CDXML/queries')
    smarts_cases = [
      ('free-sites', os.path.join(query_base, 'qrestrict_freesites_1.cdxml'),
       '[#6]-[#6]-[!#1&D{1-2}]'),
      ('implicit-hydrogens', os.path.join(query_base, 'qrestrict_implicit_hs.cdxml'),
       '[#6]-[#6]-[!#1&h0]'),
      ('ring-bond-count-as-drawn', os.path.join(query_base, 'qrestrict_ringbond_asdrawn.cdxml'),
       '[#6]-[#6]-[!#1&x0]'),
      ('ring-bond-count', os.path.join(query_base, 'qrestrict_ringbond_simple.cdxml'),
       '[#6]-[#6]-[!#1&x2]'),
      ('not-list', os.path.join(query_base, 'qatom_notlist.cdxml'),
       '[#6]-[#6]-[!#6&!#7&!#8]'),
      ('substituents-exactly', os.path.join(query_base, 'qrestrict_sub_exact_2.cdxml'),
       '[#6]-[#6]-[!#1&D2]'),
      ('substituents-up-to', os.path.join(query_base, 'qrestrict_sub_upto_2.cdxml'),
       '[#6]-[#6]-[!#1&D{0-2}]'),
      ('unsaturated-bonds', os.path.join(query_base, 'qrestrict_unsat_present.cdxml'),
       '[#6]-[#6]-[!#1&$(*=,:,#*)]'),
    ]
    prop_cases = [
      ('link-node', os.path.join(query_base, 'qlinknode_1_3.cdxml'),
       '_molLinkNodes', '1 3 2 2 1 2 3'),
      ('reaction-stereo', os.path.join(query_base, 'qrestrict_rxnstereo_inversion.cdxml'),
       'molInversionFlag', 1),
    ]
    params = Chem.CDXMLParserParams()
    params.parseQueries = True

    for _, filename, expected_smarts in smarts_cases:
      mols = Chem.MolsFromCDXMLFile(filename, params)
      self.assertEqual(len(mols), 1)
      self.assertEqual(Chem.MolToSmarts(mols[0]), expected_smarts)

    for _, filename, prop_name, expected_value in prop_cases:
      mols = Chem.MolsFromCDXMLFile(filename, params)
      self.assertEqual(len(mols), 1)
      if prop_name == '_molLinkNodes':
        self.assertTrue(mols[0].HasProp(prop_name))
        self.assertEqual(mols[0].GetProp(prop_name), expected_value)
      else:
        atom = mols[0].GetAtomWithIdx(2)
        self.assertTrue(atom.HasProp(prop_name))
        self.assertEqual(atom.GetIntProp(prop_name), expected_value)

    mols = Chem.MolsFromCDXMLFile(
      os.path.join(query_base, 'qvarattach.cdxml'),
      params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(mols[0].GetAtomWithIdx(3).GetAtomicNum(), 0)
    bond = mols[0].GetBondBetweenAtoms(1, 3)
    self.assertIsNotNone(bond)
    self.assertTrue(bond.HasProp('_MolFileBondAttach'))
    self.assertEqual(bond.GetProp('_MolFileBondAttach'), 'ANY')
    self.assertTrue(bond.HasProp('_MolFileBondEndPts'))
    self.assertEqual(bond.GetProp('_MolFileBondEndPts'), '(2 1 3)')

    mols = Chem.MolsFromCDXMLFile(
      os.path.join(query_base, 'qecp_rgroup.cdxml'),
      params)
    self.assertEqual(len(mols), 1)
    self.assertEqual(Chem.MolToSmarts(mols[0]), '[#6]-[*:1]')
    atom = mols[0].GetAtomWithIdx(1)
    self.assertEqual(atom.GetAtomMapNum(), 1)
    self.assertTrue(atom.HasProp('atomLabel'))
    self.assertEqual(atom.GetProp('atomLabel'), 'R')
          

        


if __name__ == '__main__':
  if "RDTESTCASE" in os.environ:
    suite = unittest.TestSuite()
    testcases = os.environ["RDTESTCASE"]
    for name in testcases.split(':'):
      suite.addTest(TestCase(name))

    runner = unittest.TextTestRunner()
    runner.run(suite)
  else:
    unittest.main()
