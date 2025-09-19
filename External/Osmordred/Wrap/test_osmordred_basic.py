"""
Test module for RDKit Osmordred functionality.

This module provides comprehensive tests for all Osmordred descriptor functions
available through the Python wrapper.
"""

import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdOsmordred


def to_list(result):
    """Convert RDKit vector objects to Python lists."""
    if hasattr(result, '__iter__') and not isinstance(result, (list, tuple)):
        return list(result)
    return result


class TestOsmordred(unittest.TestCase):
    """Test class for Osmordred functionality."""

    def setUp(self):
        """Set up test molecules."""
        self.ethanol = Chem.MolFromSmiles("CCO")
        self.benzene = Chem.MolFromSmiles("c1ccccc1")
        self.complex_mol = Chem.MolFromSmiles(
            "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
        )
        self.t_butyl_benzene = Chem.MolFromSmiles("CC(C)(C)c1ccc(C(C)(C)C)cc1")

    def test_basic_functionality(self):
        """Test basic Osmordred functions with simple molecules."""
        self.assertIsNotNone(self.ethanol)
        
        # Test atom counts (using version 2 as default)
        atom_counts = rdOsmordred.CalcAtomCount(self.ethanol)
        atom_counts_list = to_list(atom_counts)
        self.assertGreater(len(atom_counts_list), 0)
        
        # Test bond counts
        bond_counts = rdOsmordred.CalcBondCount(self.ethanol)
        bond_counts_list = to_list(bond_counts)
        self.assertGreater(len(bond_counts_list), 0)
        
        # Test molecular weight
        weight = rdOsmordred.CalcWeight(self.ethanol)
        weight_list = to_list(weight)
        self.assertGreater(len(weight_list), 0)
        self.assertGreater(weight_list[0], 0.0)

    def test_aromatic_functionality(self):
        """Test aromatic-related functions."""
        self.assertIsNotNone(self.benzene)
        
        # Test aromatic counts
        aromatic_counts = rdOsmordred.CalcAromatic(self.benzene)
        aromatic_counts_list = to_list(aromatic_counts)
        self.assertGreater(len(aromatic_counts_list), 0)
        
        # Test ring count
        ring_counts = rdOsmordred.CalcRingCount(self.benzene)
        ring_counts_list = to_list(ring_counts)
        self.assertGreater(len(ring_counts_list), 0)
        
        # Test aromatic atom count
        aromatic_atoms = rdOsmordred.CalcCountAromaticAtoms(self.benzene)
        self.assertIsInstance(aromatic_atoms, int)
        self.assertGreaterEqual(aromatic_atoms, 0)
        
        # Test aromatic bond count
        aromatic_bonds = rdOsmordred.CalcCountAromaticBonds(self.benzene)
        self.assertIsInstance(aromatic_bonds, int)
        self.assertGreaterEqual(aromatic_bonds, 0)

    def test_matrix_descriptors(self):
        """Test matrix-based descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test distance matrix descriptors
        dist_descs = rdOsmordred.CalcDistanceMatrix(self.ethanol)
        dist_descs_list = to_list(dist_descs)
        self.assertGreater(len(dist_descs_list), 0)
        
        dist_descs_eigen = rdOsmordred.CalcDistanceMatrixEigen(self.ethanol)
        dist_descs_eigen_list = to_list(dist_descs_eigen)
        self.assertGreater(len(dist_descs_eigen_list), 0)
        
        # Test adjacency matrix descriptors
        adj_descs = rdOsmordred.CalcAdjacencyMatrix(self.ethanol)
        adj_descs_list = to_list(adj_descs)
        self.assertGreater(len(adj_descs_list), 0)
        
        adj_descs_eigen = rdOsmordred.CalcAdjacencyMatrixEigen(self.ethanol)
        adj_descs_eigen_list = to_list(adj_descs_eigen)
        self.assertGreater(len(adj_descs_eigen_list), 0)
        
        # Test detour matrix descriptors
        detour_descs = rdOsmordred.CalcDetourMatrix(self.ethanol)
        detour_descs_list = to_list(detour_descs)
        self.assertGreater(len(detour_descs_list), 0)
        
        detour_descs_eigen = rdOsmordred.CalcDetourMatrixEigen(self.ethanol)
        detour_descs_eigen_list = to_list(detour_descs_eigen)
        self.assertGreater(len(detour_descs_eigen_list), 0)
        
        # Test Barysz matrix descriptors
        barysz_descs = rdOsmordred.CalcBaryszMatrix(self.ethanol)
        barysz_descs_list = to_list(barysz_descs)
        self.assertGreater(len(barysz_descs_list), 0)
        
        barysz_descs_eigen = rdOsmordred.CalcBaryszMatrixEigen(self.ethanol)
        barysz_descs_eigen_list = to_list(barysz_descs_eigen)
        self.assertGreater(len(barysz_descs_eigen_list), 0)

    def test_physical_properties(self):
        """Test physical property calculations."""
        self.assertIsNotNone(self.ethanol)
        
        # Test molecular weight
        weight = rdOsmordred.CalcWeight(self.ethanol)
        weight_list = to_list(weight)
        self.assertGreater(len(weight_list), 0)
        self.assertGreater(weight_list[0], 0.0)
        
        # Test van der Waals volume
        vdw_vol = rdOsmordred.CalcVdwVolumeABC(self.ethanol)
        vdw_vol_list = to_list(vdw_vol)
        self.assertGreater(vdw_vol_list[0], 0.0)
        
        # Test McGowan volume
        mcgowan_vol = rdOsmordred.CalcMcGowanVolume(self.ethanol)
        mcgowan_vol_list = to_list(mcgowan_vol)
        self.assertGreater(mcgowan_vol_list[0], 0.0)
        
        # Test polarizability
        polarizability = rdOsmordred.CalcPolarizability(self.ethanol)
        polarizability_list = to_list(polarizability)
        self.assertGreater(polarizability_list[0], 0.0)
        
        # Test SLogP
        slogp = rdOsmordred.CalcSLogP(self.ethanol)
        slogp_list = to_list(slogp)
        self.assertIsInstance(slogp_list, list)
        
        # Test LogS
        logs = rdOsmordred.CalcLogS(self.ethanol)
        logs_list = to_list(logs)
        self.assertIsInstance(logs_list, list)

    def test_topological_descriptors(self):
        """Test topological descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test Wiener index
        wiener = rdOsmordred.CalcWienerIndex(self.ethanol)
        wiener_list = to_list(wiener)
        self.assertGreater(len(wiener_list), 0)
        
        # Test Balaban J index
        balaban_j = rdOsmordred.CalcBalabanJ(self.ethanol)
        balaban_j_list = to_list(balaban_j)
        self.assertGreater(len(balaban_j_list), 0)
        
        # Test Bertz CT index
        bertz_ct = rdOsmordred.CalcBertzCT(self.ethanol)
        bertz_ct_list = to_list(bertz_ct)
        self.assertGreater(len(bertz_ct_list), 0)
        
        # Test Zagreb index
        zagreb = rdOsmordred.CalcZagrebIndex(self.ethanol)
        zagreb_list = to_list(zagreb)
        self.assertGreater(len(zagreb_list), 0)
        
        # Test eccentric connectivity index
        eccentric = rdOsmordred.CalcEccentricConnectivityIndex(self.ethanol)
        eccentric_list = to_list(eccentric)
        self.assertGreater(len(eccentric_list), 0)

    def test_chi_descriptors(self):
        """Test Chi descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test all Chi descriptors
        chi_descs = rdOsmordred.CalcChi(self.ethanol)
        chi_descs_list = to_list(chi_descs)
        self.assertGreater(len(chi_descs_list), 0)
        
        # Test individual Chi descriptors
        chipath = rdOsmordred.CalcChipath(self.ethanol)
        chipath_list = to_list(chipath)
        self.assertGreater(len(chipath_list), 0)
        
        chichain = rdOsmordred.CalcChichain(self.ethanol)
        chichain_list = to_list(chichain)
        self.assertGreater(len(chichain_list), 0)
        
        chicluster = rdOsmordred.CalcChicluster(self.ethanol)
        chicluster_list = to_list(chicluster)
        self.assertGreater(len(chicluster_list), 0)
        
        chipathcluster = rdOsmordred.CalcChipathcluster(self.ethanol)
        chipathcluster_list = to_list(chipathcluster)
        self.assertGreater(len(chipathcluster_list), 0)

    def test_estate_descriptors(self):
        """Test EState descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test EState descriptors
        estate_descs = rdOsmordred.CalcEState(self.ethanol, False)
        estate_descs_list = to_list(estate_descs)
        self.assertGreater(len(estate_descs_list), 0)
        
        # Test extended EState descriptors
        estate_descs_ext = rdOsmordred.CalcEState(self.ethanol, True)
        estate_descs_ext_list = to_list(estate_descs_ext)
        self.assertGreater(len(estate_descs_ext_list), 0)
        
        # Test BEState descriptors
        bestate_descs = rdOsmordred.CalcBEState(self.ethanol)
        bestate_descs_list = to_list(bestate_descs)
        self.assertGreater(len(bestate_descs_list), 0)
        
        # Test HEState descriptors
        hestate_descs = rdOsmordred.CalcHEState(self.ethanol)
        hestate_descs_list = to_list(hestate_descs)
        self.assertGreater(len(hestate_descs_list), 0)

    def test_count_functions(self):
        """Test count-based functions."""
        self.assertIsNotNone(self.ethanol)
        
        # Test acidic group count
        acidic_count = rdOsmordred.CalcAcidicGroupCount(self.ethanol)
        self.assertIsInstance(acidic_count, int)
        self.assertGreaterEqual(acidic_count, 0)
        
        # Test basic group count
        basic_count = rdOsmordred.CalcBasicGroupCount(self.ethanol)
        self.assertIsInstance(basic_count, int)
        self.assertGreaterEqual(basic_count, 0)
        
        # Test aromatic atom count
        aromatic_atoms = rdOsmordred.CalcCountAromaticAtoms(self.ethanol)
        self.assertIsInstance(aromatic_atoms, int)
        self.assertGreaterEqual(aromatic_atoms, 0)
        
        # Test aromatic bond count
        aromatic_bonds = rdOsmordred.CalcCountAromaticBonds(self.ethanol)
        self.assertIsInstance(aromatic_bonds, int)
        self.assertGreaterEqual(aromatic_bonds, 0)

    def test_complex_molecule(self):
        """Test functions with complex molecules."""
        self.assertIsNotNone(self.complex_mol)
        
        # Test that all functions work with complex molecules
        # Functions that need version parameter
        functions_with_version = [
        ]
        
        # Functions that need no additional parameters
        functions_simple = [
            rdOsmordred.CalcABCIndex,
            rdOsmordred.CalcAcidBase,
            rdOsmordred.CalcAromatic,
            rdOsmordred.CalcBalabanJ,
            rdOsmordred.CalcBertzCT,
            rdOsmordred.CalcBondCount,
            rdOsmordred.CalcVertexAdjacencyInformation,
            rdOsmordred.CalcWeight,
            rdOsmordred.CalcWienerIndex,
            rdOsmordred.CalcVdwVolumeABC,
            rdOsmordred.CalcTopoPSA,
            rdOsmordred.CalcSLogP,
            rdOsmordred.CalcHydrogenBond,
            rdOsmordred.CalcLogS,
            rdOsmordred.CalcLipinski,
            rdOsmordred.CalcMcGowanVolume,
            rdOsmordred.CalcPolarizability,
            rdOsmordred.CalcRotatableBond,
            rdOsmordred.CalcFragmentComplexity,
            rdOsmordred.CalcConstitutional,
            rdOsmordred.CalcTopologicalIndex,
        ]
        
        # Test functions with version parameter
        for func, version in functions_with_version:
            try:
                result = func(self.complex_mol, version)
                result_list = to_list(result)
                self.assertGreater(len(result_list), 0)
            except Exception as e:
                self.fail(f"Function {func.__name__} failed: {e}")
        
        # Test functions without additional parameters
        for func in functions_simple:
            try:
                result = func(self.complex_mol)
                result_list = to_list(result)
                self.assertGreater(len(result_list), 0)
            except Exception as e:
                self.fail(f"Function {func.__name__} failed: {e}")

    def test_edge_cases(self):
        """Test edge cases and error handling."""
        # Test with None molecule
        with self.assertRaises((TypeError, AttributeError)):
            rdOsmordred.CalcAtomCount(None)
        
        # Test with empty SMILES
        empty_mol = Chem.MolFromSmiles("")
        if empty_mol is not None:
            atom_counts = rdOsmordred.CalcAtomCount(empty_mol)
            atom_counts_list = to_list(atom_counts)
            self.assertIsInstance(atom_counts_list, list)
        
        # Test with single atom
        single_atom = Chem.MolFromSmiles("[He]")
        if single_atom is not None:
            atom_counts = rdOsmordred.CalcAtomCount(single_atom)
            atom_counts_list = to_list(atom_counts)
            self.assertIsInstance(atom_counts_list, list)
            self.assertGreater(len(atom_counts_list), 0)

    def test_consistency(self):
        """Test that multiple calls return consistent results."""
        self.assertIsNotNone(self.ethanol)
        
        # Test atom counts consistency
        atom_counts1 = rdOsmordred.CalcAtomCount(self.ethanol)
        atom_counts2 = rdOsmordred.CalcAtomCount(self.ethanol)
        atom_counts1_list = to_list(atom_counts1)
        atom_counts2_list = to_list(atom_counts2)
        self.assertEqual(atom_counts1_list, atom_counts2_list)
        
        # Test bond counts consistency
        bond_counts1 = rdOsmordred.CalcBondCount(self.ethanol)
        bond_counts2 = rdOsmordred.CalcBondCount(self.ethanol)
        bond_counts1_list = to_list(bond_counts1)
        bond_counts2_list = to_list(bond_counts2)
        self.assertEqual(bond_counts1_list, bond_counts2_list)
        
        # Test molecular weight consistency
        weight1 = rdOsmordred.CalcWeight(self.ethanol)
        weight2 = rdOsmordred.CalcWeight(self.ethanol)
        weight1_list = to_list(weight1)
        weight2_list = to_list(weight2)
        self.assertEqual(weight1_list, weight2_list)

    def test_information_content(self):
        """Test information content descriptors."""
        self.assertIsNotNone(self.ethanol)
        
        # Test information content
        info_content = rdOsmordred.CalcInformationContent(self.ethanol, 5)
        info_content_list = to_list(info_content)
        self.assertGreater(len(info_content_list), 0)
        
        # Test with different max radius
        info_content_r3 = rdOsmordred.CalcInformationContent(self.ethanol, 3)
        info_content_r3_list = to_list(info_content_r3)
        self.assertGreater(len(info_content_r3_list), 0)
        
        info_content_r7 = rdOsmordred.CalcInformationContent(self.ethanol, 7)
        info_content_r7_list = to_list(info_content_r7)
        self.assertGreater(len(info_content_r7_list), 0)

    def test_all_functions_exist(self):
        """Test that all expected functions are available."""
        expected_functions = [
            'CalcABCIndex', 'CalcAcidBase', 'CalcAromatic', 'CalcAtomCount',
            'CalcBalabanJ', 'CalcBertzCT', 'CalcBondCount', 'CalcVertexAdjacencyInformation',
            'CalcWeight', 'CalcWienerIndex', 'CalcVdwVolumeABC', 'CalcTopoPSA',
            'CalcSLogP', 'CalcHydrogenBond', 'CalcLogS', 'CalcLipinski',
            'CalcMcGowanVolume', 'CalcPolarizability', 'CalcRotatableBond',
            'CalcFragmentComplexity', 'CalcConstitutional', 'CalcTopologicalIndex',
            'CalcDetourMatrixEigen', 'CalcDetourMatrix', 'CalcDistanceMatrixEigen',
            'CalcDistanceMatrix', 'CalcAdjacencyMatrixEigen', 'CalcAdjacencyMatrix',
            'CalcCarbonTypes', 'CalcEccentricConnectivityIndex', 'CalcBaryszMatrix',
            'CalcBaryszMatrixEigen', 'CalcZagrebIndex', 'CalcMoeType',
            'CalcMolecularDistanceEdge', 'CalcEState', 'CalcWalkCount',
            'CalcTopologicalCharge', 'CalcChi', 'CalcPathCount', 'CalcKappaShapeIndex',
            'CalcRingCount', 'CalcMolecularId', 'CalcBCUT', 'CalcAutocorrelation',
            'CalcFramework', 'CalcExtendedTopochemicalAtom',
            'CalcChipath', 'CalcChichain', 'CalcChicluster', 'CalcChipathcluster',
            'CalcAcidicGroupCount', 'CalcBasicGroupCount', 'CalcCountAromaticAtoms',
            'CalcCountAromaticBonds', 'CalcBEState', 'CalcHEState',
            'CalcAlphaKappaShapeIndex', 'CalcAbrahams', 'CalcPol', 'CalcMR',
            'CalcFlexibility', 'CalcODT', 'CalcSchultz', 'CalcRNCGRPCG',
            'CalcAZV', 'CalcASV', 'CalcDSV', 'CalcAZS', 'CalcASZ', 'CalcDN2S',
            'CalcDN2I', 'CalcASI', 'CalcDSI', 'CalcASN', 'CalcDSN', 'CalcDN2N',
            'CalcANS', 'CalcANV', 'CalcAZN', 'CalcANZ', 'CalcANI', 'CalcDSZ',
            'CalcANN', 'CalcDN2Z', 'CalcANMat', 'CalcAZMat', 'CalcASMat',
            'CalcDSMat', 'CalcDN2Mat', 'CalcFrags', 'CalcAddFeatures',
            'CalcInformationContent'
        ]
        
        for func_name in expected_functions:
            self.assertTrue(hasattr(rdOsmordred, func_name),
                          f"Function {func_name} not found in Osmordred module")


if __name__ == '__main__':
    unittest.main() 
