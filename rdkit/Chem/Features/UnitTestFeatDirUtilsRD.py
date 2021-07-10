import unittest

from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.Chem.Features.FeatDirUtilsRD import GetDonor2FeatVects

class TestCase(unittest.TestCase):
    def assertListAlmostEqual(self, list1, list2, msg, tol=7):
        self.assertEqual(len(list1), len(list2), msg)
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, tol, msg)

    def setUp(self):
        #Define molecule for using in tests of GetDonor2FeatVects
        self.mol = Chem.MolFromSmiles('C=CCONC')
        emol = Chem.RWMol(self.mol)
        emol = Chem.AddHs(emol)
        emol.AddConformer(Chem.Conformer(15))
        emol.GetConformer().SetAtomPosition(0, [-2.8272, -0.2716, 0.4130])  #C
        emol.GetConformer().SetAtomPosition(1, [-1.7908, -0.1146, -0.4177]) #C
        emol.GetConformer().SetAtomPosition(2, [-0.5452, 0.6287, -0.0653])  #C
        emol.GetConformer().SetAtomPosition(3, [0.5603, -0.2584, -0.1671])  #O
        emol.GetConformer().SetAtomPosition(4, [1.7601, 0.4902, 0.1811])    #N
        emol.GetConformer().SetAtomPosition(5, [2.8427, -0.4743, 0.0559])   #C
        emol.GetConformer().SetAtomPosition(6, [-3.7065, -0.8216, 0.0959])  #H
        emol.GetConformer().SetAtomPosition(7, [-2.8190, 0.1408, 1.4159])   #H
        emol.GetConformer().SetAtomPosition(8, [-1.8508, -0.5462, -1.4133]) #H
        emol.GetConformer().SetAtomPosition(9, [-0.6006, 1.0355, 0.9521])   #H
        emol.GetConformer().SetAtomPosition(10, [-0.4226, 1.4609, -0.7693]) #H
        emol.GetConformer().SetAtomPosition(11, [1.8283, 1.1371, -0.6054])  #H
        emol.GetConformer().SetAtomPosition(12, [2.8648, -0.9271, -0.9408]) #H
        emol.GetConformer().SetAtomPosition(13, [2.7437, -1.2668, 0.8044])  #H
        emol.GetConformer().SetAtomPosition(14, [3.8013, 0.0259, 0.2227])   #H
        self.mol = Chem.Mol(emol)

    def test1_GetDonor2FeatVects(self):
        '''Case 1: two hydrogens'''
        conf = self.mol.GetConformer(-1)
        case1 = GetDonor2FeatVects(conf, [2], scale=1.5)
        pos_heavy_atom = conf.GetAtomPosition(2)

        #Check if there are two vectors
        self.assertEqual(len(case1[0]), 2, 'Incorrect number of vectors')

        #Check initial points of the vectors
        self.assertListAlmostEqual(case1[0][0][0], pos_heavy_atom,
                'Incorrect starting point of vector 1')
        self.assertListAlmostEqual(case1[0][1][0], pos_heavy_atom,
                'Incorrect starting point of vector 2')

        #Check directions of the vectors
        vec_h1 = conf.GetAtomPosition(9) - pos_heavy_atom
        vec_h2 = conf.GetAtomPosition(10) - pos_heavy_atom
        vec_1 = case1[0][0][1] - case1[0][0][0]
        vec_2 = case1[0][1][1] - case1[0][1][0]

        self.assertListAlmostEqual(vec_1.CrossProduct(vec_h1), Point3D(0,0,0),
                'Incorrect direction of vector 1')
        self.assertTrue(vec_1.DotProduct(vec_h1) > 0,
                'Incorrect direction of vector 1')
        self.assertListAlmostEqual(vec_2.CrossProduct(vec_h2), Point3D(0,0,0),
                'Incorrect direction of vector 2')
        self.assertTrue(vec_2.DotProduct(vec_h2) > 0,
                'Incorrect direction of vector 2')

        #Check length of the vectors
        self.assertAlmostEqual(vec_1.Length(), 1.5,
                msg='Incorrect length of vector 1')

        self.assertAlmostEqual(vec_2.Length(), 1.5,
                msg='Incorrect length of vector 2')

    def test2_1_GetDonor2FeatVects(self):
        '''Case 2.1: one hydrogen with sp2 arrangement'''
        conf = self.mol.GetConformer(-1)

        case21 = GetDonor2FeatVects(conf, [1], scale=1.5)
        pos_heavy_atom = conf.GetAtomPosition(1)

        #Check if there is one vector
        self.assertEqual(len(case21[0]), 1, 'Incorrect number of vectors')

        #Check initial point of the vector
        self.assertListAlmostEqual(case21[0][0][0], pos_heavy_atom,
                'Incorrect starting point of vector')

        #Check direction of the vector
        vec_h = conf.GetAtomPosition(8) - (pos_heavy_atom)
        vec = case21[0][0][1] - case21[0][0][0]
        self.assertListAlmostEqual(vec.CrossProduct(vec_h), Point3D(0,0,0),
                'Incorrect direction of vector')
        self.assertTrue(vec.DotProduct(vec_h) > 0,
                'Incorrect direction of vector')

        #Check length of the vector
        self.assertAlmostEqual(vec.Length(), 1.5,
                msg='Incorrect length of vector')

    def test2_2_GetDonor2FeatVects(self):
        #Case 2.2: one hydrogen with sp3 arrangement
        conf = self.mol.GetConformer(-1)
        case22 = GetDonor2FeatVects(conf, [4], scale=1.5)
        pos_heavy_atom = conf.GetAtomPosition(4)

        #Check if there are two vectors
        self.assertEqual(len(case22[0]), 2, 'Incorrect number of vectors')

        #Check initial points of the vectors
        self.assertListAlmostEqual(case22[0][0][0], pos_heavy_atom,
                'Incorrect starting point of vector 1')
        self.assertListAlmostEqual(case22[0][1][0], pos_heavy_atom,
                'Incorrect starting point of vector 2')

        #Check directions of the vectors
        vec_h = conf.GetAtomPosition(11) - pos_heavy_atom

        vec_nbr1 = conf.GetAtomPosition(3) - pos_heavy_atom
        vec_nbr1.Normalize()
        vec_nbr2 = conf.GetAtomPosition(5) - pos_heavy_atom
        vec_nbr2.Normalize()
        avg_vec = (vec_nbr1 + vec_nbr2)

        vec_1 = case22[0][0][1] - case22[0][0][0]
        vec_2 = case22[0][1][1] - case22[0][1][0]

        self.assertListAlmostEqual(vec_1.CrossProduct(vec_h), Point3D(0,0,0),
                'Incorrect direction of vector 1')
        self.assertTrue(vec_1.DotProduct(vec_h) > 0,
                'Incorrect direction of vector 1')
        self.assertListAlmostEqual(vec_2.CrossProduct(avg_vec), Point3D(0,0,0),
                'Incorrect direction of vector 2')
        self.assertTrue(vec_2.DotProduct(avg_vec) < 0,
                'Incorrect direction of vector 2')

        #Check length of the vectors
        self.assertAlmostEqual(vec_1.Length(), 1.5,
                msg='Incorrect length of vector 1')

        self.assertAlmostEqual(vec_2.Length(), 1.5,
                msg='Incorrect length of vector 2')

    def test3_GetDonor2FeatVects(self):
        '''Case 3: no hydrogens'''
        conf = self.mol.GetConformer(-1)
        case3 = GetDonor2FeatVects(conf, [3], scale=1.5)
        pos_heavy_atom = conf.GetAtomPosition(3)

        #Check if there is one vector
        self.assertEqual(len(case3[0]), 1, 'Incorrect number of vectors')

        #Check initial point of the vector
        self.assertListAlmostEqual(case3[0][0][0], pos_heavy_atom,
                'Incorrect starting point of vector')

        #Check direction of the vector
        vec_nbr1 = conf.GetAtomPosition(2) - pos_heavy_atom
        vec_nbr1.Normalize()
        vec_nbr2 = conf.GetAtomPosition(4) - pos_heavy_atom
        vec_nbr2.Normalize()
        avg_vec = (vec_nbr1 + vec_nbr2)
        vec = case3[0][0][1] - case3[0][0][0]
        self.assertListAlmostEqual(vec.CrossProduct(avg_vec), Point3D(0,0,0),
                'Incorrect direction of vector')
        self.assertTrue(vec.DotProduct(avg_vec) < 0,
                'Incorrect direction of vector')

        #Check length of the vector
        self.assertAlmostEqual(vec.Length(), 1.5,
                msg='Incorrect length of vector')

if __name__ == '__main__':
    unittest.main()
