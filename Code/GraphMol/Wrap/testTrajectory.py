from __future__ import print_function
import os,sys
import unittest

from rdkit import RDConfig

def feq(v1, v2, tol=1.0e-4):
    return abs(v1-v2) < tol

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def testSnapshot(self):
      s = Snapshot([])
      e = False
      try:
        s.GetPoint2D(12)
      except:
        e = True
      self.assertTrue(e)
      s = Snapshot([0.0, 0.0, 0.0])
      e = False
      try:
        s.GetPoint2D(0)
      except:
        e = True
      self.assertTrue(e)

    def testTrajectory2D(self):
      dim = 2
      np = 10
      ns = 5
      traj = Trajectory(dim, np)
      self.assertEqual(traj.Dimension(), dim)
      self.assertEqual(traj.NumPoints(), np)
      c = []
      for i in range(np * dim):
        c.append(float(i))
      for i in range(ns):
        traj.AddSnapshot(Snapshot(c, float(i)))
      self.assertEqual(len(traj), ns)
      e = False
      try:
        traj.GetSnapshot(ns)
      except:
        e = True
      self.assertTrue(e)
      e = False
      try:
        traj.GetSnapshot(0).GetPoint2D(np)
      except:
        e = True
      self.assertTrue(e)
      for i in range(np):
        self.assertAlmostEqual(traj.GetSnapshot(0).GetPoint2D(i).x, float(i * dim))
        self.assertAlmostEqual(traj.GetSnapshot(0).GetPoint2D(i).y, float(i * dim + 1))
        e = False
        try:
          self.assertAlmostEqual(traj.GetSnapshot(0).GetPoint3D(i).z, 0.0)
        except:
          e = True;
        self.assertFalse(e)
      for i in range(ns):
        self.assertAlmostEqual(traj.GetSnapshot(i).GetEnergy(), float(i))
      traj.RemoveSnapshot(0)
      self.assertEqual(len(traj), ns - 1)
      for i in range(ns - 1):
        self.assertAlmostEqual(traj.GetSnapshot(i).GetEnergy(), float(i + 1))
      traj.InsertSnapshot(0, Snapshot(c, 999.0))
      self.assertEqual(len(traj), ns)
      copySnapshot = Snapshot(traj.GetSnapshot(0))
      traj.AddSnapshot(copySnapshot)
      self.assertEqual(len(traj), ns + 1)
      self.assertAlmostEqual(traj.GetSnapshot(0).GetEnergy(), 999.0)
      self.assertAlmostEqual(traj.GetSnapshot(1).GetEnergy(), 1.0)
      self.assertAlmostEqual(traj.GetSnapshot(len(traj) - 1).GetEnergy(), 999.0)
      traj2 = Trajectory(traj);
      self.assertEqual(len(traj), len(traj2))

    def testTrajectory3D(self):
      dim = 3
      np = 10
      ns = 5
      traj = Trajectory(dim, np)
      self.assertEqual(traj.Dimension(), dim)
      self.assertEqual(traj.NumPoints(), np)
      c = []
      for i in range(np * dim):
        c.append(float(i))
      for i in range(ns):
        traj.AddSnapshot(Snapshot(c, float(i)))
      self.assertEqual(len(traj), ns)
      e = False
      try:
        traj.GetSnapshot(ns)
      except:
        e = True
      self.assertTrue(e)
      e = False
      try:
        traj.GetSnapshot(0).GetPoint2D(np)
      except:
        e = True
      self.assertTrue(e)
      for i in range(np):
        self.assertAlmostEqual(traj.GetSnapshot(0).GetPoint3D(i).x, float(i * dim))
        self.assertAlmostEqual(traj.GetSnapshot(0).GetPoint3D(i).y, float(i * dim + 1))
        self.assertAlmostEqual(traj.GetSnapshot(0).GetPoint3D(i).z, float(i * dim + 2))
        if (not i):
          e = False
          try:
            traj.GetSnapshot(0).GetPoint2D(i)
          except:
            e = True
          self.assertTrue(e)
      for i in range(ns):
        self.assertAlmostEqual(traj.GetSnapshot(i).GetEnergy(), float(i))
      traj.RemoveSnapshot(0)
      self.assertEqual(len(traj), ns - 1)
      for i in range(ns - 1):
        self.assertAlmostEqual(traj.GetSnapshot(i).GetEnergy(), float(i + 1))
      traj.InsertSnapshot(0, Snapshot(c, 999.0))
      self.assertEqual(len(traj), ns)
      copySnapshot = Snapshot(traj.GetSnapshot(0))
      traj.AddSnapshot(copySnapshot)
      self.assertEqual(len(traj), ns + 1)
      self.assertAlmostEqual(traj.GetSnapshot(0).GetEnergy(), 999.0)
      self.assertAlmostEqual(traj.GetSnapshot(1).GetEnergy(), 1.0)
      self.assertAlmostEqual(traj.GetSnapshot(len(traj) - 1).GetEnergy(), 999.0)
      traj2 = Trajectory(traj);
      self.assertEqual(len(traj), len(traj2))

    def testReadAmber(self):
      rdbase = os.environ['RDBASE']
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords_bad.trx')
      traj = Trajectory(2, 0)
      ok = False
      try:
        ReadAmberTrajectory(fName, traj)
      except:
        ok = True
      self.assertTrue(ok)
      traj = Trajectory(3, 3)
      ok = False
      try:
        ReadAmberTrajectory(fName, traj)
      except:
        ok = True
      self.assertTrue(ok)
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords_bad2.trx')
      ok = False
      try:
        traj = Trajectory(3, 3)
        ReadAmberTrajectory(fName, traj)
      except:
        ok = True
      self.assertTrue(ok)
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords.trx')
      traj = Trajectory(3, 3)
      ReadAmberTrajectory(fName, traj)
      self.assertEqual(len(traj), 1)
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords2.trx')
      traj = Trajectory(3, 3)
      ReadAmberTrajectory(fName, traj)
      self.assertEqual(len(traj), 2)

    def testReadAmberPython(self):
      # reimplemented the Amber trajectory reader in Python
      # let's check we get the same data as the C++ reader
      # (test for building a trajectory out of Snapshots from Python)
      rdbase = os.environ['RDBASE']
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords2.trx')
      traj = Trajectory(3, 3)
      nCoords = traj.NumPoints() * 3
      nSnapshots = 0
      hnd = open(fName, 'r')
      line = hnd.readline()
      lineNum = 0
      c = []
      i = 0
      while (line):
        lineNum += 1
        if (lineNum > 1):
          tok = line.strip().split()
          j = 0
          while ((i < nCoords) and (j < len(tok))):
            c.append(float(tok[j]))
            j += 1
            i += 1
          if (i == nCoords):
            nSnapshots += 1
            traj.AddSnapshot(Snapshot(c))
            c = []
            i = 0
            line = ' '.join(tok[j:]) + ' '
          else:
            line = ''
        else:
          line = ''
        line += hnd.readline()
      hnd.close()
      self.assertEqual(i, 0)
      self.assertEqual(nSnapshots, 2)
      traj2 = Trajectory(3, 3)
      ReadAmberTrajectory(fName, traj2)
      self.assertEqual(len(traj), len(traj2))
      self.assertEqual(traj.NumPoints(), traj2.NumPoints())
      for snapshotNum in range(len(traj)):
        for pointNum in range(traj.NumPoints()):
          for i in range(3):
            self.assertEqual(traj.GetSnapshot(snapshotNum).GetPoint3D(pointNum)[i],
              traj2.GetSnapshot(snapshotNum).GetPoint3D(pointNum)[i])

    def testReadGromos(self):
      rdbase = os.environ['RDBASE']
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords_bad.trc')
      traj = Trajectory(2, 0)
      ok = False
      try:
        ReadGromosTrajectory(fName, traj)
      except:
        ok = True
      self.assertTrue(ok)
      traj = Trajectory(3, 3)
      ok = False
      try:
        ReadGromosTrajectory(fName, traj)
      except:
        ok = True
      self.assertTrue(ok)
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords_bad2.trc')
      ok = False
      try:
        traj = Trajectory(3, 3)
        ReadGromosTrajectory(fName, traj)
      except:
        ok = True
      self.assertTrue(ok)
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords.trc')
      traj = Trajectory(3, 3)
      ReadGromosTrajectory(fName, traj)
      self.assertEqual(len(traj), 1)
      fName = os.path.join(rdbase, 'Code', 'GraphMol', 'test_data', 'water_coords2.trc')
      traj = Trajectory(3, 3)
      ReadGromosTrajectory(fName, traj)
      self.assertEqual(len(traj), 2)

if __name__=='__main__':
    print("Testing Trajectory wrapper")
    unittest.main()
