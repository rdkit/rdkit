import rdGeometry

def bench():
  p1 = rdGeometry.Point2D(1.0, 2.0)
  p2 = rdGeometry.Point2D(3.0, 4.0)
  for _ in range(1000000):
    p3 = p1 + p2
    p4 = p3 * 2.0
    p5 = p4 - p1
  return p5

def bench2():
  for _ in range(1000000):
    p1 = rdGeometry.Point2D(1.0, 2.0)
    p2 = rdGeometry.Point2D(3.0, 4.0)
    p3 = p1 + p2
    p4 = p3 * 2.0
    p5 = p4 - p1
  return p5

def bench3():
  p1 = rdGeometry.Point2D(1.0, 2.0)
  for _ in range(1000000):
    p1 *= 12.
    _ = p1.Length()
  return 1

import pickle
def bench4():
  p1 = rdGeometry.Point2D(1.0, 2.0)
  for _ in range(1000000):
    p2 = pickle.loads(pickle.dumps(p1))
  return 1