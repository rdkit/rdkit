
from Numeric import *
import linalg

print(linalg.GetFirstElement(array([1, 0, 2], Int)))
print(linalg.GetFirstElement(array([1, 0, 2], Float)))
print(linalg.GetFirstElement(array([1, 0, 2], Float64)))
print(linalg.GetFirstElement(array([1, 0, 2], Complex)))  # will be zero
