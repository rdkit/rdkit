#from DataStructs import cDataStructs
import moduleA
import moduleB
import moduleC
from DataStructs import cDataStructs

print("*****************************")
print("Testing self print for classA")
A = moduleA.classA()
A.printA()
print("*****************************")
print()

print("*****************************")
print("Testing self print for classC")
C = moduleC.classC()
C.printC()
print("*****************************")
print()

print("*****************************")
print("Testing crossing C to B")
moduleB.testCrossC(C)
print("*****************************")
print()

cc = moduleB.GetClassC()
cc.printC()

fp = cDataStructs.ExplicitBitVect(20)

print("created fp")

fp = moduleB.GetEBV()

#print("*****************************")
#print("Testing crossing A to B")
#moduleB.testCrossA(A)
#print("*****************************")
