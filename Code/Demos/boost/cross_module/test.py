from __future__ import print_function
import moda
import modb

ca = moda.ClassA()
print('ca:',ca.Get4())

cb = modb.ClassB()
print('cb:',cb.Get3())

newca = cb.ReturnOther()
print('new:',newca.Get4())

print('arg:',cb.AcceptOther(newca))
