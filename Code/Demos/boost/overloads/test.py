import overloads

assert overloads.f1(1)==4
assert overloads.f1(1,4)==5
assert overloads.f1(1.)==4.
assert overloads.f1(1.,5)==6.

assert overloads.f2(1)==4
assert overloads.f2(1,12)==13



