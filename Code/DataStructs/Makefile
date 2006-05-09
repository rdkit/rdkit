include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) -I$(RDBASE)/Code/DataStructs -I$(RDBASE)/Code $(BOOSTINC) -Wall
SOURCES=BitVect.cpp SparseBitVect.cpp ExplicitBitVect.cpp Utils.cpp \
	base64.cpp BitOps.cpp DiscreteDistMat.cpp DiscreteValueVect.cpp
LIB=libDataStructs.a

$(LIB): $(OBJS)
	ar -rv $@ $^

testExecs/test.exe: testDatastructs.o $(LIB) 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(RDGENERAL)

clean:
	rm -f $(OBJS) $(LIB) testDatastructs.o testExecs/test.exe $(DEPENDS)

regrs: testExecs/test.exe

testExecs/bench.exe: bench.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(RDGENERAL)

include $(DEPENDS)
