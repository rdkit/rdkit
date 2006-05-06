include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS) -Wall
SOURCES=Crippen.cpp MolDescriptors.cpp
LIB=libDescriptors.a

$(LIB): $(OBJS)
	ar -r $@ $?

testExecs/testDescriptors.exe: test.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(RDSMILES) $(RDKIT) $(RDGENERAL) $(RDGEOMETRY) $(RDSUBSTRUCT) $(VFLIB)

clean:
	rm -f $(OBJS) $(LIB) test.o testExecs/testDescriptors.exe $(DEPENDS)

regrs: testExecs/testDescriptors.exe

include $(DEPENDS)


