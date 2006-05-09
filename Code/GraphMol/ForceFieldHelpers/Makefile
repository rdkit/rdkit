include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS) -Wall
SOURCES=UFF/AtomTyper.cpp UFF/Builder.cpp

libForceFieldHelpers.a: $(OBJS)
	ar -rv $@ $?

testExecs/testHelpers.exe: testHelpers.o libForceFieldHelpers.a 
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lForceFieldHelpers $(RDSMILES) $(RDSUBSTRUCT) $(RDDGEOMHELPERS) $(RDDISTGEOM) $(VFLIB) $(RDFILEPARSE) $(RDKIT) $(RDFORCEFIELDS) $(RDOPTIMIZER) $(RDEIGENSOLVERS) $(RDGENERAL) $(RDGEOMETRY)

clean:
	rm -f $(OBJS) libForceFieldHelpers.a testHelpers.o testExecs/testHelpers.exe $(DEPENDS)

regrs: testExecs/testHelpers.exe

include $(DEPENDS)

