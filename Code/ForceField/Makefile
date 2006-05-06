include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDBASE)/Code -Wall
SOURCES=ForceField.cpp UFF/AngleBend.cpp UFF/BondStretch.cpp UFF/Nonbonded.cpp UFF/TorsionAngle.cpp UFF/DistanceConstraint.cpp UFF/Params.cpp 
LIB=libForceFields.a

$(LIB): $(OBJS)
	ar -rv $@ $?

testExecs/testForceField.exe: testForceField.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(RDKIT) $(RDOPTIMIZER) $(RDGEOMETRY) $(RDGENERAL) 

clean:
	rm -f $(OBJS) $(LIB) testForceField.o testExecs/testForceField.exe $(DEPENDS)

regrs: testExecs/testForceField.exe

include $(DEPENDS)

