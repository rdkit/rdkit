include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) $(VFLIBINC) -I$(RDBASE)/Code -Wall
SOURCES=AlignMolecules.cpp 
LIB=libMolAlign.a

$(LIB): $(OBJS)
	ar -rv $@ $?

testExecs/testMolAlign.exe: testMolAlign.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lMolAlign $(RDMOLTRANSFORMS) $(RDFILEPARSE) \
            $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB) $(RDDGEOMHELPERS) $(RDDISTGEOM) \
            $(RDFORCEFIELDHELPERS) $(RDFORCEFIELDS) $(RDOPTIMIZER) $(RDKIT) \
            $(RDALIGNMENT) $(RDEIGENSOLVERS) $(RDGENERAL) $(RDGEOMETRY) 

clean:
	rm -f $(OBJS) $(LIB) testMolAlign.o testExecs/testMolAlign.exe $(DEPENDS)

regrs: testExecs/testMolAlign.exe

include $(DEPENDS)

