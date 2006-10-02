include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDBASE)/Code -Wall
SOURCES=BoundsMatrixBuilder.cpp Embedder.cpp
LIB=libDistGeomHelpers.a

$(LIB): $(OBJS)
	ar -rv $@ $?

testExecs/main.exe: testDgeomHelpers.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lDistGeomHelpers $(RDFILEPARSE) $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB) $(RDKIT) $(RDDISTGEOM) $(RDFORCEFIELDHELPERS) $(RDFORCEFIELDS) $(RDOPTIMIZER) $(RDALIGNMENT) $(RDEIGENSOLVERS) $(RDGENERAL) $(RDGEOMETRY) 

clean:
	rm -f $(OBJS) $(LIB) testDgeomHelpers.o testExecs/main.exe $(DEPENDS)

regrs: testExecs/main.exe

benchDGeom.exe: benchDGeom.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lDistGeomHelpers $(RDFILEPARSE) $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB) $(RDKIT) $(RDDISTGEOM) $(RDFORCEFIELDHELPERS) $(RDFORCEFIELDS) $(RDOPTIMIZER) $(RDMOLALIGN) $(RDMOLTRANSFORMS) $(RDALIGNMENT) $(RDEIGENSOLVERS) $(RDGENERAL) $(RDGEOMETRY) 

include $(DEPENDS)


