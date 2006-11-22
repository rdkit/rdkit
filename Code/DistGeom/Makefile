include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDBASE)/Code -Wall
SOURCES=DistGeomUtils.cpp TriangleSmooth.cpp DistViolationContrib.cpp ChiralViolationContrib.cpp
LIB=libDistGeom.a

$(LIB): $(OBJS)
	ar -rv $@ $?

testExecs/main.exe: testDistGeom.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(RDFORCEFIELDS) $(RDOPTIMIZER) $(RDEIGENSOLVERS) $(RDGENERAL) $(RDGEOMETRY)

clean:
	rm -f $(OBJS) $(LIB) testDistGeom.o testExecs/main.exe $(DEPENDS)

regrs: testExecs/main.exe

include $(DEPENDS)

