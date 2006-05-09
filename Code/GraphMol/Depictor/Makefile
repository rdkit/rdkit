include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS)

LDFLAGS=$(RDDEPICTOR) $(RDSMILES) $(RDFILEPARSELIB) $(RDKIT) $(RDQUERY) $(RDFILEPARSELIB) $(RDBITVECS) $(RDGENERAL) $(RDGEOMETRY)
SOURCES=RDDepictor.cpp EmbeddedFrag.cpp DepictUtils.cpp


static: $(RDDEPICTORLIB_S)

$(RDDEPICTORLIB): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $(RDDEPICTORLIB) $(OBJS)

$(RDDEPICTORLIB_S): $(OBJS) 
	ar -rv $@ $?

install: $(RDDEPICTORLIB)
	cp $(RDDEPICTORLIB) $(RDBASE)/bin

testExecs/testDepictor.exe: testDepictor.o $(RDDEPICTORLIB)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

benchDepictor.exe: benchDepictor.o $(RDDEPICTORLIB_S)
	$(CXX) -o $@ $(CXXFLAGS) $< $(RDDEPICTOR_S) $(RDKIT_S) $(RDGENERAL_S)


regrs: testExecs/testDepictor.exe

clean:
	rm -f $(OBJS) testDepictor.o $(RDDEPICTORLIB) $(RDDEPICTORLIB_S) $(DEPENDS)

include $(DEPENDS)
