include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS)

LDFLAGS=$(RDCHEMTRANSFORMS) $(RDSUBSTRUCT) $(RDSMILES) $(RDFILEPARSELIB) $(RDKIT) $(RDQUERY) $(RDFILEPARSELIB) $(RDBITVECS) $(RDGENERAL) $(RDGEOMETRY) $(VFLIB)
SOURCES=ChemTransforms.cpp

$(RDCHEMTRANSFORMSLIB): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $?

install: $(RDCHEMTRANSFORMSLIB)
	cp $(RDCHEMTRANSFORMSLIB) $(RDBASE)/bin

testExecs/main.exe: testChemTransforms.o $(RDCHEMTRANSFORMSLIB)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

regrs: testExecs/main.exe

clean:
	rm -f $(OBJS) testChemTransforms.o $(RDCHEMTRANSFORMSLIB) $(DEPENDS)

include $(DEPENDS)
