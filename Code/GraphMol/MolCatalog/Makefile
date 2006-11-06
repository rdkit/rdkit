include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) $(VFLIBINC)
LDFLAGS=$(RDCATALOGS) $(RDSMILES) $(RDFILEPARSE) $(RDKIT) $(RDBITVECS) $(VFLIB) $(RDGENERAL) $(RDGEOMETRY)
SOURCES=MolCatalogEntry.cpp MolCatalogParams.cpp

$(RDMOLCATALOGLIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^

static: $(RDMOLCATALOGLIB_S)

$(RDMOLCATALOGLIB_S): $(OBJS)
	ar -rv $@ $?

install: $(RDMOLCATALOGLIB)
	cp $(RDMOLCATALOGLIB) $(RDBASE)/bin

testExecs/test1.exe: test1.o $(RDMOLCATALOGLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDMOLCATALOG) $(LDFLAGS)

clean:
	rm -f $(OBJS) test1.o $(RDMOLCATALOGLIB) $(RDMOLCATALOGLIB_S) testExecs/test1.exe $(DEPENDS)

regrs: testExecs/test1.exe

runregrs: testExecs/test1.exe
	./testExecs/test1.exe

include $(DEPENDS)
