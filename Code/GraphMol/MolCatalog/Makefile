include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) $(VFLIBINC)
LDFLAGS=$(RDCATALOGS) $(RDSMILES) $(RDFILEPARSE) $(RDKIT) $(RDBITVECS) $(VFLIB) $(RDGENERAL) $(RDGEOMETRY)
SOURCES=MolCatalogEntry.cpp MolCatalogParams.cpp

$(RDMOLCATALOG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^

static: $(RDMOLCATALOG_S)

$(RDMOLCATALOG_S): $(OBJS)
	ar -rv $@ $?

install: $(RDMOLCATALOG)
	cp $(RDMOLCATALOG) $(RDBASE)/bin

testExecs/test1.exe: test1.o $(RDMOLCATALOG)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDMOLCATALOG) $(LDFLAGS)

clean:
	rm -f $(OBJS) test1.o $(RDMOLCATALOG) $(RDMOLCATALOG_S) testExecs/test1.exe $(DEPENDS)

regrs: testExecs/test1.exe

runregrs: testExecs/test1.exe
	./testExecs/test1.exe

include $(DEPENDS)
