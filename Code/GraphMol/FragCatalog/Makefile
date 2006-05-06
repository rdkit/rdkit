include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) $(VFLIBINC)

SOURCES=FragCatalogUtils.cpp FragCatGenerator.cpp FragCatalogEntry.cpp FragCatParams.cpp \
        FragFPGenerator.cpp

$(RDFRAGCATLIB): $(OBJS)
	g++ $(CXXFLAGS) $(SOFLAGS) -o $@ $^

static: $(RDFRAGCATLIB_S)

$(RDFRAGCATLIB_S): $(OBJS)
	ar -rv $@ $?

install: $(RDFRAGCATLIB)
	cp $(RDFRAGCATLIB) $(RDBASE)/bin

testExecs/test1.exe: test1.o $(RDFRAGCATLIB)
	g++ $(CXXFLAGS) -o $@ $< \
	$(RDFRAGCAT) $(RDCATALOGS) $(RDSUBGRAPHS) $(RDSUBSTRUCT) \
	$(RDSMILES) $(RDFILEPARSE)\
	$(RDKIT) $(RDBITVECS) $(LAPACKLIB) $(LDFLAGS) $(VFLIB) $(RDGENERAL) $(RDGEOMETRY)

clean:
	rm -f $(OBJS) test1.o $(RDFRAGCATLIB) $(RDFRAGCATLIB_S) testExecs/test1.exe $(DEPENDS)

regrs: testExecs/test1.exe

runregrs: testExecs/test1.exe
	./testExecs/test1.exe

include $(DEPENDS)
