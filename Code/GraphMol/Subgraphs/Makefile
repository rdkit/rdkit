include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) -Wno-deprecated $(LAPACKINC)

LDFLAGS= $(RDSUBGRAPHS) $(RDKIT) $(RDBITVECS) $(LAPACKLIB) $(RDGENERAL) $(RDGEOMETRY)

SOURCES=Subgraphs.cpp SubgraphUtils.cpp

libs: $(RDSUBGRAPHSLIB) $(RDSUBGRAPHSLIB_S)

static: $(RDSUBGRAPHSLIB_S)

$(RDSUBGRAPHSLIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^
$(RDSUBGRAPHSLIB_S): $(OBJS)
	ar -rv $@ $?


install: $(RDSUBGRAPHSLIB)
	cp $(RDSUBGRAPHSLIB) $(RDBASE)/bin

testExecs/test1.exe: test1.o $(RDSUBGRAPHSLIB)
		$(CXX) -o $@ $(CXXFLAGS) $< $(RDSMILES) $(LDFLAGS)

testExecs/test2.exe: test2.o $(RDSUBGRAPHSLIB)
		$(CXX) -o $@ $(CXXFLAGS) $< $(RDSMILES) $(LDFLAGS) 

regrs: testExecs/test1.exe testExecs/test2.exe


runregrs: testExecs/test1.exe testExecs/test2.exe
	./testExecs/test1.exe
	./testExecs/test2.exe

clean:
	rm -f $(OBJS) test1.o test2.o ./testExecs/test1.exe test2.exe $(RDSUBGRAPHSLIB) $(RDSUBGRAPHSLIB_S) $(DEPENDS)

include $(DEPENDS)

