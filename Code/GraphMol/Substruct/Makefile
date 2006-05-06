include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS) $(VFLIBINC)

SOURCES=SubstructMatch.cpp SubstructUtils.cpp

LDFLAGS=$(RDSUBSTRUCT) $(RDSMILES) $(RDKIT) $(VFLIB) $(RDGENERAL) $(RDGEOMETRY)

LIB=$(RDSUBSTRUCTLIB)
LIB_S=$(RDSUBSTRUCTLIB_S)

libs: $(LIB) $(LIB_S)

static: $(LIB_S)

$(LIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^

$(LIB_S): $(OBJS)
	ar -rv $@ $?

install: $(LIB)
	cp $(LIB) $(RDBASE)/bin

clean:
	rm -f $(LIB) $(LIB_S) $(OBJS) test1.o test2.o testExecs/test1.exe testExecs/test2 $(DEPENDS)

CXXFLAGS=$(BASECXXFLAGS) $(VFLIBINC)

testExecs/test1.exe: test1.o $(LIB) $(LIBHEADERS)
	$(CXX) -o $@ $< $(LDFLAGS)

cmd_match.exe: cmd_match.o $(LIB) $(LIBHEADERS)
	$(CXX) -o $@ $< $(LDFLAGS)

testExecs/test2.exe: test2.o $(LIB) $(LIBHEADERS)
	$(CXX) -o $@ $< $(LDFLAGS)  $(RDCDXML) $(XERCESLIB)

regrs: testExecs/test1.exe

runregrs: testExecs/test1.exe
	./testExecs/test1.exe

bench.exe: bench.o $(LIB) $(LIBHEADERS)
	$(CXX) -o $@ $< $(CXXFLAGS) $(RDSUBSTRUCT_S) $(RDFILEPARSE_S) $(RDSMILES_S) $(RDKIT_S) $(RDGENERAL_S) $(VFLIB)


include $(DEPENDS)
