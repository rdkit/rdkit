include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS)

LDFLAGS= $(RDFILEPARSE) $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB) $(RDKIT) $(RDBITVECS) \
         $(RDGENERAL) $(RDGEOMETRY)
LDFLAGS2= $(RDFILEPARSE) $(RDSMILES) $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

SOURCES=MolFileParser.cpp MolFileStereochem.cpp MolFileWriter.cpp \
  SDFileParser.cpp SDMolSupplier.cpp SmilesMolSupplier.cpp \
  SmilesWriter.cpp SDWriter.cpp TDTWriter.cpp TDTMolSupplier.cpp \
  TplFileParser.cpp

INSTLIB=$(RDBASE)/bin/libFileParsers.so

libs: $(RDFILEPARSELIB) $(RDFILEPARSELIB_S)

install: $(INSTLIB)

$(INSTLIB): $(RDFILEPARSELIB)
	cp $(RDFILEPARSELIB) $(RDBASE)/bin

static: $(RDFILEPARSELIB_S)

$(RDFILEPARSELIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^
$(RDFILEPARSELIB_S): $(OBJS)
	ar -rv $@ $?

testExecs/test1.exe: test1.o $(RDFILEPARSELIB)
		$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

testExecs/testTplParser.exe: testTpls.o $(RDFILEPARSELIB)
		$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

bench.exe: bench.o $(RDFILEPARSELIB_S)
	$(CXX) -o $@ $(CXXFLAGS) $< $(RDFILEPARSE_S) $(RDKIT_S) $(RDGENERAL_S)

testExecs/testMolSupplier.exe: testMolSupplier.o $(RDFILEPARSELIB_S)
		$(CXX) -o $@ $(CXXFLAGS) $< $(RDFILEPARSE_S) $(RDSMILES) $(RDKIT_S) \
                          $(RDGENERAL) $(RDGEOMETRY) $(RDDEPICTOR)

testExecs/testMolWriter.exe: testMolWriter.o $(RDFILEPARSELIB)
		$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS2)

regrs: testExecs/test1.exe testExecs/testMolSupplier.exe testExecs/testMolWriter.exe testExecs/testTplParser.exe

runregrs: testExecs/test1.exe testExecs/testMolSupplier.exe  testExecs/testMolWriter.exe
	./testExecs/test1.exe
	./testExecs/testMolSupplier.exe
	./testExecs/testMolWriter.exe

clean:
	rm -f $(OBJS) main.o testMolSupplier.o testExecs/testMolSupplier.exe testMolWriter.o testExecs/testMolWriter.exe testExecs/testTplParser.exe $(RDFILEPARSELIB) $(RDFILEPARSELIB_S) $(DEPENDS)

include $(DEPENDS)
