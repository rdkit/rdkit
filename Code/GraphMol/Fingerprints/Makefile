include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) -Wno-deprecated $(LAPACKINC)

LDFLAGS= $(RDSMILES) $(RDSUBGRAPHS) $(RDKIT) $(RDBITVECS) $(RDGENERAL) $(RDGEOMETRY)

SOURCES=Fingerprints.cpp

INSTLIB=$(RDBASE)/bin/libFingerprints.so

libs: $(RDFINGERPRINTSLIB) $(RDFINGERPRINTSLIB_S)

static: $(RDFINGERPRINTSLIB_S)


$(RDFINGERPRINTSLIB): $(OBJS)
	g++ $(CXXFLAGS) $(SOFLAGS) -o $@ $^
$(RDFINGERPRINTSLIB_S): $(OBJS)
	ar -rv $@ $?

install: $(INSTLIB) 

$(INSTLIB): $(RDFINGERPRINTSLIB)
	cp $(RDFINGERPRINTSLIB) $(RDBASE)/bin


testExecs/test1.exe: test1.o $(RDFINGERPRINTSLIB)
		g++ -o $@ $(CXXFLAGS) $< $(RDFINGERPRINTSLIB) $(LDFLAGS) 

bench.exe: bench.o $(RDFINGERPRINTSLIB)
		g++ -o $@ $(CXXFLAGS) $< $(RDFINGERPRINTSLIB) $(RDFILEPARSE) $(LDFLAGS)

regrs: testExecs/test1.exe

runregrs: testExecs/test1.exe
	./testExecs/test1.exe


clean:
	rm -f $(OBJS) test1.o testExecs/test1.exe $(RDFINGERPRINTSLIB) $(RDFINGERPRINTSLIB_S) $(DEPENDS)


include $(DEPENDS)

