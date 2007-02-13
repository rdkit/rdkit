include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS) $(LAPACKINC)

SOURCES=Atom.cpp QueryAtom.cpp QueryBond.cpp Bond.cpp MolOps.cpp FindRings.cpp ROMol.cpp RWMol.cpp \
	PeriodicTable.cpp atomic_data.cpp QueryOps.cpp MolPickler.cpp Canon.cpp AtomIterators.cpp \
	BondIterators.cpp Aromaticity.cpp Kekulize.cpp MolDiscriminators.cpp ConjugHybrid.cpp AddHs.cpp \
	RankAtoms.cpp Matrices.cpp Chirality.cpp RingInfo.cpp Conformer.cpp

libs: $(RDGRAPHLIB) $(RDGRAPHLIB_S)

static: $(RDGRAPHLIB_S)

$(RDGRAPHLIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^

$(RDGRAPHLIB_S): $(OBJS)
	ar -rv $@ $^

install: $(RDGRAPHLIB)
	cp $(RDGRAPHLIB) $(RDBASE)/bin

clean:
	rm -f $(RDGRAPHLIB) $(OBJS) test1.o querytest.o test1.exe querytest.exe molopstest.o \
	testExecs/MolOpsTest.exe bulktest.o bulktest.exe testExecs/itertest.exe \
	itertest.o sanitTest.exe testExecs/testCanon.exe \
	sanitTest.o $(DEPENDS)

testExecs/test1.exe: test1.o $(RDGRAPHLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

testExecs/cptest.exe: cptest.o $(RDGRAPHLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

testExecs/querytest.exe: querytest.o $(RDGRAPHLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

testExecs/MolOpsTest.exe: molopstest.cpp $(RDGRAPHLIB) $(RDSMILESLIB)
	$(CXX) $(DEBUGCXXFLAGS) $(LAPACKINC) -o $@ $< \
	$(RDFILEPARSE) $(RDSMILES) $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

bulktest.exe: bulktest.o $(RDGRAPHLIB) $(RDSMILESLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDSMILES) $(RDKIT) $(RDBITVECS)  $(RDGENERAL) $(RDGEOMETRY)

leakfinder.exe: leakfinder.o $(RDGRAPHLIB) $(RDSMILESLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDSMILES) $(RDKIT) $(RDBITVECS)  $(RDGENERAL) $(RDGEOMETRY)

testExecs/itertest.exe: itertest.o $(RDGRAPHLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDSMILES) $(RDKIT) $(RDBITVECS) $(RDGENERAL) $(RDGEOMETRY)

sanitTest.exe: sanitTest.o $(RDGRAPHLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDSMILES) $(RDKIT) $(RDBITVECS) $(RDGENERAL) $(RDGEOMETRY)

testCanon.o: testCanon.cpp
	$(CXX) $(CXXFLAGS) $(VFLIBINC) -c testCanon.cpp

testExecs/testCanon.exe: testCanon.o $(RDGRAPHLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDSMILES) \
	$(RDSUBSTRUCT) $(RDSMILES) $(RDKIT) $(VFLIB) $(RDGENERAL) $(RDGEOMETRY)

testExecs/testPickler.exe: testPickler.o $(RDGRAPHLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDSMILES) \
	$(RDFILEPARSE) $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

benchPickler.exe: benchPickler.o $(RDGRAPHLIB_S)
	$(CXX) $(CXXFLAGS) -o $@ $< $(RDKIT_S) $(RDGENERAL_S)

benchSearch.exe: benchSearch.o $(RDGRAPHLIB_S)
	$(CXX) $(CXXFLAGS) -o $@ $< \
	$(RDSMILES_S) $(RDSUBSTRUCT_S)  $(RDKIT_S) $(RDGENERAL_S) $(VFLIB)

allclean:
	$(MAKE) clean
	$(MAKE) --directory=FileParsers clean	
	$(MAKE) --directory=SmilesParse clean
	$(MAKE) --directory=Substruct clean
	$(MAKE) --directory=Subgraphs clean
	$(MAKE) --directory=Fingerprints clean
	$(MAKE) --directory=FragCatalog clean
	$(MAKE) --directory=Depictor clean
	$(MAKE) --directory=PartialCharges clean
	$(MAKE) --directory=ForceFieldHelpers clean
	$(MAKE) --directory=DistGeomHelpers clean
	$(MAKE) --directory=MolChemicalFeatures clean
	$(MAKE) --directory=MolAlign clean
	$(MAKE) --directory=MolTransforms clean
	$(MAKE) --directory=ShapeHelpers clean
	$(MAKE) --directory=ChemTransforms clean
	$(MAKE) --directory=MolCatalog clean

all: $(RDGRAPHLIB) $(RDGRAPHLIB_S)
	$(MAKE) --directory=SmilesParse install
	$(MAKE) --directory=Substruct install
	$(MAKE) --directory=FileParsers install
	$(MAKE) --directory=Subgraphs install
	$(MAKE) --directory=Fingerprints install
	$(MAKE) --directory=FragCatalog install
	$(MAKE) --directory=Depictor install
	$(MAKE) --directory=PartialCharges 	
	$(MAKE) --directory=MolTransforms install
	$(MAKE) --directory=ForceFieldHelpers
	$(MAKE) --directory=DistGeomHelpers
	$(MAKE) --directory=Descriptors
	$(MAKE) --directory=MolChemicalFeatures
	$(MAKE) --directory=MolAlign
	$(MAKE) --directory=ShapeHelpers
	$(MAKE) --directory=ChemTransforms install
	$(MAKE) --directory=MolCatalog install

allstatic: $(RDGRAPHLIB_S)
	$(MAKE) --directory=SmilesParse static
	$(MAKE) --directory=Substruct static
	$(MAKE) --directory=FileParsers static
	$(MAKE) --directory=Subgraphs static
	$(MAKE) --directory=Fingerprints static
	$(MAKE) --directory=FragCatalog static
	$(MAKE) --directory=Depictor static
	$(MAKE) --directory=PartialCharges 	
	$(MAKE) --directory=MolTransforms static
	$(MAKE) --directory=ForceFieldHelpers
	$(MAKE) --directory=DistGeomHelpers
	$(MAKE) --directory=Descriptors
	$(MAKE) --directory=MolChemicalFeatures
	$(MAKE) --directory=MolAlign
	$(MAKE) --directory=ShapeHelpers
	$(MAKE) --directory=ChemTransforms static
	$(MAKE) --directory=MolCatalog static

wrappers: $(RDGRAPHLIB)
	for d in . FragCatalog PartialCharges Depictor DistGeomHelpers ForceFieldHelpers \
		MolChemicalFeatures MolAlign MolTransforms ShapeHelpers MolCatalog;\
	do								\
		cd $(RDBASE)/Code/GraphMol/$$d/Wrap;						\
		python setup.py build --force && python setup.py install --install-lib=$(RDBASE)/Python; \
	done

regrs: testExecs/test1.exe testExecs/querytest.exe testExecs/MolOpsTest.exe testExecs/itertest.exe testExecs/cptest.exe testExecs/testCanon.exe testExecs/testPickler.exe

runregrs: regrs
	./testExecs/test1.exe
	./testExecs/querytest.exe	
	./testExecs/MolOpsTest.exe
	./testExecs/itertest.exe
	./testExecs/cptest.exe
	./testExecs/testCanon.exe
	./testExecs/testPickler.exe

allregrs: 
	$(MAKE) regrs
	$(MAKE) --directory=SmilesParse regrs
	$(MAKE) --directory=Substruct regrs
	$(MAKE) --directory=Subgraphs regrs
	$(MAKE) --directory=FileParsers regrs
	$(MAKE) --directory=Fingerprints regrs
	$(MAKE) --directory=FragCatalog regrs
	$(MAKE) --directory=Depictor regrs
	$(MAKE) --directory=MolTransforms regrs
	$(MAKE) --directory=ForceFieldHelpers regrs
	$(MAKE) --directory=DistGeomHelpers regrs
	$(MAKE) --directory=Descriptors regrs
	$(MAKE) --directory=MolChemicalFeatures regrs
	$(MAKE) --directory=MolAlign regrs
	$(MAKE) --directory=ShapeHelpers regrs
	$(MAKE) --directory=ChemTransforms regrs
	$(MAKE) --directory=MolCatalog regrs

runallregrs: 
	$(MAKE) allregrs
	$(MAKE) --directory=SmilesParse runregrs
	$(MAKE) --directory=Substruct runregrs
	$(MAKE) --directory=Subgraphs runregrs
	$(MAKE) --directory=FileParsers runregrs
	$(MAKE) --directory=Fingerprints runregrs
	$(MAKE) --directory=FragCatalog runregrs

include $(DEPENDS)
