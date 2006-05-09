include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDBASE)/Code $(VFLIBINC) -Wall
SOURCES=MolChemicalFeature.cpp MolChemicalFeatureDef.cpp MolChemicalFeatureFactory.cpp FeatureParser.cpp
LIB=libMolChemicalFeatures.a

$(LIB): $(OBJS)
	ar -rv $@ $? 

testExecs/testFeatures.exe: testFeatures.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lMolChemicalFeatures $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB)\
            $(RDKIT) $(RDGENERAL) $(RDGEOMETRY) 

clean:
	rm -f $(OBJS) $(LIB) testFeatures.o testExecs/testFeatures.exe $(DEPENDS)

regrs: testExecs/testFeatures.exe

benchFeatures.exe: benchFeatures.o libMolChemicalFeatures.a 
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lMolChemicalFeatures $(RDSMILES) $(RDFILEPARSE) $(RDSUBSTRUCT) $(VFLIB) $(RDKIT) $(RDGENERAL) $(RDGEOMETRY) 

include $(DEPENDS)

