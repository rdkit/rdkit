include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS)

SOURCES=FreeChemicalFeature.cpp
LIB=libChemicalFeatures.a
$(LIB): $(OBJS)
	ar -rv $@ $?

regrs: testExecs/main.exe

testExecs/main.exe: testChemicalFeatures.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ testChemicalFeatures.o $(RDCHEMFEATS) $(RDGENERAL)

clean:
	rm -f $(OBJS) $(LIB) $(DEPENDS) testExecs/main.exe testChemicalFeatures.o

include $(DEPENDS)
