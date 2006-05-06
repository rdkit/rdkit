include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS)
LDFLAGS= $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

SOURCES=MolTransforms.cpp
LIB=libMolTransforms.a
SLIB=libMolTransforms.so
DEST=$(RDBASE)/bin/$(SLIB)

libs: $(LIB) $(SLIB)

$(LIB): $(OBJS)
	ar -rv $@ $?

static: $(LIB)

$(SLIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^ $(RDEIGENSOLVERS)

testExecs/testMolTransform.exe: test1.o $(LIB)
	g++ -o $@ $(CXXFLAGS) $< -L. -lMolTransforms $(RDFILEPARSE) $(RDSMILES) $(RDKIT) $(RDEIGENSOLVERS) $(RDGEOMETRY) $(RDGENERAL) $(LDFLAGS)

clean:
	rm -f $(OBJS) test1.o testExecs/testMolTransform.exe $(SLIB) $(LIB) $(DEPENDS)

regrs: testExecs/testMolTransform.exe libMolTransforms.a 

$(DEST): $(SLIB)
	cp $(SLIB) $(RDBASE)/bin

install: $(DEST)

include $(DEPENDS)
