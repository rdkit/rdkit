include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS)
LDFLAGS= $(RDSMILES) $(RDSUBSTRUCT) $(RDKIT) $(RDGENERAL) $(RDGEOMETRY)

SOURCES=Reaction.cpp MDLParser.cpp DaylightParser.cpp
LIB=$(RDCHEMREACTION)
LIB_S=$(RDCHEMREACTION_S)


libs: $(LIB) $(LIB_S)

$(LIB_S): $(OBJS)
	ar -rv $@ $?

static: $(LIB_S)

$(LIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^ 

testExecs/testReaction.exe: testReaction.o $(LIB)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIB) $(LDFLAGS) $(RDSUBSTRUCT) $(RDFILEPARSE) $(RDSMILES) $(RDKIT) $(RDGEOMETRY) $(RDGENERAL) $(VFLIB)

clean:
	rm -f $(OBJS) testReaction.o testExecs/testReaction.exe $(LIB) $(LIB_S) $(DEPENDS)

regrs: testExecs/testReaction.exe $(LIB_S)

install: $(LIB)
	cp $(LIB) $(RDBASE)/bin

include $(DEPENDS)
