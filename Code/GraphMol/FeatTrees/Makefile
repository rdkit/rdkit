include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDBASE)/Code -Wall
SOURCES=FeatTreeUtils.cpp FeatTree.cpp
LIB=libFeatTrees.a

$(LIB): $(OBJS)
	ar -r $@ $(OBJS)

testExecs/testFeatTrees.exe: testFeatTrees.o $(LIB)
	g++ $(CXXFLAGS) -o $@ $< -L. -lFeatTrees $(RDSMILES) $(VFLIB) $(RDKIT) $(RDEIGENSOLVERS) $(RDGENERAL) $(RDGEOMETRY) 

clean:
	rm -f $(OBJS) $(LIB) testFeatTrees.o testExecs/testFeatTrees.exe $(DEPENDS)

regrs: testExecs/testFeatTrees.exe

include $(DEPENDS)


