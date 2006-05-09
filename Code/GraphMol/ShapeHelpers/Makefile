include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) $(VFLIBINC) -I$(RDBASE)/Code -Wall
SOURCES=ShapeEncoder.cpp ShapeUtils.cpp
LIB=libShapeHelpers.a

$(LIB): $(OBJS)
	ar -r $@ $?

testExecs/main.exe: testShapeHelpers.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -lShapeHelpers $(RDMOLALIGN) $(RDMOLTRANSFORMS) $(RDFILEPARSE) $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB) $(RDKIT) $(RDEIGENSOLVERS) $(RDALIGNMENT) $(RDGEOMETRY) $(RDBITVECS) $(RDGENERAL) 

clean:
	rm -f $(OBJS) $(LIB) testShapeHelpers.o testExecs/main.exe $(DEPENDS)

regrs: testExecs/main.exe

include $(DEPENDS)

