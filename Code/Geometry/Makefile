include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(CXXOPTFLAGS) $(BOOSTINC) -I$(RDCODEDIR)


SOURCES=point.cpp Transform2D.cpp Transform3D.cpp UniformGrid3D.cpp GridUtils.cpp
LIB=libRDGeometry.a
SLIB=libRDGeometry.so
DEST=$(RDBASE)/bin/$(SLIB)

all: $(SLIB) $(LIB)

static: $(LIB)

$(SLIB): $(OBJS)
	$(MAKE) --directory $(RDBASE)/Code/DataStructs
	$(MAKE) --directory $(RDBASE)/Code/RDGeneral static
	g++ $(CXXFLAGS) $(SOFLAGS) -o $@ $^ $(RDBITVECS) $(RDGENERAL_S)

$(LIB): $(OBJS)
	ar -rv $@ $^

testExecs/testTransforms.exe: testTransforms.o $(SLIB)
	g++ -o $@ $(CXXFLAGS) $^ $(RDBITVECS) $(RDGENERAL) 

testExecs/testGrid.exe : testGrid.o $(LIB)
	g++ -o $@ $(CXXFLAGS) $^ $(RDBITVECS) $(RDGENERAL) 

clean:
	rm -f testTransforms.o $(OBJS) $(SLIB) $(LIB) $(DEPENDS)

regrs: testExecs/testTransforms.exe testExecs/testGrid.exe

runregrs: regrs
	testExecs/testTransforms.exe	
	testExecs/testGrid.exe	

$(DEST): $(SLIB)
	cp $< $@

install: $(DEST)

include $(DEPENDS)


