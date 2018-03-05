include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) -DSUPPORT_COMPRESSED_IO
LDFLAGS= $(RDFILEPARSE) $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB) $(RDKIT) $(RDBITVECS) \
         $(RDGENERAL) $(RDGEOMETRY) $(RDDEPICTOR)

SOURCES=iotest.cpp

testExecs/iotest.exe: iotest.o $(RDFILEPARSELIB)
		$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS) -L$(BOOSTHOME)/lib -lboost_iostreams-gcc-mt

clean:
	rm -f iotest.o testExecs/iotest.exe

include $(DEPENDS)
