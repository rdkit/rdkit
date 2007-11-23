include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS)
LDFLAGS= -L. $(RDGENERAL)

SOURCES=test.cpp

testExecs/test.exe: test.o 
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) 

regrs: testExecs/test.exe

clean:
	rm -f $(OBJS) test.o $(DEPENDS)

include $(DEPENDS)
