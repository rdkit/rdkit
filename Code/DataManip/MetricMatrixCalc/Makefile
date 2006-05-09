include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDCODEDIR)

regrs: testExecs/testMatCalc.exe

testExecs/testMatCalc.exe: testMatCalc.o
	$(CXX) -o $@ $(CXXFLAGS) $^ $(RDGENERAL)

clean:
	rm -f  testMatCalc.o testExecs/testMatCalc.exe $(DEPENDS)

include $(DEPENDS)
