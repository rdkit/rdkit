include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) $(LAPACKINC) -I$(RDBASE)/Code -Wall

testExecs/main.exe: testMatrices.o
	g++ $(CXXFLAGS) -o $@ $^ $(RDGENERAL)

clean:
	rm -f testMatrices.o testExecs/main.exe $(DEPENDS)

regrs: testExecs/main.exe

include $(DEPENDS)

