include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) $(LAPACKINC) -I$(RDBASE)/Code -Wall
SOURCES=PowerEigenSolver.cpp 
LIB=libEigenSolvers.a

$(LIB): $(OBJS)
	ar -rv $@ $^

testExecs/main.exe: testEigenSolvers.o $(LIB)
	$(MAKE) --directory $(RDBASE)/External/Lapack++ libs
	$(CXX) $(CXXFLAGS) -o $@ $^ $(RDGENERAL) $(LAPACKLIB)

clean:
	rm -f $(OBJS) $(LIB) testEigenSolvers.o testExecs/main.exe $(DEPENDS)

regrs: testExecs/main.exe
include $(DEPENDS)


