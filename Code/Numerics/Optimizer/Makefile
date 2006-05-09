include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDBASE)/Code -Wall
SOURCES=BFGSOpt.cpp LinearSearch.cpp
LIB=libOptimizer.a

libOptimizer.a: $(OBJS)
	ar -rv $@ $?

# NOTE: The -O1 here is associated with Issue 241:
LinearSearch.o: LinearSearch.cpp $(HDRS)
	$(CXX) -c $(CXXFLAGS) -O1 LinearSearch.cpp


testExecs/testOptimizer.exe: testOptimizer.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^  $(RDGENERAL)

clean:
	rm -f $(OBJS) $(LIB) main.o testExecs/testOptimizer.exe $(DEPENDS)

regrs: testExecs/testOptimizer.exe
include $(DEPENDS)


