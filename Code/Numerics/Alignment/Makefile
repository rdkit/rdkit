include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) $(LAPACKINC) -I$(RDBASE)/Code -Wall
SOURCES=AlignPoints.cpp 
LIB=libAlignment.a

$(LIB): $(OBJS)
	ar -rv $@ $?

testExecs/main.exe: testAlignment.o $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(RDGENERAL) $(RDGEOMETRY)

clean:
	rm -f $(OBJS) $(LIB) testAlignment.o testExecs/main.exe $(DEPENDS)

regrs: testExecs/main.exe
include $(DEPENDS)


