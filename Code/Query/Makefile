include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS)
LDFLAGS= -L. $(RDGENERAL) -lboost_serialization-gcc$(GCCVERSION)

SOURCES=Query.cpp
LIB=libQuery.a

$(LIB): $(OBJS)
	ar -rv $@ $?

testExecs/test.exe: test.o $(LIB)
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) 

regrs: testExecs/test.exe


clean:
	rm -f $(OBJS) test.o $(LIB) $(DEPENDS)

include $(DEPENDS)
