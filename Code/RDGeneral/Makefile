include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(CXXOPTFLAGS) $(BOOSTINC) -I$(RDCODEDIR)


SOURCES=Invariant.cpp types.cpp utils.cpp Dict.cpp RDLog.cpp
LIB=libRDGeneral_s.a
SLIB=libRDGeneral.so 
DEST=$(RDBASE)/bin/$(SLIB)

all: $(LIB) $(SLIB)

static: $(LIB)


$(SLIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^ $(BOOSTLOGLIB)

$(LIB): $(OBJS)
	ar -rv $@ $?

# Dict.o needs to be produced without optimization, otherwise some important
# definitions get compiled out. 
Dict.o: Dict.cpp
	$(CXX) -g -c $(BOOSTINC) -I$(RDCODEDIR) Dict.cpp

testExecs/testDict.exe: testDict.o $(SLIB)
	$(CXX) -o $@ $(CXXFLAGS) testDict.o $(RDGENERAL)

clean:
	rm -f testExecs/testDict.exe testDict.o $(OBJS) $(SLIB) $(LIB) $(DEPENDS)

regrs: testExecs/testDict.exe 

runregrs: regrs
	testExecs/testDict.exe

$(DEST): $(SLIB)
	cp $^ $@

install: $(DEST)

include $(DEPENDS)


