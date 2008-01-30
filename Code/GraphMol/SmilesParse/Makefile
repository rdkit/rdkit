include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) #-DVERBOSE_CANON

LDFLAGS= $(RDSMILES) $(RDKIT) $(RDBITVECS) $(RDGENERAL) $(RDGEOMETRY)

SOURCES=lex.yysmiles.cpp smiles.tab.cpp InputFiller.cpp SmilesParse.cpp SmilesParseOps.cpp SmilesWrite.cpp smarts.tab.cpp lex.yysmarts.cpp SmartsWrite.cpp

libs: $(RDSMILESLIB) $(RDSMILESLIB_S)

static: $(RDSMILESLIB_S)


$(RDSMILESLIB): $(OBJS)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) -o $@ $^

$(RDSMILESLIB_S): $(OBJS)
	ar -rv $@ $?


INSTLIB=$(RDBASE)/bin/libSmilesParse.so

$(INSTLIB): $(RDSMILESLIB)
	cp $(RDSMILESLIB) $(RDBASE)/bin

install: $(INSTLIB)

main.exe: main.o $(RDSMILESLIB_S) $(RDGRAPHLIB)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

smamain.exe: smamain.o $(INSTLIB) $(RDGRAPHLIB)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

testExecs/test.exe: test.o $(INSTLIB) $(RDGRAPHLIB)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

smatest.o: smatest.cpp
	$(CXX) -c -g $(BOOSTINC) $(RDKITINCFLAGS) $(VFLIBINC) $<

testExecs/smatest.exe: smatest.o $(INSTLIB) $(RDGRAPHLIB)
	$(CXX) -o $@ $(CXXFLAGS) $< $(RDSUBSTRUCT)  $(LDFLAGS)  $(VFLIB)

clean:
	rm -f $(OBJS) $(RDSMILESLIB) $(RDSMILELIB_S) lex.yysmiles.cpp smiles.tab.cpp smiles.tab.h main.exe \
        main.o smatest.o test2.exe test2.o \
	testExecs/smatest.exe testExecs/test.exe $(DEPENDS)

lex.yysmiles.cpp: smiles.flex smiles.tab.hpp
	flex -Pyysmiles_ -olex.yysmiles.cpp smiles.flex

smiles.tab.cpp smiles.tab.hpp: smiles.y 
	bison -p yysmiles_ -d -o smiles.tab.cpp smiles.y

lex.yysmarts.cpp: smarts.flex smarts.tab.hpp
	flex -Pyysmarts_ -olex.yysmarts.cpp smarts.flex

smarts.tab.cpp smarts.tab.hpp: smarts.y 
	bison -p yysmarts_  -d -o smarts.tab.cpp smarts.y

regrs: testExecs/test.exe testExecs/smatest.exe

runregrs: testExecs/test.exe testExecs/smatest.exe
	./testExecs/test.exe
	./testExecs/smatest.exe

$(RDGRAPHLIB):
	$(MAKE) --directory=$(RDCODEDIR)/GraphMol libs

include $(DEPENDS)
