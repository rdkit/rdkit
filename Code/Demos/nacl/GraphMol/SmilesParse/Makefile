CLIENT:=$(filter-out *%,$(subst /googleclient,/googleclient *,$(CURDIR)))
PLATFORM:=linux
ifeq (nacl,$(filter nacl,$(MAKECMDGOALS)))
CC = $(NACL_BIN_PATH)/nacl-gcc
NACL_INCLUDE = -I$(GOOGLE_CLIENT)/third_party/npapi/files/include
else
 $(error This sample only works with nacl builds.)
endif

CCFLAGS=-Wall -Wno-unused-function -Wno-deprecated -I$(RDBASE)/Code -I$(BOOSTHOME)
include $(CLIENT)/native_client/common/Makefile.mk

objs=lex.yysmiles.o  smiles.tab.o lex.yysmarts.o smarts.tab.o InputFiller.o SmilesParse.o SmilesParseOps.o SmilesWrite.o SmartsWrite.o


lexes=lex.yysmiles.o lex.yysmarts.o
smiles.tab.cpp smiles.tab.hpp : smiles.yy
	bison -p yysmiles_ -d -o smiles.tab.cpp $<
lex.yysmiles.cpp: smiles.ll smiles.tab.hpp
	flex -Pyysmiles_ -o $@ $<
smarts.tab.cpp smarts.tab.hpp : smarts.yy
	bison -p yysmarts_ -d -o smarts.tab.cpp $<
lex.yysmarts.cpp: smarts.ll smarts.tab.hpp
	flex -Pyysmarts_ -o $@ $<


$(objs): %.o: %.cpp
	$(CC) -c $(CCFLAGS) $< -o $@



LIB=$(RDBASE)/bin/libSmilesParse.a
$(LIB): $(objs)
	ar -r $(LIB) $(objs)

lib : $(LIB)