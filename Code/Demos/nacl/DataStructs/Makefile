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

.cpp.o:

objs=BitVect.o SparseBitVect.o ExplicitBitVect.o Utils.o base64.o BitOps.o DiscreteDistMat.o DiscreteValueVect.o

$(objs): %.o: %.cpp
	$(CC) -c $(CCFLAGS) $< -o $@

LIB=$(RDBASE)/bin/libDataStructs.a
$(LIB): $(objs)
	ar -r $(LIB) $(objs)

lib : $(LIB)