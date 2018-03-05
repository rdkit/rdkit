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

objs=point.o Transform2D.o Transform3D.o UniformGrid3D.o GridUtils.o

$(objs): %.o: %.cpp
	$(CC) -c $(CCFLAGS) $< -o $@

LIB=$(RDBASE)/bin/libRDGeometry.a
$(LIB): $(objs)
	ar -r $(LIB) $(objs)

lib : $(LIB)