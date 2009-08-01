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

objs=Atom.o QueryAtom.o QueryBond.o Bond.o MolOps.o \
  FindRings.o ROMol.o RWMol.o PeriodicTable.o atomic_data.o \
  QueryOps.o MolPickler.o Canon.o AtomIterators.o BondIterators.o \
  Aromaticity.o Kekulize.o MolDiscriminators.o ConjugHybrid.o \
  AddHs.o RankAtoms.o Matrices.o Chirality.o RingInfo.o \
  Conformer.o 

$(objs): %.o: %.cpp
	$(CC) -c $(CCFLAGS) $< -o $@

LIB=$(RDBASE)/bin/libGraphMol.a
$(LIB): $(objs)
	ar -r $(LIB) $(objs)

lib : $(LIB)