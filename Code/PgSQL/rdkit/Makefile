# -------------------------
#     Variables used and default values:
# USE_INCHI=0      # enables InChI functions; requires rdkit built with inchi support
# USE_AVALON=0     # enables avalon fingerprint; requires rdkit built with avalon support
# USE_POPCOUNT=0   # enables use of the CPU's popcount instruction
# USE_THREADS=0    # links against boost.system; required with non-ancient boost versions if inchi is enabled or the rdkit is built with threadsafe SSS
# STATIC_LINK=1    # link against the static RDKit libraries
# -------------------------


RDKIT=$(RDBASE)


ifeq ($(USE_INCHI),1)
  INCHILIBS=-lRDInchiLib -lInchi
  INCHIFLAGS=-I${RDKIT}/External -DBUILD_INCHI_SUPPORT  
  INCHIREGRESS=inchi
  USE_THREADS=1
endif
ifeq ($(USE_AVALON),1)
  AVALONLIBS=-lAvalonLib -lavalon_clib
  AVALONFLAGS=-I${RDKIT}/External -DBUILD_AVALON_SUPPORT  
  AVALONREGRESS=avalon
endif

ifeq ($(USE_POPCOUNT),1)
  POPCOUNTFLAGS=-DUSE_BUILTIN_POPCOUNT -msse4.2
endif

ifeq ($(USE_THREADS),1)
  THREADLIBS=-L${BOOSTHOME}/lib -lboost_system
endif

RDKLIBS       = ${AVALONLIBS} ${INCHILIBS} -lChemReactions -lChemTransforms -lFileParsers -lSmilesParse -lFingerprints -lSubgraphs -lDescriptors -lPartialCharges -lSubstructMatch  -lGraphMol -lDataStructs -lDepictor -lRDGeometryLib -lRDGeneral 

ifeq ($(STATIC_LINK),0)
  SHLIB_LINK += -L${RDKIT}/lib -Wl,-rpath,'${RDKIT}/lib' ${RDKLIBS}
else
  SHLIB_LINK += -L${RDKIT}/lib $(addsuffix _static,${RDKLIBS})
endif

ifndef BOOSTHOME
  BOOSTHOME=/usr/local
endif
SHLIB_LINK += -pthread

SHLIB_LINK += ${THREADLIBS}

PG_CPPFLAGS = -I${BOOSTHOME}/include -I${RDKIT}/Code -DRDKITVER='"007300"' ${AVALONFLAGS} ${INCHIFLAGS} ${POPCOUNTFLAGS}

CPLUSPLUSFLAGS = $(filter-out -fexcess-precision=standard -Wmissing-prototypes -Wdeclaration-after-statement, $(CFLAGS))
CPLUSPLUSFLAGS += -Wno-unused-function 
#CPLUSPLUSFLAGS += -march=native
CPLUSPLUSFLAGS += $(PG_CPPFLAGS)
CPLUSPLUSFLAGS := $(filter-out -fexcess-precision=standard,$(CPLUSPLUSFLAGS))

EXTENSION  = rdkit
EXTVERSION = $(shell grep default_version $(EXTENSION).control | sed -e "s/default_version[[:space:]]*=[[:space:]]*'\([^']*\)'/\1/")
PG_CONFIG  = pg_config
MODULE_big = rdkit
OBJS       = rdkit_io.o mol_op.o bfp_op.o sfp_op.o rxn_op.o rdkit_gist.o low_gist.o guc.o cache.o adapter.o
PGXS       := $(shell $(PG_CONFIG) --pgxs)
PG91 = $(shell $(PG_CONFIG) --version | grep -qE " 8\.| 9\.0" && echo no || echo yes)

ifeq ($(PG91),yes)
all: $(EXTENSION)--$(EXTVERSION).sql

$(EXTENSION)--$(EXTVERSION).sql: $(EXTENSION).sql91.in
	cp $< $@

REGRESS    = rdkit-91 props btree molgist bfpgist-91 sfpgist slfpgist fps reaction ${INCHIREGRESS} ${AVALONREGRESS}
DATA = $(EXTENSION)--$(EXTVERSION).sql
EXTRA_CLEAN = $(EXTENSION)--$(EXTVERSION).sql
else
DATA_built = rdkit.sql
DATA = uninstall_rdkit.sql
REGRESS = rdkit-pre91 props btree molgist bfpgist-pre91 sfpgist slfpgist fps reaction ${INCHIREGRESS}
endif
include $(PGXS)

# save our current CC value:
OCC := $(CC)
# and then change the defn of CC so that we'll link with g++
CC = $(CXX)
# we build .c files using this original defn:
%.o : %.c
	$(OCC) $(CPPFLAGS) -fPIC -c -o $@ $< 

%.o : %.cpp
	$(CXX) $(CPLUSPLUSFLAGS) $(CPPFLAGS) -fPIC -c -o $@ $< 
