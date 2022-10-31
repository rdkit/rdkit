# -------------------------
#     Variables used and default values:
# USE_INCHI enables InChI functions; requires rdkit built with inchi support:
#USE_INCHI=0
# USE_AVALON enables the avalon fingerprint; requires rdkit built with avalon support:
#USE_AVALON=0
# USE_THREADS links against boost.system; required with non-ancient boost versions if inchi is enabled or the rdkit is built with threadsafe SSS:
#USE_THREADS=0
# STATIC_LINK links against the static RDKit libraries:
#STATIC_LINK=1
# -------------------------


RDKIT=$(RDBASE)


ifeq ($(USE_INCHI),1)
  INCHILIBS=-lRDKitRDInchiLib -lRDKitInchi
  INCHIFLAGS=-I${RDKIT}/External -DBUILD_INCHI_SUPPORT
  INCHIREGRESS=inchi
  USE_THREADS=1
endif
ifeq ($(USE_AVALON),1)
  AVALONLIBS=-lAvalonLib -lavalon_clib
  AVALONFLAGS=-I${RDKIT}/External -DBUILD_AVALON_SUPPORT
  AVALONREGRESS=avalon
endif

POPCOUNTFLAGS=-DRDK_OPTIMIZE_POPCNT -mpopcnt

ifeq ($(USE_THREADS),1)
  THREADLIBS=-L${BOOSTHOME}/lib -lboost_thread -lboost_system
endif

RDKLIBS       = ${AVALONLIBS} ${INCHILIBS} -lRDKitMolDraw2D -lRDKitMolTransforms -lRDKitMolHash -lRDKitFMCS -lRDKitChemReactions -lRDKitChemTransforms -lRDKitFileParsers -lRDKitSmilesParse -lRDKitFingerprints -lRDKitSubgraphs -lRDKitDescriptors -lRDKitPartialCharges -lRDKitSubstructMatch  -lRDKitGraphMol -lRDKitEigenSolvers -lRDKitDataStructs -lRDKitDepictor -lRDKitRDGeometryLib -lRDKitRDGeneral

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
OBJS       = rdkit_io.o mol_op.o bfp_op.o sfp_op.o rxn_op.o rdkit_gist.o bfp_gist.o bfp_gin.o low_gist.o guc.o cache.o adapter.o bitstring.o
PGXS       := $(shell $(PG_CONFIG) --pgxs)
PG95 = $(shell $(PG_CONFIG) --version | grep -qE " 9\.[1-4]" && echo no || echo yes)
PG92 = $(shell $(PG_CONFIG) --version | grep -qE " 9\.[1-2]" && echo yes || echo no)

all: $(EXTENSION)--$(EXTVERSION).sql

ifeq ($(PG95),yes)
$(EXTENSION)--$(EXTVERSION).sql: $(EXTENSION).sql.in
	sed -e's/@RDKIT_PG_BFP_GIST_FETCH@/FUNCTION    9   gbfp_fetch (internal),/;' \
-e's/@RDKIT_GIN_BFP_TRICONSISTENT@/FUNCTION    6   gin_bfp_triconsistent(internal, int2, bfp, int4, internal, internal, internal),/' \
 $< > $@
	cp $< $@
else
ifeq ($(PG92),no)
$(EXTENSION)--$(EXTVERSION).sql: $(EXTENSION).sql.in
	sed -e's/@RDKIT_PG_BFP_GIST_FETCH@//;' \
-e's/@RDKIT_GIN_BFP_TRICONSISTENT@/FUNCTION    6   gin_bfp_triconsistent(internal, int2, bfp, int4, internal, internal, internal),/' \
$< > $@
else
$(EXTENSION)--$(EXTVERSION).sql: $(EXTENSION).sql.in
	sed -e's/@RDKIT_PG_BFP_GIST_FETCH@//;' \
-e's/@RDKIT_GIN_BFP_TRICONSISTENT@//;' \
$< > $@
endif
endif

REGRESS    = rdkit-91 props btree molgist bfpgist-91 bfpgin sfpgist slfpgist fps reaction ${INCHIREGRESS} ${AVALONREGRESS}
DATA = $(EXTENSION)--$(EXTVERSION).sql
EXTRA_CLEAN = $(EXTENSION)--$(EXTVERSION).sql
include $(PGXS)

# save our current CC value:
OCC := $(CC)
# and then change the defn of CC so that we'll link with g++
CC = $(CXX)
# we build .c files using this original defn:
%.o : %.c
	$(OCC) $(CFLAGS) $(CPPFLAGS) -fPIC -c -o $@ $<

%.o : %.cpp
	$(CXX) $(CPLUSPLUSFLAGS) $(CPPFLAGS) -fPIC -c -o $@ $<


COMPILE.cxx.bc = $(CLANG) -std=c++17 -xc++ -Wno-ignored-attributes $(BITCODE_CXXFLAGS) $(CPPFLAGS) -emit-llvm -c

%.bc : %.cpp
	$(COMPILE.cxx.bc) -o $@ $<
	$(LLVM_BINPATH)/opt -module-summary -f $@ -o $@
