RDKIT=$(RDBASE)
USE_INCHI=0
STATIC_LINK=1

ifeq ($(USE_INCHI),1)
  INCHILIBS=-lRDInchiLib -lInchi
  INCHIFLAGS=-I${RDKIT}/External -DBUILD_INCHI_SUPPORT  
  INCHIREGRESS=inchi
endif

RDKLIBS       = ${INCHILIBS} -lChemTransforms -lFileParsers -lSmilesParse -lFingerprints -lSubgraphs -lSubstructMatch  -lDescriptors -lPartialCharges -lGraphMol -lDataStructs -lRDGeometryLib -lRDGeneral 

ifeq ($(STATIC_LINK),1)
  SHLIB_LINK += -L${RDKIT}/lib $(addsuffix _static,${RDKLIBS})
else
  SHLIB_LINK += -L${RDKIT}/lib -Wl,-rpath,'${RDKIT}/lib' ${RDKLIBS}
endif

SHLIB_LINK += -pthread
ifndef BOOSTHOME
  BOOSTHOME=/usr/local/include
endif

PG_CPPFLAGS = -I${BOOSTHOME} -I${RDKIT}/Code -DRDKITVER='"007000"' ${INCHIFLAGS} #-DUSE_BUILTIN_POPCOUNT -msse4.2

CPLUSPLUSFLAGS = $(filter-out -fexcess-precision=standard -Wmissing-prototypes -Wdeclaration-after-statement, $(CFLAGS))
CPLUSPLUSFLAGS += -Wno-unused-function 
#CPLUSPLUSFLAGS += -march=native
CPLUSPLUSFLAGS += $(PG_CPPFLAGS)
CPLUSPLUSFLAGS := $(filter-out -fexcess-precision=standard,$(CPLUSPLUSFLAGS))

EXTENSION  = rdkit
EXTVERSION = $(shell grep default_version $(EXTENSION).control | sed -e "s/default_version[[:space:]]*=[[:space:]]*'\([^']*\)'/\1/")
PG_CONFIG  = pg_config
MODULE_big = rdkit
OBJS       = rdkit_io.o mol_op.o bfp_op.o sfp_op.o rdkit_gist.o low_gist.o guc.o cache.o adapter.o
PGXS       := $(shell $(PG_CONFIG) --pgxs)
PG91 = $(shell $(PG_CONFIG) --version | grep -qE " 8\.| 9\.0" && echo no || echo yes)

ifeq ($(PG91),yes)
all: $(EXTENSION)--$(EXTVERSION).sql

$(EXTENSION)--$(EXTVERSION).sql: $(EXTENSION).sql91.in
	cp $< $@

REGRESS    = rdkit-91 props btree molgist bfpgist-91 sfpgist slfpgist fps ${INCHIREGRESS}
DATA = $(EXTENSION)--$(EXTVERSION).sql
EXTRA_CLEAN = $(EXTENSION)--$(EXTVERSION).sql
else
DATA_built = rdkit.sql
DATA = uninstall_rdkit.sql
REGRESS = rdkit-pre91 props btree molgist bfpgist-pre91 sfpgist slfpgist fps ${INCHIREGRESS}
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
