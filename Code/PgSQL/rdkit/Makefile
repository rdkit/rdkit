RDKIT=$(RDBASE)
USE_INCHI=0

ifeq ($(USE_INCHI),1)
  INCHISTATICLIBS=-lRDInchiLib_static -lInchi_static 
  INCHILIBS=-lRDInchiLib -lInchi
  INCHIFLAGS=-I${RDKIT}/External -DBUILD_INCHI_SUPPORT  
  INCHIREGRESS=inchi
endif

RDKSTATICLIBS=-lFileParsers_static -lSmilesParse_static -lFingerprints_static -lSubgraphs_static -lSubstructMatch_static  -lDescriptors_static -lPartialCharges_static -lGraphMol_static -lDataStructs_static -lRDGeometryLib_static -lRDGeneral_static
RDKLIBS       = -lDescriptors -lGraphMol -lSmilesParse  -lFingerprints -lSubstructMatch -lDescriptors -lPartialCharges -lSubgraphs -lDataStructs -lRDGeometryLib -lRDGeneral

SHLIB_LINK += -L${RDKIT}/lib ${INCHISTATICLIBS} ${RDKSTATICLIBS}
#SHLIB_LINK += -L${RDKIT}/lib -Wl,-rpath,'${RDKIT}/lib' ${INCHILIBS} ${RDKLIBS} 

ifndef BOOSTHOME
  BOOSTHOME=/usr/local/include
endif
  
PG_CPPFLAGS = -I${BOOSTHOME} -I${RDKIT}/Code -DRDKITVER='"003100"' ${INCHIFLAGS}

CPLUSPLUSFLAGS = $(filter-out -Wmissing-prototypes -Wdeclaration-after-statement, $(CFLAGS))
CPLUSPLUSFLAGS += -Wno-deprecated -Wno-unused-function
CPLUSPLUSFLAGS += -march=native
CPLUSPLUSFLAGS += $(PG_CPPFLAGS)

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
	$(OCC) $(CPPFLAGS) -march=native -fPIC -c -o $@ $< 

%.o : %.cpp
	$(CXX) $(CPLUSPLUSFLAGS) $(CPPFLAGS)  -c -o $@ $< 
