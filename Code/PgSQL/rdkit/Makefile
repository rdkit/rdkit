RDKIT=$(RDBASE)


#SHLIB_LINK += -L${RDKIT}/lib -Wl,-rpath,'${RDKIT}/lib' -lGraphMol_static -lSmilesParse_static  -lFingerprints_static -lSubstructMatch_static -lDescriptors_static -lSubgraphs_static -lDataStructs_static -lRDGeometryLib_static -lRDGeneral_static
SHLIB_LINK += -L${RDKIT}/lib -Wl,-rpath,'${RDKIT}/lib' -lGraphMol -lSmilesParse  -lFingerprints -lSubstructMatch -lDescriptors -lSubgraphs -lDataStructs -lRDGeometryLib -lRDGeneral

ifndef BOOSTHOME
  BOOSTHOME=/usr/local/include
endif
PG_CPPFLAGS = -I${BOOSTHOME} -I${RDKIT}/Code -DRDKITVER='"002000"' 

MODULE_big = rdkit
OBJS = rdkit_io.o mol_op.o bfp_op.o sfp_op.o rdkit_gist.o low_gist.o guc.o cache.o adapter.o 
PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)

DATA_built = rdkit.sql
DATA = uninstall_rdkit.sql
REGRESS = rdkit props btree molgist bfpgist sfpgist slfpgist fps 
#REGRESS = rdkit

include $(PGXS)

CPLUSPLUSFLAGS = $(filter-out -Wmissing-prototypes -Wdeclaration-after-statement, $(CFLAGS))
CPLUSPLUSFLAGS += -Wno-deprecated -Wno-unused-function
CPLUSPLUSFLAGS += $(PG_CPPFLAGS)

%.o : %.cpp
	$(CC) $(CPLUSPLUSFLAGS) $(CPPFLAGS)  -c -o $@ $< 
