RDKIT=$(RDBASE)

subdir = contrib/rdkit
top_builddir = ../..
include $(top_builddir)/src/Makefile.global

SHLIB_LINK += -L${RDKIT}/lib  -L${RDKIT}/bin -Wl,-rpath,'${RDKIT}/bin' -lGraphMol -lSmilesParse  -lFingerprints -lSubstructMatch -lDescriptors -lSubgraphs -lDataStructs -lRDGeometryLib -lRDGeneral 
PG_CPPFLAGS = -I$(BOOSTHOME) -I${RDKIT}/Code -DRDKITVER='"001900"' 

MODULE_big = rdkit
OBJS = rdkit_io.o mol_op.o bfp_op.o sfp_op.o rdkit_gist.o low_gist.o guc.o cache.o adapter.o 

DATA_built = rdkit.sql
DATA = uninstall_rdkit.sql
REGRESS = rdkit props btree molgist bfpgist sfpgist slfpgist fps 



include $(top_srcdir)/contrib/contrib-global.mk

CPLUSPLUSFLAGS = $(filter-out -Wmissing-prototypes -Wdeclaration-after-statement, $(CFLAGS))
CPLUSPLUSFLAGS += -Wno-deprecated -Wno-unused-function


%.o : %.cpp
	$(CC) $(CPLUSPLUSFLAGS) $(CPPFLAGS)  -c -o $@ $< 
