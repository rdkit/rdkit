include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make
CXXFLAGS=$(BASECXXFLAGS) $(VFLIBINC)

SOURCES=Catalog.cpp CatalogParams.cpp CatalogEntry.cpp
LIB=$(RDCATALOGLIB)

$(LIB): $(OBJS) 
	ar -rv $@ $?
clean:
	rm -f $(OBJS) $(LIB) $(DEPENDS)

include $(DEPENDS)

