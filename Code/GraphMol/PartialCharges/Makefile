include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(BASECXXFLAGS)

SOURCES=GasteigerCharges.cpp GasteigerParams.cpp
LIB=libPartialCharges.a

$(LIB): $(OBJS)
	ar -rv $@ $?

clean:
	rm -f $(OBJS) $(LIB) $(DEPENDS)

include $(DEPENDS)
