include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make

CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDCODEDIR)

SOURCES=DistPicker.cpp MaxMinPicker.cpp HierarchicalClusterPicker.cpp
LIB=libSimDivPickers.a

$(LIB): $(OBJS)
	ar -rv $@ $?

clean: 
	rm -f $(OBJS) $(LIB) $(DEPENDS)

include $(DEPENDS)
