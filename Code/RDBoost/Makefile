include $(RDBASE)/Code/rdvars.make
include $(RDBASE)/Code/rdrules.make


INCS= -I$(RDBASE)/Code $(BOOSTINC) $(PYTHONINC)
CXXFLAGS=$(RDOPTFLAGS) -DBOOST_PYTHON_DYNAMIC_LIB -ftemplate-depth-100 -Wall $(INCS)
SOURCES=Wrap.cpp
LIB=libRDBoost.a

libRDBoost.a: $(OBJS)
	ar -rv $@ $?

wrapper: $(RDBASE)/Python/rdBase.so

$(RDBASE)/Python/rdBase.so: libRDBoost.a RDBase.cpp
	python setup.py install --install-lib=$(RDBASE)/Python


clean:
	rm -f $(OBJS) $(LIB) $(DEPENDS)
	rm -rf build

include $(DEPENDS)
