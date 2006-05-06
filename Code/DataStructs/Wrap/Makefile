
SOURCES=DataStructs.cpp DataStructs.h wrap_BitOps.cpp wrap_ExplicitBV.cpp wrap_helpers.h \
	wrap_SparseBV.cpp wrap_Utils.cpp

$(RDBASE)/Python/DataStructs/cDataStructs.so: $(SOURCES) ../libDataStructs.a
	python setup.py install --install-lib=$(RDBASE)/Python

clean:
	if test -e $(RDBASE)/Python/DataStructs/cDataStructs.so; then rm $(RDBASE)/Python/DataStructs/cDataStructs.so; fi
	if test -e build; then rm -rf build; fi
