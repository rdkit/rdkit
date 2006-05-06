
SOURCES=HierarchicalClusterPicker.cpp MaxMinPicker.cpp rdSimDivPickers.cpp 

$(RDBASE)/Python/SimDivFilters/rdSimDivPickers.so: $(SOURCES) ../libSimDivPickers.a
	python setup.py install --install-lib=$(RDBASE)/Python

clean:
	if test -e $(RDBASE)/Python/SimDivFilters/rdSimDivPickers.so; then rm $(RDBASE)/Python/SimDivFilters/rdSimDivPickers.so; fi
	if test -e build; then rm -rf build; fi;
