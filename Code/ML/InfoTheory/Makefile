include $(RDBASE)/Code/rdvars.make
CXXFLAGS=$(RDOPTFLAGS) $(BOOSTINC) -I$(RDCODEDIR)

OBJS=InfoBitRanker.o 

HDRS=InfoBitRanker.h InfoGainFuncs.h CorrMatGenerator.h

libInfoTheory.a: $(OBJS) $(HDRS)
	ar -r libInfoTheory.a $(OBJS)

clean: 
	rm -f $(OBJS) libInfoTheory.a
	rm -rf build

$(RDBASE)/Python/ML/InfoTheory/cEntropy.so: libInfoTheory.a
	python setup.py install --install-lib=$(RDBASE)/Python

wrapper: $(RDBASE)/Python/ML/InfoTheory/cEntropy.so
