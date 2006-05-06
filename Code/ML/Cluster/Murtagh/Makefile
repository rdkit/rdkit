CXX=g++
CC=gcc
FC=g77
CFLAGS=$(RDOPTFLAGS)
CXXFLAGS=$(CFLAGS)
FFLAGS=-O3

libhc.a: hcdriver.o hc.o
	ar -r libhc.a hcdriver.o hc.o


$(RDBASE)/Python/ML/Cluster/Clustering.so: libhc.a
	python setup.py install --install-lib=$(RDBASE)/Python

wrapper: $(RDBASE)/Python/ML/Cluster/Clustering.so



