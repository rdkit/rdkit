RM = rm -rf
INCDIR = $(RDBASE)/Code
LIBDIR = $(RDBASE)/lib
FLAGS = -O2 -I$(INCDIR) -W -Wall -Wno-unused-function
LIBS = $(OBJS) \
       -L$(LIBDIR) -lFileParsers -lSmilesParse -lGraphMol -lRDGeneral \
       -lboost_thread -lboost_date_time -lz\
       -Wl,-rpath,$(LIBDIR)

INCS = 
OBJS = 

all:		smi2pdb mol2pdb pdb2smi pdb2mol pdb2pdb

smi2pdb:	smi2pdb.cpp $(OBJS) $(INCS)
		$(CXX) $(FLAGS) -o smi2pdb smi2pdb.cpp $(LIBS)

mol2pdb:	mol2pdb.cpp $(OBJS) $(INCS)
		$(CXX) $(FLAGS) -o mol2pdb mol2pdb.cpp $(LIBS)

pdb2smi:	pdb2smi.cpp $(OBJS) $(INCS)
		$(CXX) $(FLAGS) -o pdb2smi pdb2smi.cpp $(LIBS)

pdb2mol:	pdb2mol.cpp $(OBJS) $(INCS)
		$(CXX) $(FLAGS) -o pdb2mol pdb2mol.cpp $(LIBS)

pdb2pdb:	pdb2pdb.cpp $(OBJS) $(INCS)
		$(CXX) $(FLAGS) -o pdb2pdb pdb2pdb.cpp $(LIBS)

# OBJS
PDBWriter.o:	PDBWriter.cpp PDBWriter.h
		$(CXX) $(FLAGS) -c PDBWriter.cpp

PDBParser.o:	PDBParser.cpp PDBParser.h
		$(CXX) $(FLAGS) -c PDBParser.cpp

PDBSupplier.o:	PDBSupplier.cpp PDBSupplier.h PDBParser.h
		$(CXX) $(FLAGS) -c PDBSupplier.cpp

clean:
		$(RM) smi2pdb mol2pdb pdb2smi pdb2mol pdb2pdb
		$(RM) $(OBJS)

