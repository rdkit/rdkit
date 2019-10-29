CXXFLAGS=-g -I$(RDBASE)/External/$(BOOSTBASE) -I$(RDBASE)/Code/GraphMol

dict.exe: dict.o $(RDBASE)/Code/GraphMol/Invariant/Invariant.o
	g++ -o dict.exe $(CFLAGS) dict.o $(RDBASE)/Code/GraphMol/Invariant/Invariant.o