include $(RDBASE)/Code/rdvars.make
CXXFLAGS=$(BASECXXFLAGS)

LDFLAGS= $(RDFILEPARSE) $(RDSMILES) $(RDSUBSTRUCT) $(VFLIB) $(RDKIT) $(RDBITVECS) \
         $(RDMOLTRANSFORMS) $(RDGENERAL) $(RDGEOMETRY)

OBJS=main.o TemplEnumTools.o

main.exe: $(OBJS)
	g++ -o main.exe $(CXXFLAGS) $(OBJS) $(LDFLAGS)

SIMPLEOBJS=simple.o TemplEnumTools.o
simple.exe: $(SIMPLEOBJS)
	g++ -o simple.exe $(CXXFLAGS) $(SIMPLEOBJS) $(LDFLAGS)

clean:
	rm -f main.exe $(OBJS)

