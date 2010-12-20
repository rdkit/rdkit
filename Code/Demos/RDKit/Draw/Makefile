

LDFLAGS = -L${RDBASE}/lib  -lSmilesParse_static  \
-lDepictor_static -lSubstructMatch_static -lSubgraphs_static -lGraphMol_static -lDataStructs_static -lRDGeometryLib_static -lRDGeneral_static

ifndef BOOSTHOME
  BOOSTHOME=/usr/local/include
endif
CPPFLAGS = -I${BOOSTHOME} -I${RDBASE}/Code -Wno-deprecated

%.o : %.cpp MolDrawing.h
	$(CXX) $(CPPFLAGS)  -c -o $@ $< 

demo.exe: demo.o 
	$(CXX) -o $@ $< $(LDFLAGS)
