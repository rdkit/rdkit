

LDFLAGS = -g -L${RDBASE}/lib  -lSmilesParse_static  \
-lDepictor_static -lSubstructMatch_static -lSubgraphs_static -lGraphMol_static -lDataStructs_static -lRDGeometryLib_static -lRDGeneral_static

ifndef BOOSTHOME
  BOOSTHOME=/usr/local/include
endif
CPPFLAGS = -g -I${BOOSTHOME} -I${RDBASE}/Code -Wno-deprecated $(shell pkg-config --cflags cairo)

%.o : %.cpp ${RDBASE}/Code/GraphMol/MolDrawing/MolDrawing.h ${RDBASE}/Code/GraphMol/MolDrawing/DrawingToCairo.h ${RDBASE}/Code/GraphMol/MolDrawing/DrawingToSVG.h
	$(CXX) $(CPPFLAGS)  -c -o $@ $< 

demo.exe: demo.o 
	$(CXX)  -o $@ $< $(LDFLAGS) $(shell pkg-config --libs cairo)
