
LDFLAGS = -L${RDBASE}/lib  -lFileParsers_static  \
-lGraphMol_static -lDataStructs_static -lEigenSolvers_static -lRDGeometryLib_static -lRDGeneral_static \
-lboost_thread -lboost_system

ifndef BOOSTHOME
  BOOSTHOME=/usr/local/include
endif
CPPFLAGS = -O2 -I${BOOSTHOME} -I${RDBASE}/Code -Wno-deprecated -I/usr/local/include/eigen3 -DUSE_EIGEN2

%.o : %.cpp PBFRDKit.h
	$(CXX) $(CPPFLAGS)  -c -o $@ $< 

demo.exe: demo.o PBFRDKit.o 
	$(CXX)  -o $@ $< PBFRDKit.o $(LDFLAGS) 
