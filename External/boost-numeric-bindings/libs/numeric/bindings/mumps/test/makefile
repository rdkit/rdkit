CXXFLAGS =                 -I/home/karlm/local_disk/svn/glas/trunk/glas -I/home/karlm/local_disk/svn/boost/ -I/home/karlm/local_disk/svn/boost-sandbox/ -I/home/karlm/local_disk/ext/mumps/MUMPS_4.6.4/include
LDFLAGS =        
LIBS =  -L/home/karlm/local_disk/svn/glas/trunk/glas/libs/ -lglas -L/home/karlm/local_disk/svn/boost/libs/ -lboost_programoptions -lboost_system -lboost_filesystem   -L/home/karlm/local_disk/ext/mumps/MUMPS_4.6.4/lib -lsmumps -lcmumps -ldmumps -lzmumps -lpord -L/home/karlm/local_disk/ext/mumps/MUMPS_4.6.4/libseq -lmpiseq -lpthread -L/home/karlm/local_disk/ext/lapack/ -llapack-gfortran.4.2 -llapack_tools-gfortran.4.2 -L/home/karlm/local_disk/ext/lib -lf77blas -latlas -lg2c -L/cw/cc/gcc/lib -lstdc++ -lgfortran

MUMPS_UBLAS_DEPS = 

mumps_ublas_d.o : mumps_ublas.cpp $(MUMPS_UBLAS_DEPS)
	@/cw/cc/gcc/bin/g++-4.2 -c -g $(CXXFLAGS) -o mumps_ublas_d.o mumps_ublas.cpp

mumps_ublas.o : mumps_ublas.cpp $(MUMPS_UBLAS_DEPS)
	@/cw/cc/gcc/bin/g++-4.2 -c -DNDEBUG -O3 $(CXXFLAGS) -o mumps_ublas.o mumps_ublas.cpp


mumps_ublas_d: mumps_ublas_d.o
	@/cw/cc/gcc/bin/g++-4.2 -o mumps_ublas_d $(LDFLAGS) mumps_ublas_d.o $(LIBS)

mumps_ublas: mumps_ublas.o
	@/cw/cc/gcc/bin/g++-4.2 -o mumps_ublas $(LDFLAGS) mumps_ublas.o $(LIBS)
