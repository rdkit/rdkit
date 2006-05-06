#      LAPACK++ (V. 1.1)
#      (C) 1992-1996 All Rights Reserved.
#
#       Test Code for Lapack++ drivers

LAPACKPP_DIR = ../..
include ../../makefile.def


TESTS = tblasd++  

all:  $(TESTS)


tblasd++ : tblasd++.o
	$(CPP) $(CPPFLAGS) -o tblasd++ tblasd++.o $(LAMATRIXPP_LIB) \
	 $(COMPLEX_LIB)  $(BLASPP_LIB) $(BLAS_LIB) \
		 $(FLIBS)  

 
clean:
	rm -f *.o *.out core $(TESTS)

wipe:
	rm -f *.o a.out core $(TESTS)

.cpp.o: 
	$(CPP) $(CPPFLAGS) -c $*.cpp
