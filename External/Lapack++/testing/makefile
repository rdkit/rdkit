#      LAPACK++ (v. 1.1 Beta)
#      (C) 1992-1996 All Rights Reserved.
#
#   	Test Code for Lapack++ drivers 

LAPACKPP_DIR = ..
include ../makefile.def

##############################################################################


# tSpdSolve and tSymmSolve not fully working
#TESTS = tGenSolve 
TESTS = tGenSolve tSpdSolve  tutil tEigSolve

all :  $(TESTS)


tGenSolve: tGenSolve.o
	$(CPP) $(CPPFLAGS) -o tGenSolve tGenSolve.o $(ALL_LAPACKPP_LIB)  -lm

tSpdSolve: tSpdSolve.o
	$(CPP) $(CPPFLAGS)  -o tSpdSolve tSpdSolve.o $(ALL_LAPACKPP_LIB) -lm

tSymmSolve: tSymmSolve.o
	$(CPP) $(CPPFLAGS) -o tSymmSolve tSymmSolve.o $(ALL_LAPACKPP_LIB)  -lm

tutil: tutil.o
	$(CPP) $(CPPFLAGS) -o tutil tutil.o $(ALL_LAPACKPP_LIB) -lm

tEigSolve: tEigSolve.o
	$(CPP) $(CPPFLAGS) -o tEigSolve tEigSolve.o $(ALL_LAPACKPP_LIB)  -lm

tGenSolve.o  : $(LAPACKPP_INC)/lapack.h \
	$(LAPACKPP_INC)/gmd.h $(LAPACKPP_INC)/lavd.h \
	$(LAPACKPP_INC)/laslv.h  $(LAPACKPP_INC)/vd.h \
	$(LAPACKPP_INC)/genmd.h \
	$(LAPACKPP_INC)/laexcp.h $(LAPACKPP_INC)/blas++.h \
	$(LAPACKPP_INC)/blas1++.h \
	$(LAPACKPP_INC)/blas2++.h \
	$(LAPACKPP_INC)/blas3++.h 

tSpdSolve.o  : $(LAPACKPP_INC)/lapack.h $(LAPACKPP_INC)/spdmd.h \
	$(LAPACKPP_INC)/gmd.h $(LAPACKPP_INC)/lavd.h \
	$(LAPACKPP_INC)/laslv.h  $(LAPACKPP_INC)/vd.h \
	$(LAPACKPP_INC)/genmd.h \
	$(LAPACKPP_INC)/laexcp.h $(LAPACKPP_INC)/blas++.h \
	$(LAPACKPP_INC)/blas1++.h \
	$(LAPACKPP_INC)/blas2++.h \
	$(LAPACKPP_INC)/blas3++.h 


tSymmSolve.o  : $(LAPACKPP_INC)/lapack.h \
	$(LAPACKPP_INC)/gmd.h $(LAPACKPP_INC)/lavd.h \
	$(LAPACKPP_INC)/laslv.h  $(LAPACKPP_INC)/vd.h \
	$(LAPACKPP_INC)/genmd.h \
	$(LAPACKPP_INC)/laexcp.h $(LAPACKPP_INC)/blas++.h \
	$(LAPACKPP_INC)/blas1++.h \
	$(LAPACKPP_INC)/blas2++.h \
	$(LAPACKPP_INC)/blas3++.h 

clean:
	rm -f *.o out core *.out $(TESTS)

wipe:
	rm -f *.o out *.out core  $(TESTS)

.cpp.o: 
	$(CPP) $(PROFILE) $(CPPFLAGS) -c $<

.c.o:
	$(cc) $(PROFILE) -c $<

.f.o:
	$(F77) $(PROFILE) -c $<

