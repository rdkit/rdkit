#      LAPACK++ (V. 0.9 Beta)
#      (C) 1992-1993 All Rights Reserved.

#   	Matrix Testing Makefile

LAPACKPP_DIR = ../..
include ../../makefile.def


##############################################################################

TESTS = tgd tld tuld tud tuud tvi ttdd tvd tstd trvd tcvd tbd
TESTS_COMPLEX = tgc tvc

all : $(TESTS)


tgd : tgd.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tgd tgd.o $(LAMATRIXPP_LIB)

tgc : tgc.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tgc tgc.o $(LAMATRIXPP_LIB) -lm $(COMPLEX_LIB)

tld : tld.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tld tld.o $(LAMATRIXPP_LIB)

tuld : tuld.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tuld tuld.o $(LAMATRIXPP_LIB)

tud : tud.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tud tud.o $(LAMATRIXPP_LIB)

# needs to be fixed, currently dumping core
tbd : tbd.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tbd tbd.o $(LAMATRIXPP_LIB)

tuud : tuud.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tuud tuud.o $(LAMATRIXPP_LIB)

tsd : tsd.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tsd tsd.o $(LAMATRIXPP_LIB)


tvd : tvd.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tvd tvd.o $(LAMATRIXPP_LIB)

tvc : tvc.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tvc tvc.o $(LAMATRIXPP_LIB) -lm $(COMPLEX_LIB)

tsb : tsb.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tsb tsb.o $(LAMATRIXPP_LIB)

# also untested.  do not use
tsbd : tsbd.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tsbd tsbd.o $(LAMATRIXPP_LIB)


tgi : tgi.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tgi tgi.o $(LAMATRIXPP_LIB)

tvi : tvi.o $(LAMATRIXPP_LIB)
	$(CPP) $(CPPFLAGS) -o tvi tvi.o $(LAMATRIXPP_LIB)

ttdd : ttdd.o $(LAMATRIXPP_LIB) 
	$(CPP) $(CPPFLAGS) -o ttdd ttdd.o $(LAMATRIXPP_LIB) 

tstd : tstd.o $(LAMATRIXPP_LIB) 
	$(CPP) $(CPPFLAGS) -o tstd tstd.o $(LAMATRIXPP_LIB) 

tcvd : tcvd.o $(LAMATRIXPP_LIB) 
	$(CPP) $(CPPFLAGS) -o tcvd tcvd.o $(LAMATRIXPP_LIB) 
 
trvd : trvd.o $(LAMATRIXPP_LIB) 
	$(CPP) $(CPPFLAGS) -o trvd trvd.o $(LAMATRIXPP_LIB) 
 
 
templatev: templatev.o
	$(CPP) $(CPPFLAGS) -o templatev templatev.o 




#**********
clean:
	rm -f *.o out core a.out  *.out $(TESTS)

wipe:
	rm -f *.o out core a.out  *.out $(TESTS)
	
.cpp.o: 
	$(CPP) $(PROFILE) $(CPPFLAGS) -c $<

.c.o:
	$(CC) $(PROFILE) -c $<

.f.o:
	$(F77) $(PROFILE) -c $<
