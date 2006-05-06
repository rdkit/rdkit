#      LAPACK++ (V. 0.9 Beta)
#      (C) 1992-1993 All Rights Reserved.

#   	matrix Class Makefile

LAPACKPP_DIR =../..
include ../../makefile.def

##############################################################################

OBJS =  \
	bmd.o \
	gmd.o \
	gmf.o \
	gmi.o \
	gmli.o \
	lavd.o \
	lavi.o \
	lavli.o\
	ltgmd.o \
	spdbmd.o \
	spdmd.o \
	spdtrmd.o \
	sybmd.o \
	symd.o \
	sytrmd.o \
	trmd.o \
	ultgmd.o \
	utgmd.o \
	uutgmd.o \
	vd.o \
	vf.o \
	vi.o \
	vli.o 

COMPLEX_OBJS =  \
	gmc.o \
	vc.o \
	lavc.o 

LAMATRIXPP_OBJS = $(OBJS) $(COMPLEX_OBJS)
LAMATRIXPP_OBJS = $(OBJS)

$(LAMATRIXPP_LIB) : $(LAMATRIXPP_OBJS)
	rm -f $(LAMATRIXPP_LIB)
	ar cr $(LAMATRIXPP_LIB) $(LAMATRIXPP_OBJS)
	case x$(HASRANLIB) in xt ) echo ranlib; ranlib $(LAMATRIXPP_LIB) ;; esac


vd.o: $(LAPACKPP_INC)/vd.h
vi.o:  $(LAPACKPP_INC)/vi.h

gmd.o.o: vd.o $(LAPACKPP_INC)/laindex.h  $(LAPACKPP_INC)/gmd.h
gmf.o: vf.o $(LAPACKPP_INC)/laindex.h  $(LAPACKPP_INC)/gmf.h
gmi.o: vi.o $(LAPACKPP_INC)/laindex.h  $(LAPACKPP_INC)/gmi.h
lavd.o: vd.o gmd.o  $(LAPACKPP_INC)/vd.h
lavf.o: vf.o gmf.o $(LAPACKPP_INC)/vf.h
lavi.o: vi.o gmi.o $(LAPACKPP_INC)/vi.h

ltgmd.o: gmd.o  $(LAPACKPP_INC)/ltgmd.h
ultgmd.o: gmd.o  $(LAPACKPP_INC)/ultgmd.h
utgmd.o: gmd.o   $(LAPACKPP_INC)/utgmd.h
uutgmd.o: gmd.o  $(LAPACKPP_INC)/uutgmd.h
symd.o: ltgmd.o $(LAPACKPP_INC)/symd.h
spdmd.o: ltgmd.o $(LAPACKPP_INC)/spdmd.h
symd.o: ltgmd.o $(LAPACKPP_INC)/symd.h
bmd.o: gmd.o $(LAPACKPP_INC)/bmd.h
trmd.o: vd.o $(LAPACKPP_INC)/trmd.h
sytrmd.o: vd.o  $(LAPACKPP_INC)/sytrmd.h
sybmd.o: gmd.o $(LAPACKPP_INC)/sybmd.h
spdtrmd.o: vi.o vd.o  $(LAPACKPP_INC)/spdtrmd.h
spdbmd.o: vd.o $(LAPACKPP_INC)/spdbmd.h


#**********
clean:
	rm -f *.o core a.out out

wipe:
	rm -f *.o core a.out out  $(LAMATRIXPP_LIB)

.cpp.o: 
	$(CPP) $(PROFILE) $(CPPFLAGS) -c $<

.c.o:
	$(cc) $(PROFILE) -c $<

.f.o:
	$(F77) $(PROFILE) -c $<
