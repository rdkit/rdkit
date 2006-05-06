#      LAPACK++ (V. 1.1 Beta)
#      (C) 1992-1996 All Rights Reserved.
#
#   	Lapack++ drivers and factorizations

LAPACKPP_DIR = ..
include ../makefile.def

##############################################################################

LAPACKPP_OBJS = LaLUFactorDouble.o
LAPACKPP_OBJS = 


LAPACKPP_DRIVERS =	linslv.o eigslv.o

LAPACKPP_UTILS = systime.o  dopla.o dtimmg.o   lautil.o


$(LAPACKPP_LIB) : $(LAPACKPP_OBJS) $(LAPACKPP_UTILS) $(LAPACKPP_DRIVERS)
	rm -f $(LAPACKPP_LIB)
	ar cr $(LAPACKPP_LIB) $(LAPACKPP_OBJS) $(LAPACKPP_UTILS) $(LAPACKPP_DRIVERS)
	case x$(HASRANLIB) in xt ) echo ranlib; \
		ranlib $(LAPACKPP_LIB);; esac

$(LAPACKPP_UTILS) : 


# High-level drivers
#
LaLinearSolveDouble.o : $(LAPACKPP_INC)/laexcp.h  \
	$(LAPACKPP_INC)/LaUtil.h $(LAPACKPP_INC)/lapack.h  \
	$(LAPACKPP_INC)/gmd.h \
	$(LAPACKPP_INC)/lavd.h $(LAPACKPP_INC)/spdmd.h \
	$(LAPACKPP_INC)/ltgmd.h \
	$(LAPACKPP_INC)/lavi.h \
	$(LAPACKPP_INC)/linslv.h



clean:
	rm -f *.o core a.out out

wipe:
	rm -f *.o core a.out out $(LAPACKPP_LIB)

.cpp.o: 
	$(CPP) $(PROFILE) $(CPPFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) $(PROFILE) -c $<

.f.o:
	$(F77) $(PROFILE) -c $<
