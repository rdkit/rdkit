
#      LAPACK++ (V. 1.1)
#      (C) 1992-1996 All Rights Reserved.

#   	BLAS++ Library Makefile

LAPACKPP_DIR = ../..
include ../../makefile.def

##############################################################################

BLASPP_OBJS =  blas1++.o blas2++.o blas3++.o


$(BLASPP_LIB) : $(BLASPP_OBJS)
	rm -f $(BLASPP_LIB)
	ar cr $(BLASPP_LIB) $(BLASPP_OBJS)
	case x$(HASRANLIB) in xt ) echo ranlib; ranlib $(BLASPP_LIB) ;; esac
 



#**********
clean:
	rm -f *.o core a.out out

wipe:
	rm -f *.o core a.out out $(BLASPP_LIB)

.cpp.o: 
	$(CPP) $(PROFILE) $(CPPFLAGS) -c $<

.c.o:
	$(cc) $(PROFILE) -c $<

.f.o:
	$(F77) $(PROFILE) -c $<
