#      LAPACK++ (V. 1.1)
#      (C) 1992-1996 All Rights Reserved.
#
#		L A P A C K + +     M A K E F I L E
#

LAPACKPP_DIR = .
include makefile.def

error:
	@echo "Usage: make all"
	@echo "        "
	@echo "  Make sure the system-specific makefile.def has been edited "
	@echo "  to reflect your system configuration."
	@echo "        "

LIBS =   $(LAMATRIXPP_LIB) $(BLASPP_LIB) $(LAPACKPP_LIB)  
TESTS = matrix_testing blas++_testing lapack++_testing 

test: $(TESTS)
libs: $(LIBS)
all: $(LIBS) $(TESTS) 

$(LAMATRIXPP_LIB) :
	cd $(LAPACKPP_DIR)/matrix/src; make

# (optionally) make the Matrix Tester

matrix_testing :
	cd $(LAPACKPP_DIR)/matrix/testing; make;
	cd $(LAPACKPP_DIR)/matrix/testing; matrix_test


# C++ interface to the BLAS

$(BLASPP_LIB):
	cd $(LAPACKPP_DIR)/blaspp/src; make;

blas++_testing: 	
	cd $(LAPACKPP_DIR)/blaspp/testing; make;
	cd $(LAPACKPP_DIR)/blaspp/testing; blas++_test;



# make LAPACKPP_LIB
$(LAPACKPP_LIB) :
	cd $(LAPACKPP_DIR)/src; make


# (optionally) test the Lapack++ Driver Routines
#
# NOTE: this requires that the Lapack testing module 
# 	    (whose path is specified as $(TMGLIB)) be defined.
lapack++_testing :
	cd $(LAPACKPP_DIR)/testing; make
	cd $(LAPACKPP_DIR)/testing; lapack++_test

nolibs:
	cd $(LAPACKPP_LIB_DIR); rm *.a

clean: 
	cd $(LAPACKPP_DIR)/matrix/src; make clean;
	cd $(LAPACKPP_DIR)/matrix/testing; make clean;
	cd $(LAPACKPP_DIR)/blaspp/src; make clean;
	cd $(LAPACKPP_DIR)/blaspp/testing; make clean;
	cd $(LAPACKPP_DIR)/src; make clean;
	cd $(LAPACKPP_DIR)/testing; make clean;

# remove everything (libraries, executables) except input source

wipe: 
	cd $(LAPACKPP_DIR)/matrix/src; make wipe;
	cd $(LAPACKPP_DIR)/matrix/testing; make wipe;
	cd $(LAPACKPP_DIR)/blaspp/src; make wipe;
	cd $(LAPACKPP_DIR)/blaspp/testing; make wipe;
	cd $(LAPACKPP_DIR)/src; make wipe;
	cd $(LAPACKPP_DIR)/testing; make wipe;

