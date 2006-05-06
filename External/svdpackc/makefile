#######################################################################
#                     SVDPACKC (Ver 1.0) Makefile                     #
#######################################################################
#                                                                     #
#      las1: Single-Vector Lanczos SVD via 2-Cyclic Eigensystems      #
#      las2: Single-Vector Lanczos SVD via   A'A    Eigensystems      #
#      bls1: Block         Lanczos SVD via 2-Cyclic Eigensystems      #
#      bls2: Block         Lanczos SVD via   A'A    Eigensystems      #
#      sis1: Subspace Iteration    SVD via 2_Cyclic Eigensystems      #
#      sis2: Subspace Iteration    SVD via   A'A    Eigensystems      #
#      tms1: Trace Minimization    SVD via 2_Cyclic Eigensystems      #
#      tms2: Trace Minimization    SVD via   A'A    Eigensystems      #
#                                                                     #
#######################################################################
 
CC    =	gcc
CFLAGS=	-O4 -c
#CFLAGS=	-g -c 
LIB=	-lm
TIMER=	timersun.o

all:	las1 las2 bls1 bls2 sis1 sis2 tms1 tms2

las1.o:	las1.h
las2.o:	las2.h

bls1.o:	bls1.h
bls2.o:	bls2.h

sis1.o:	sisg.h sisc.h
sis2.o:	sisg.h sisc.h
sis2-rd.o: sisg.h sisc.h

tms1.o:	tmsg.h tmsc.h
tms2.o:	tmsg.h tmsc.h

.c.o:   $*.c 
	${CC} ${CFLAGS} $*.c


las1:	las1.o ${TIMER}
	${CC} -o $@ las1.o ${TIMER} ${LIB}

las2:	las2.o ${TIMER}
	${CC} -o $@ las2.o ${TIMER} ${LIB}

bls1:	bls1.o ${TIMER}
	${CC} -o $@ bls1.o ${TIMER} ${LIB}

bls2:	bls2.o ${TIMER}
	${CC} -o $@ bls2.o ${TIMER} ${LIB}

sis1:	sis1.o ${TIMER}
	${CC} -o $@ sis1.o ${TIMER} ${LIB}

sis2:	sis2.o ${TIMER}
	${CC} -o $@ sis2.o ${TIMER} ${LIB}

sis2-rd: sis2-rd.o ${TIMER}
	${CC} -o $@ sis2-rd.o ${TIMER} ${LIB}

tms1:	tms1.o ${TIMER}
	${CC} -o $@ tms1.o ${TIMER} ${LIB}

tms2:	tms2.o ${TIMER}
	${CC} -o $@ tms2.o ${TIMER} ${LIB}

clean:
	las1.o las2.o bls1.o bls2.o sis1.o sis2.o tms1.o tms2.o \
	timersun.o timermac.o las1 las2 bls1 bls2 sis1 sis2 tms1 tms2
