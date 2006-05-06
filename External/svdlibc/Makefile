CC=gcc
#CFLAGS=-Wall -O4
#CFLAGS=-g
#CC=icc
CFLAGS=-O3
LIBS=-lm

OBJ=svdlib.o svdutil.o las2.o

svd: Makefile main.o libsvd.a
	${CC} ${CFLAGS} -o svd main.o libsvd.a ${LIBS}
main.o: Makefile main.c svdlib.h
	${CC} ${CFLAGS} -c main.c

libsvd.a: ${OBJ}
	ar cr $@ ${OBJ}
	ranlib $@

svdlib.o: Makefile svdlib.h svdlib.c
	${CC} ${CFLAGS} -c svdlib.c
svdutil.o: Makefile svdutil.c svdutil.h
	${CC} ${CFLAGS} -c svdutil.c
las2.o: Makefile las2.c svdlib.h svdutil.h
	${CC} ${CFLAGS} -c las2.c
clean: 
	rm *.o
