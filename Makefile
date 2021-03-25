CC	= g++-9
SRC	= test.cpp 
CFLAGS	= -std=c++17 -O3 -march=native -mtune=native
LDFLAGS	= -lgmp -lmpfr -lntl -fopenmp

all:
	${CC} ${CFLAGS} ${SRC} ${LDFLAGS} -I/usr/include/eigen3


