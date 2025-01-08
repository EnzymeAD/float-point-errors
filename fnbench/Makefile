.SUFFIXES: .c .o

CC = clang
CFLAGS = -O3 -Wall -DREORDER
LFLAGS = -lm

TSRC  = src/meshOpt.c src/fcn3e_am.c src/elem3.c src/tet.c src/mesh3.c src/mesh.c src/opt3.c src/pre3.c 

tetOpt: ${TSRC}
	${CC} ${CFLAGS} -Iinclude -o tetOpt ${TSRC} ${LFLAGS}

clean:
	/bin/rm -f tetOpt
	/bin/rm -f core core.* *.o *.il *~
