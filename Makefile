CC = gcc
CFLAGS = -std=c99 -g -Wall -MMD
OBJECTS = main.o matrix.o linear_program.o
DEPENDS = ${OBJECTS:.o=.d}
EXEC = lp-solve

${EXEC} : ${OBJECTS}
	${CC} ${OBJECTS} -o ${EXEC}

-include ${DEPENDS}

.PHONY : clean
clean:
	rm -rf ${DEPENDS} ${OBJECTS} ${EXEC}
