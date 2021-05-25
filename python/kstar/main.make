# Makefile for Experimental Equilibrium Rebuild
SHELL=/bin/csh

NCLIB = -lblas -llapack
SPECIAL   = -L/usr/lib64

F90=gfortran

F90_SRCS = input_params.f90 set_boundary.f90 f_solver.f90 \
           critical.f90 main.f90 toms526.f90

F_SRCS = 

OBJS =  ${F90_SRCS:.f90=.o} ${F_SRCS:.f=.o}
 
OPT =

LIBS = # -L. -lssl2mt # -L./TOMS526 -ltoms5

all: ${OBJS}
	${F90} $(SPECIAL) -o main ${OBJS} ${LIBS} $(NCLIB)

main.o: main.f90 input_params.o set_boundary.o
	${F90} -c ${OPT} main.f90 

critical.o: critical.f90 input_params.o
	${F90} -c  ${OPT} critical.f90

set_boundary.o: set_boundary.f90 
	${F90} -c  ${OPT} set_boundary.f90

f_solver.o: f_solver.f90 input_params.o
	${F90} -c  ${OPT} f_solver.f90

toms526.o: toms526.f90 # input_params.o
	${F90} -c  ${OPT} toms526.f90

input_params.o: input_params.f90
	${F90} -c $?

clean:
	rm -f main core.* *.o *.mod 
