# Makefile for Experimental Equilibrium Rebuild
SHELL=/bin/csh

NCLIB = -lblas -llapack
SPECIAL   = -L/usr/lib64

F90=gfortran
F2PY=f2py

#F90_MAIN_SRCS = main.f90 input_params.f90

F90_SRCS = set_boundary.f90 f_solver.f90 \
           critical.f90 toms526.f90

F2PY_SRCS = pyfeq.f90 pyfeq_mod.f90


OBJS =  ${F90_SRCS:.f90=.o} ${F_SRCS:.f=.o}
 
OPT = -fPIC

LIBS = # -L. -lssl2mt # -L./TOMS526 -ltoms5

all: ${OBJS}
#	${F90} $(SPECIAL) -o main ${MAIN_OBJS} ${OBJS} ${LIBS} $(NCLIB)
	${F2PY} -c $(SPECIAL) $(NCLIB) -m pyfeq ${F2PY_SRCS} ${OBJS} ${LIBS} 

#main.o: main.f90 pyfeq_mod.o
#	${F90} -c  ${OPT} main.f90

critical.o: critical.f90 pyfeq_mod.o
	${F90} -c  ${OPT} critical.f90

set_boundary.o: set_boundary.f90 pyfeq_mod.o
	${F90} -c  ${OPT} set_boundary.f90

f_solver.o: f_solver.f90 pyfeq_mod.o
	${F90} -c  ${OPT} f_solver.f90

toms526.o: toms526.f90 # input_params.o
	${F90} -c  ${OPT} toms526.f90

pyfeq_mod.o: pyfeq_mod.f90
	${F90} -c $?

clean:
	rm -f core.* *.o *.mod
