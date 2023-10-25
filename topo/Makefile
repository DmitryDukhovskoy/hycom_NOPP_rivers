.SUFFIXES: .o .f90 .F90  # any suffix used in suffix rule should be declared

  FC      = /opt/intel/13/bin/ifort
  FFLAGS  = -convert big_endian -O3
#  IOLIBS = -L/opt/hpc/intel13/lib64 -lnetcdf -lnetcdff -limf -lm
  IOLIBS = -L/opt/hpc/intel13/lib64 -limf -lm
  IOINCS = -I/opt/hpc/intel13/include

# ---------------------
# LIBRARIES
# ---------------------
  LIBS  = $(IOLIBS)

  INCS  = $(IOINCS)

PRGRM  = topo_GLB2ARC.F90

MODS   = all_variables.F90 utils.F90 

SRCS   = $(MODS) $(PRGRM)

OBJS   = $(SRCS:.F90=.o) 

EXEC   = subTOPO.x

.F90.o:
	$(FC) -c $(FFLAGS) $(INCS) $<

$(EXEC):  $(OBJS)
	$(FC) $(OBJS) $(LIBS) -o $(EXEC)

#.F90.o:
#	$(FC) -c $(FFLAGS) $(INCS) $*.F90
#.F90.o:
#	$(FC) -c $(FFLAGS) $(INCS) $(SRCS)

all: $(EXEC) 

clean:
	/bin/rm -f *.o *.mod ${EXEC}

#all: clean make
