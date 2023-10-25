.SUFFIXES: .o .f90 .F90  # any suffix used in suffix rule should be declared

  FC      = /opt/intel/12/bin/ifort
  FFLAGS  = -convert big_endian -O3
#  IOLIBS = -L/opt/hpc/intel13/lib64 -lnetcdf -lnetcdff -limf -lm
  IOLIBS = -L/opt/hpc/intel13/lib64 -limf -lm
  IOINCS = -I/opt/hpc/intel13/include

# ---------------------
# LIBRARIES
# ---------------------
  LIBS  = $(IOLIBS)

  INCS  = $(IOINCS)

PRGRM  = Fill_OBnans.F90

MODS   = all_variables.F90 

SRCS   = $(MODS) $(PRGRM)

OBJS   = $(SRCS:.F90=.o) 

EXEC   = fillnans.x

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
