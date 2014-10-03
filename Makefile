# Set the compiler flags
#FFLAGS = -fast#-fast #-arch pn4 -tpp7 -tune pn4 -ax
#FFLAGS = -r8 -fast -fastsse -O3 -Mipa=fast -Minline -Munroll	#PGI optimised
#FFLAGS = -Mbounds -g 				#PGI Debug
#FFLAGS = -O3 -fast                            	#Intel
#FFLAGS = -fpe0 -nothreads -traceback -fltconsistency -CB -g  #Intel Debug
 
FFLAGS = -O3

# Set some of the build parameters
TARGET = lare3d

#Uncomment the following line to use Qmono viscosity
#QMONO = -DQMONO

#Uncomment the following line to run in single precision
#QSINGLE = -DSINGLE

#Uncomment the following line to use fourth order RK scheme for resistive update
#QFOURTHORDER = -DFOURTHTORDER


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODULEFLAG = -module
OPFLAGS = $(QMONO)  $(QSINGLE) $(QFIRSTORDER)
FC = mpif90 $(OPFLAGS)
PREPROFLAGS = $(NONMPIIO)

OBJFILES = shared_data.o mpi_routines.o openboundary.o mpiboundary.o boundary.o normalise.o conduct.o diagnostics.o setup.o lagran.o  \
 remap.o xremap.o yremap.o zremap.o initial_conditions.o\
 output_cartesian.o iocontrol.o output.o iocommon.o input.o inputfunctions.o\
 input_cartesian.o neutral.o control.o\
 welcome.o lare3d.o
FULLTARGET = $(BINDIR)/$(TARGET)

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(OBJDIR):$(SRCDIR)/core:$(SRCDIR)/io/

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES))

clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ $(SRCDIR)/core/*~ $(SRCDIR)/io/*~ *.log

tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

touch:
	@touch src/* ; touch src/core/* 

datatidy:
	@rm -rf Data/*

visit:
	@cd VisIT; ./build

visitclean:
	@cd VisIT; make clean; \
	  rm -rf .depend *.d *Info.C *Info.h CMake* cmake* Makefile

.PHONY: clean tidy touch datatidy visit visitclean

# All the dependencies
shared_data.o:shared_data.F90
mpi_routines.o:mpi_routines.f90 shared_data.o 
normalise.o:normalise.f90 shared_data.o
setup.o:setup.F90 shared_data.o normalise.o iocommon.o iocontrol.o input.o input_cartesian.o
mpiboundary.o: mpiboundary.f90 shared_data.o
openboundary.o: openboundary.f90 shared_data.o
boundary.o:boundary.f90 shared_data.o mpiboundary.o 
xremap.o:xremap.f90 shared_data.o boundary.o
yremap.o:yremap.f90 shared_data.o boundary.o
zremap.o:zremap.f90 shared_data.o boundary.o
diagnostics.o:diagnostics.F90 shared_data.o boundary.o normalise.o output_cartesian.o output.o iocontrol.o 
iocommon.o:iocommon.f90 shared_data.o
output.o:output.f90 shared_data.o iocommon.o
output_cartesian.o: output_cartesian.f90 shared_data.o iocommon.o output.o
iocontrol.o: iocontrol.f90 shared_data.o iocommon.o output.o input.o
input.o: input.f90 shared_data.o iocommon.o inputfunctions.o
inputfunctions.o: inputfunctions.f90 shared_data.o iocommon.o  
input_cartesian.o: input_cartesian.f90 iocommon.o inputfunctions.o
conduct.o:conduct.f90 shared_data.o boundary.o
lagran.o:lagran.F90 shared_data.o boundary.o diagnostics.o neutral.o conduct.o
remap.o:remap.f90 shared_data.o xremap.o yremap.o zremap.o
initial_conditions.o:initial_conditions.f90 shared_data.o normalise.o neutral.o
neutral.o: neutral.f90 shared_data.o boundary.o normalise.o
control.o: control.f90 shared_data.o normalise.o
welcome.o: welcome.f90 shared_data.o
lare3d.o:lare3d.f90 shared_data.o setup.o boundary.o diagnostics.o lagran.o remap.o mpi_routines.o welcome.o initial_conditions.o openboundary.o control.o  
