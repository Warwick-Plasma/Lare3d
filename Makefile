# Set the compiler flags
FFLAGS = -O3 

# Set some of the build parameters
TARGET = lare3d

# Set any precompiler options
#NONMPIIO = -DNONMPIIO   # Uncomment for non-MPI file access.

# Uncomment one of the following lines if on a cluster.
#COPSON = -compiler intel
#MHDCLUSTER = -f90=pgf90 -DMHDCLUSTER -fpic

#Uncomment the following line to use Qmono viscosity
#QMONO = -DQ_MONO


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODULEFLAG = -module
MACHINEFLAGS = $(COPSON) $(MHDCLUSTER)
OPFLAGS = $(QMONO)
FC = mpif90 $(MACHINEFLAGS) $(OPFLAGS)
PREPROFLAGS = $(NONMPIIO)

OBJFILES = shared_data.o mpi_routines.o openboundary.o mpiboundary.o boundary.o diagnostics.o  setup.o lagran.o  \
 remap.o xremap.o yremap.o zremap.o strings.o initial_conditions.o deck_control_block.o deck_boundaries_block.o deck.o welcome.o lare3d.o
FULLTARGET = $(BINDIR)/$(TARGET)

#vpath %.f90 $(SRCDIR)
#vpath %.o $(OBJDIR)
VPATH = $(SRCDIR):$(OBJDIR)

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES))

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: visit
visit:
	@cd VisIT;xml2makefile -clobber l3dv2.xml;make

	
# All the dependencies
shared_data.o:shared_data.F90
mpi_routines.o:mpi_routines.f90 shared_data.o 
setup.o:setup.F90 shared_data.o
openboundary.o:openboundary.f90 shared_data.o
mpiboundary.o:mpiboundary.f90 shared_data.o
boundary.o:boundary.f90 openboundary.o mpiboundary.o shared_data.o 
xremap.o:xremap.f90 shared_data.o boundary.o
yremap.o:yremap.f90 shared_data.o boundary.o
zremap.o:zremap.f90 shared_data.o boundary.o
diagnostics.o:diagnostics.F90 shared_data.o boundary.o
lagran.o:lagran.F90 shared_data.o boundary.o diagnostics.o
remap.o:remap.f90 shared_data.o xremap.o yremap.o zremap.o
strings.o:strings.f90 shared_data.o
initial_conditions.o:initial_conditions.f90 strings.o shared_data.o
deck_control_block.o:deck_control_block.f90 strings.o shared_data.o
deck_boundaries_block.o:deck_boundaries_block.f90 strings.o shared_data.o
deck.o:deck.f90 strings.o initial_conditions.o shared_data.o deck_control_block.o deck_boundaries_block.o
welcome.o:welcome.f90 shared_data.o
lare3d.o:lare3d.f90 shared_data.o setup.o boundary.o diagnostics.o lagran.o remap.o mpi_routines.o deck.o initial_conditions.o welcome.o
