#############################################################################
#
# Primitive makefile for MFDn_transitions
#
# Instructions:
#
#   - Uncomment appropriate compiler "configuration" below.
#
#   - To enable Tz-changing operators, uncomment appropriate line below.

$(warning This makefile is deprecated. Use 'cmake -B build' instead.)

################################################################
# common definitions
################################################################

# Recommended default settings for compiler preprocessor values
DEFMxN= -DMaxNumTBMEops=8 -DMaxNumKets=8 -DMAX_PATHLENGTH=1024
# Note: there is no MaxNumBras, there is only one bra state
#

# To enable Tz changing operators (|Delta_Tz| = 1 or 2), uncomment following line:
## DEFTz= -DDeltaTz

DEFS= $(DEFMxN) $(DEFTz) -DVCS_REVISION=\"old-makefile\"
#
# For intel optimization report
RPT= -list -qopt-report=1 -qopt-report-phase=loop,openmp,vec

################################################################
# configuration: NERSC Intel
################################################################

FCOMPILE= ftn
FFLAGS= $(DEFS) -O3 -ip -qopenmp $(RPT) -cpp

# debugging options
#FFLAGS= $(DEFS) -g -check all -warn -qopenmp -cpp
#FFLAGS= $(DEFS) -g -check bounds -warn -qopenmp -cpp

################################################################
# configuration: vanilla Intel
################################################################

FCOMPILE=mpiifort
FFLAGS= $(DEFS) -xHOST -Ofast -ip -qopenmp -cpp $(RPT) -no-wrap-margin

# debugging options
#FFLAGS= $(DEFS) -g -check bounds -warn -qopenmp -cpp -traceback -no-wrap-margin
#
#FCOMPILE= gfortran-mp-4.9
#FCOMPILE= gfortran-mp-8
#FFLAGS= $(DEFS) -fopenmp -cpp
#FFLAGS= $(DEFS) -fcheck=all -fopenmp -cpp

################################################################
# configuration: vanilla GCC
################################################################

## FCOMPILE= mpif90
## GENCMPLFLAGS= $(DEFS) -O3 -ffast-math -funroll-loops -fopenmp
## FFLAGS= $(GENCMPLFLAGS) $(TARGET_ARCH) -cpp
## F90FLAGS= $(GENCMPLFLAGS) $(TARGET_ARCH) -cpp

################################################################
# Intel MKL configuration -- from ndconfig
################################################################
# Note: MKLROOT not needed under Ubuntu native package.
MKL_CXXFLAGS := -m64 -I${MKLROOT}/include
MKL_LDFLAGS := -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed
MKL_LDLIBS := -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
################################################################
#
LDFLAGS+=$(MKL_LDFLAGS)
LDLIBS+=$(MKL_LDLIBS)

#############################################################################
#
# Either use module_TBME_Tz0
# or use module_TBME_Tz12 (for Tz changing operators, not yet tested)
#
objects := module_Transitions_input.o module_nodeinfo.o module_SPbasis.o \
	module_Sparsity.o module_Wig3J.o module_MjStates.o module_Wavefunctions.o \
	module_TBME_Tz0.o module_TBME_Tz12.o module_TBops.o module_OBops.o \
	subrts_general.o subrts_bitrep.o subrts_read_MBbasis.o subrts_MBnonzero.o \
	subrts_OBDME.o subrts_EvalOBDME.o subrts_inputTBME.o subrts_EvalTBTrans.o \
	subrts_ApplyTBTrans.o

################################################################
#
# build rules for object files
#
################################################################

COMPILE.f90 = $(FCOMPILE) $(FFLAGS) -c
LINK.f = $(FCOMPILE) $(FFLAGS) $(LDFLAGS)

#  command to compile (not link):
%.o: %.f90
		$(COMPILE.f90) -o $@ $<

################################################################
#
# build rule for mfdn executable
#
################################################################

ifdef DEFTz
xtransitions-DeltaTz: 	$(addprefix src/,$(objects)) MFDn_Transitions.o
		$(LINK.f) $^ $(LDLIBS) -o $@
default:
	xtransitions-DeltaTz
else
xtransitions: 	$(addprefix src/,$(objects)) MFDn_Transitions.o
		$(LINK.f) $^ $(LDLIBS) -o $@
default:
	xtransitions
endif

xapply: 	$(addprefix src/,$(objects)) MFDn_apply_operator.o
		$(LINK.f) $^ $(LDLIBS) -o $@

clean:
	rm -f src/*.o src/*.oo *.mod src/*__genmod.f90 *.lst src/*.optrpt xtransitions* xapply*

################################################################
