########################################################################
# Fast Auxiliary Space Preconditioners (FASP)
# 
# This is the Makefile for ASCPR test! 
#
########################################################################

#==============================================================================#
# User compilers                                                               #
# FASP has been tested with many different compilers; see README for details.  #
#==============================================================================#
CC  = icc 
CPP = icpc 
FC  = ifort -nofor_main
AR  = ar ruc

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-O2

#==============================================================================#
# OpenMP Support                                                               #
# Uncomment the following line to turn on the OpenMP support                   #
#==============================================================================#
BOPT+=-fopenmp

########################################################################
# Root directory for FASPSOLVER package
########################################################################
FASPDIR = ../faspsolver
INCLUDE = -I$(FASPDIR)/include -I./include -I$(FASPDIR)/base/src -I./src 
FASPLIB = $(FASPDIR)/lib/libfasp.a 


########################################################################
# External libraries for FASP package 
########################################################################
# We recommend using Pardiso.
# (1) If you wish to use Pardiso as a direct solver, uncomment the following
MKL_DIR = /opt/intel/mkl
MKL_LIB = -L $(MKL_DIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
MKL_INC = -I $(MKL_DIR)/include


# (2) If you wish to use MUMPS as a direct solver, uncomment the following
# MUMPS_DIR = [mumps_root_dir]
# MUMPS_INC = -I$(MUMPS_DIR)/libseq -I$(MUMPS_DIR)/include
# MUMPS_LIB = -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -lpord -L$(MUMPS_DIR)/libseq -lmpiseq -lblas


# (3) If you wish to use SuperLU as a direct solver, uncomment the following
# two lines and give the location of SuperLU head and lib:
# SUPERLU_INC = -I/usr/include/SuperLU_4.0
# SUPERLU_LIB = /usr/lib/libsuperlu_4.0.a


# (4) If you want to use UMFPACK as a direct solver, uncomment the following
# three lines as well as the BLASLIB
# UMFPACK_DIR=/dir/to/UMFPACK
# UMFPACK_INC=-I$(UMFPACK_DIR)/Include -I$(UMFPACK_DIR)/../AMD/Include -I$(UMFPACK_DIR)/../UFconfig 
# UMFPACK_LIB=-L$(UMFPACK_DIR)/../lib -lsuitesparseconfig -lumfpack -lamd -lcholmod -lcolamd -lcamd -lccolamd -L/opt/local/lib -lmetis


#==============================================================================#
# User preprocessing definitions                                               #
#==============================================================================#
CDEFS=
CDEFS+=-DWITH_PARDISO=1 -DWITH_MUMPS=0 -DWith_SuperLU=0 -DWith_UMFPACK=0 

COPTS=$(BOPT)
CINCLUDES=$(INCLUDE)
CINCLUDES+=$(MKL_INC) $(MUMPS_INC) $(SUPERLU_INC) $(UMFPACK_INC) 
CFLAGS=$(CDEFS) $(COPTS) $(CINCLUDES)

FOPTS=$(BOPT)
FDEFS=$(CDEFS)
FINCLUDES=$(CINCLUDES)
FFLAGS=$(FDEFS) $(FOPTS) $(FINCLUDES)

########################################################################
# Link options
########################################################################
LINKOPTS=$(BOPT)

LIBS=$(TESTLIB) $(FASPLIB)
LIBS+=$(MKL_LIB) $(MUMPS_LIB) $(SUPERLU_LIB) $(UMFPACK_LIB)

CLFLAGS=-lstdc++ $(LINKOPTS) $(LIBS)
FLFLAGS=-lm $(LINKOPTS) $(LIBS)

CSRCDIR = ./src
FSRCDIR = ./src
TESTLIB = ./lib/libfaspcpr.a

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .cc .cpp .for .f .f77 .f90 .f95
#
FSRC := $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.for))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f77))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f90))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f95))
CSRC := $(foreach dir,$(CSRCDIR),$(wildcard $(CSRCDIR)/*.c))
CSRC += $(foreach dir,$(EXTRDIR),$(wildcard $(EXTRDIR)/*.c))
#
OBJSF := $(patsubst %.for,%.o,$(FSRC))
OBJSF += $(patsubst %.f,%.o,$(FSRC))
OBJSF += $(patsubst %.f77,%.o,$(FSRC))
OBJSF += $(patsubst %.f90,%.o,$(FSRC))
OBJSF += $(patsubst %.f95,%.o,$(FSRC))
OBJSC := $(patsubst %.c,%.o,$(CSRC))
#
.for.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F object $@'
	@$(AR) $(TESTLIB) $@
#
.f90.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F90 object $@'
	@$(AR) $(TESTLIB) $@
#
.f95.o:
	@$(FC) -c $< -o $@ $(FFLAGS)
	@echo 'Building F95 object $@'
	@$(AR) $(TESTLIB) $@
#
.c.o:
	@$(CC) -c $< -o $@ $(CFLAGS)
	@echo 'Building C object $@'
	@$(AR) $(TESTLIB) $@
#
.cpp.o:
	@$(CPP) -c $< -o $@ $(CFLAGS)
	@echo 'Building CPP object $@'
	@$(AR) $(TESTLIB) $@
#
########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(TESTLIB)

########################################################################
# Link
########################################################################

all: $(ALLPROG) ascpr

Default:
	regression

headers: 
	cat $(CSRCDIR)/*.c \
	| awk -v name="faspcpr_functs.h" -f ./util/mkheaders.awk > ./include/faspcpr_functs.h

$(TESTLIB): $(OBJSC) $(OBJSF)
	@ranlib $(TESTLIB)
	@echo 'Generating library $@'

lib: $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)
	@echo 'Generating library $@'

########################################################################
# Some test problems
########################################################################

ascpr:
	@$(CC) $(CFLAGS) -c main/ascpr.c -o main/ascpr.o
	@$(FC) $(LOPT) main/ascpr.o $(FLFLAGS) -o test_ascpr.ex
	@echo 'Building executable $@'


########################################################################
# Clean up
########################################################################

.PHONY : clean distclean help

clean:
	@rm -f $(CSRCDIR)/*.o
	@rm -f $(FSRCDIR)/*.o
	@rm -f main/*.o

distclean:
	@make clean
	@rm -f lib/*.a
	@rm -f *~ *.ex *.out
	@rm -f $(CSRCDIR)/*~
	@rm -f $(FSRCDIR)/*~
	@rm -f main/*~

help:
	@echo "======================================================"
	@echo " Fast Auxiliary Space Preconditioners (FASP)"
	@echo "======================================================"
	@echo " "
	@echo " make            : build all exe files "
	@echo " make headers    : build the header file automatically"
	@echo " make clean      : clean all obj files "
	@echo " make distclean  : clean all obj, exe, bak, out files "
	@echo " make help       : show this screen "
	@echo " "
