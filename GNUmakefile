AMREX_HOME ?= ../amrex

DEBUG	= FALSE
DEBUG	= TRUE

DIM	= 3

COMP    = gcc

USE_MPI   = TRUE
USE_OMP   = TRUE
USE_CUDA  = TRUE
USE_HIP   = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
