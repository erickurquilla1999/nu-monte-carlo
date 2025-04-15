AMREX_HOME ?= ../amrex

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

TINY_PROFILE  = TRUE
USE_PARTICLES = TRUE

PRECISION     = DOUBLE

USE_MPI   = TRUE
USE_OMP   = TRUE
USE_CUDA  = TRUE
USE_HIP   = FALSE

###################################################

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules